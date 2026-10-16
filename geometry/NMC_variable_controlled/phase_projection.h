#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

// Interface-weighted conservative projection, not a phase-field discretization.
// sum() must be a global all-reduce for distributed lattices. Pure phases are
// not used as redistribution capacity: if no diffuse interface remains, fail.
struct PhaseProjectionStats {
  double massBefore{},massAfter{},l1{},changed{},rawMin{},rawMax{};
  int iterations{};
};
template<class Sum>
PhaseProjectionStats projectPhase(std::vector<double>& values,Sum sum) {
  PhaseProjectionStats s;
  s.rawMin=std::numeric_limits<double>::infinity();s.rawMax=-s.rawMin;
  double invalid=0,count=double(values.size());
  std::vector<double> original=values;
  for(double v:values){invalid+=!std::isfinite(v);s.massBefore+=v;s.rawMin=std::min(s.rawMin,v);s.rawMax=std::max(s.rawMax,v);}
  if(sum(invalid)>0)throw std::runtime_error("Nonfinite phase before conservative projection");
  s.massBefore=sum(s.massBefore);count=sum(count);
  const double tol=1e-12*std::max(1.,std::abs(s.massBefore));
  if(s.massBefore < -tol || s.massBefore > count+tol)throw std::runtime_error("Phase mass incompatible with bounded conservative projection");
  double boundedMass=0;
  for(double& v:values){v=std::clamp(v,0.,1.);boundedMass+=v;}
  double remaining=s.massBefore-sum(boundedMass);
  while(std::abs(remaining)>tol){
    if(++s.iterations>128)throw std::runtime_error("Conservative phase projection did not converge");
    double weight=0;for(double v:values)weight+=v*(1-v);
    weight=sum(weight);
    if(weight<=0)throw std::runtime_error("No diffuse-interface capacity for conservative phase projection");
    const double factor=remaining/weight;
    double change=0;
    for(double& v:values){double next=std::clamp(v+factor*v*(1-v),0.,1.);change+=next-v;v=next;}
    remaining-=sum(change);
  }
  for(size_t i=0;i<values.size();++i){s.massAfter+=values[i];s.l1+=std::abs(values[i]-original[i]);s.changed+=std::abs(values[i]-original[i])>1e-12;}
  s.massAfter=sum(s.massAfter);s.l1=sum(s.l1);s.changed=sum(s.changed);
  if(std::abs(s.massAfter-s.massBefore)>4*tol)throw std::runtime_error("Conservative phase projection mass residual exceeds tolerance");
  return s;
}
