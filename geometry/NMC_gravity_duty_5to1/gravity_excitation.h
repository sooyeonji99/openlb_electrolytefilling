#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
struct GravityExcitation {
  static constexpr double referenceG=9.81;
  double frequency=100,scale=1;
  double minFrequency=10,sweepPeriod=.12;
  double onTime=.005,restTime=.001,endTime=.12;
  void validate()const {
    if(!std::isfinite(frequency)||frequency<=0||!std::isfinite(scale)||scale<0||!std::isfinite(peak()))
      throw std::runtime_error("frequency must be positive and gravity-scale must be nonnegative and finite");
    if(!std::isfinite(minFrequency)||minFrequency<=0||minFrequency>frequency||
       !std::isfinite(sweepPeriod)||sweepPeriod<=0)
      throw std::runtime_error("frequency bounds must satisfy 0 < min <= max; sweep period must be positive and finite");
    if(!std::isfinite(onTime)||onTime<=0||!std::isfinite(restTime)||restTime<0||!std::isfinite(period())||
       !std::isfinite(endTime)||endTime<0)
      throw std::runtime_error("excitation-on-s must be positive; rest and end times must be nonnegative and finite");
  }
  double peak()const{return scale*referenceG;}
  double period()const{return onTime+restTime;}
  double tolerance(double time)const {
    return 64*std::numeric_limits<double>::epsilon()*std::max(period(),std::abs(time));
  }
  // Snap boundaries to avoid floating point drift at integer lattice steps.
  double localTime(double time)const {
    double local=std::fmod(time,period());
    if(local<tolerance(time)||period()-local<tolerance(time))local=0;
    return local;
  }
  bool active(double time)const {
    if(!std::isfinite(time)||time<0||time>=endTime-tolerance(time))return false;
    return restTime==0||localTime(time)<onTime-tolerance(time);
  }
  double activeTime(double time)const {
    time=std::clamp(time,0.,endTime);
    double sets=std::floor((time+tolerance(time))/period());
    return sets*onTime+std::clamp(time-sets*period(),0.,onTime);
  }
  double frequencyAt(double time)const {
    const double angle=2*std::acos(-1.)*time/sweepPeriod;
    return minFrequency+(frequency-minFrequency)*(1-std::cos(angle))/2;
  }
  // Integrate instantaneous frequency: using 2*pi*f(t)*t would add t*f'(t).
  double phase(double time)const {
    const double middle=(frequency+minFrequency)/2;
    const double halfRange=(frequency-minFrequency)/2;
    return 2*std::acos(-1.)*middle*time-halfRange*sweepPeriod*std::sin(2*std::acos(-1.)*time/sweepPeriod);
  }
  double at(double time)const{return active(time)?peak()*std::sin(phase(time)):0.;}
};
