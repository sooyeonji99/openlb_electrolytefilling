#include "force_sweep.h"
#include <cassert>
#include <fstream>
#include <iomanip>
#include <vector>
int main(){
 ForceSweep q;q.validate();const double A=.00734,pi=std::acos(-1.);
 assert(std::abs(q.duration()-1.06)<1e-14);
 assert(!q.at(0,A).active);assert(!q.at(.029,A).active);
 assert(q.at(.03,A).active&&std::abs(q.at(.03,A).frequency-10)<1e-12);
 assert(std::abs(q.at(.28,A).frequency-55)<1e-12);
 assert(std::abs(q.at(.53,A).frequency-100)<1e-12);
 assert(std::abs(q.at(.78,A).frequency-55)<1e-12);
 assert(q.at(1.029999,A).frequency>10);
 assert(!q.at(1.03,A).active&&!q.at(1.06,A).active);
 assert(std::abs(q.at(.53,A).phase/(2*pi)-27.5)<1e-12);
 assert(std::abs(q.at(1.03-1e-9,A).phase/(2*pi)-55)<1e-6);
 for(double t:{.08,.28,.5,.58,.78,1.0}){
  const double h=1e-7;
  auto s=q.at(t,A),a=q.at(t-h,A),b=q.at(t+h,A);
  const double dphase=(b.phase-a.phase)/(2*h);
  const double accel=(b.velocity-a.velocity)/(2*h);
  assert(std::abs(dphase-2*pi*s.frequency)<1e-5);
  assert(std::abs(-accel-s.inertial)<1e-3);
 }
 auto left=q.at(.53-1e-9,A),right=q.at(.53+1e-9,A);
 assert(std::abs(left.phase-right.phase)<2e-6);
 bool rejected=false;try{auto bad=q;bad.up=0;bad.validate();}catch(...){rejected=true;}assert(rejected);
 std::ofstream out("results/forcing_schedule_full.csv");
 out<<std::setprecision(17)<<"time_s,segment,frequency_hz,frequency_rate_hz_s,phase_rad,displacement_m,velocity_m_s,inertial_acceleration_m_s2,effective_acceleration_m_s2,pressure_pa\n";
 for(int i=0;i<=1060;++i){double t=i*.001;auto s=q.at(t,A);out<<t<<','<<s.segment<<','<<s.frequency<<','<<s.frequencyRate<<','<<s.phase<<','<<s.displacement<<','<<s.velocity<<','<<s.inertial<<','<<-9.81+s.inertial<<",2500\n";}
}
