#include "force_sweep.h"
#include <cassert>
#include <fstream>
#include <iomanip>
int main(){
 ForceSweep q;q.validate();const double accel=9.81,pi=std::acos(-1.);
 assert(std::abs(q.duration()-1.06)<1e-14);
 assert(!q.at(0,accel).active&&!q.at(.029,accel).active);
 assert(q.at(.03,accel).active&&std::abs(q.at(.03,accel).frequency-10)<1e-12);
 assert(std::abs(q.at(.28,accel).frequency-55)<1e-12);
 assert(std::abs(q.at(.53,accel).frequency-100)<1e-12);
 assert(std::abs(q.at(.78,accel).frequency-55)<1e-12);
 assert(!q.at(1.03,accel).active&&!q.at(1.06,accel).active);
 assert(std::abs(q.at(.53,accel).phase/(2*pi)-27.5)<1e-12);
 assert(std::abs(q.at(1.03-1e-9,accel).phase/(2*pi)-55)<1e-6);
 for(double t:{.08,.28,.5,.58,.78,1.0}){
  const double h=1e-7;
  auto s=q.at(t,accel),a=q.at(t-h,accel),b=q.at(t+h,accel);
  assert(std::abs((b.phase-a.phase)/(2*h)-2*pi*s.frequency)<1e-5);
  assert(std::abs(s.inertial-accel*std::sin(s.phase))<1e-12);
  assert(std::abs(q.at(t,2*accel).inertial-2*s.inertial)<1e-12);
  assert(q.at(t,0).inertial==0);
 }
 double peak=0,trough=0;
 for(int i=0;i<=106000;++i){auto s=q.at(i*1e-5,accel);assert(std::abs(s.inertial)<=accel+1e-12);peak=std::max(peak,s.inertial);trough=std::min(trough,s.inertial);}
 assert(peak>9.8099&&trough<-9.8099);
 assert(std::abs(q.at(.53-1e-9,accel).phase-q.at(.53+1e-9,accel).phase)<2e-6);
 bool bad=false;try{auto s=q;s.up=0;s.validate();}catch(...){bad=true;}assert(bad);
 std::ofstream out("results/forcing_schedule_full.csv");
 out<<std::setprecision(17)<<"time_s,segment,frequency_hz,frequency_rate_hz_s,phase_rad,excitation_g,inertial_acceleration_m_s2,effective_acceleration_m_s2,pressure_pa\n";
 for(int i=0;i<=10600;++i){double t=i*.0001;auto s=q.at(t,accel);out<<t<<','<<s.segment<<','<<s.frequency<<','<<s.frequencyRate<<','<<s.phase<<",1,"<<s.inertial<<','<<-9.81+s.inertial<<",2500\n";}
}
