#include "gravity_excitation.h"
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
int main() {
  GravityExcitation e;
  e.validate();
  const double dt=(.5075-.5)*5e-6*5e-6/(3e-6);
  // Independent integer-step oracle: 80,000 on + 16,000 rest, 20 sets.
  long activeSteps=0;
  const double pi=std::acos(-1.);
  double integratedPhase=0;
  for(long i=0;i<=1920000;++i) {
    const double t=i*dt;
    const bool expected=i<1920000 && i%96000<80000;
    assert(e.active(t)==expected);
    activeSteps+=e.active(t);
    if(i) {
      const double midpoint=(i-.5)*dt;
      const double f=55-45*std::cos(2*pi*midpoint/.12);
      integratedPhase+=2*pi*f*dt;
    }
    assert(std::abs(e.phase(t)-integratedPhase)<1e-8);
    const double expectedAcceleration=expected?9.81*std::sin(integratedPhase):0;
    assert(std::abs(e.at(t)-expectedAcceleration)<1e-7);
    assert(e.frequencyAt(t)>=10-1e-12 && e.frequencyAt(t)<=100+1e-12);
    if(!expected)assert(e.at(t)==0);
  }
  assert(activeSteps==1600000);
  assert(std::abs(e.activeTime(.12)-.100)<1e-14);
  assert(std::abs(e.activeTime(.0055)-.005)<1e-14);
  assert(e.active(.004999) && !e.active(.005));
  assert(!e.active(.005999) && e.active(.006));
  assert(!e.active(.119) && !e.active(.12));
  assert(std::abs(e.frequencyAt(0)-10)<1e-12);
  assert(std::abs(e.frequencyAt(.03)-55)<1e-12);
  assert(std::abs(e.frequencyAt(.06)-100)<1e-12);
  assert(std::abs(e.frequencyAt(.12)-10)<1e-12);
  assert(std::abs(e.phase(.12)/(2*pi)-6.6)<1e-12);
  for(double t:{.001,.006,.03,.06,.09,.119}) {
    const double h=1e-7;
    double derivative=(e.phase(t+h)-e.phase(t-h))/(2*h)/(2*pi);
    assert(std::abs(derivative-e.frequencyAt(t))<1e-7);
  }
  const double unit=e.at(.0025);
  for(double scale:{0.,.5,1.,2.}) {
    e.scale=scale;e.validate();
    assert(std::abs(e.at(.0025)-scale*unit)<1e-12);
    assert(e.at(.0055)==0);
  }
  e.scale=1;e.restTime=0;e.minFrequency=100;e.validate();
  assert(e.active(.005) && e.active(.006));
  assert(std::abs(e.at(.0075)+9.81)<1e-12);
  for(int invalid=0;invalid<7;++invalid) {
    GravityExcitation bad;
    if(invalid==0)bad.scale=-1;
    if(invalid==1)bad.onTime=0;
    if(invalid==2)bad.restTime=-1;
    if(invalid==3)bad.restTime=std::numeric_limits<double>::infinity();
    if(invalid==4)bad.minFrequency=101;
    if(invalid==5)bad.minFrequency=0;
    if(invalid==6)bad.sweepPeriod=0;
    bool caught=false;try{bad.validate();}catch(const std::runtime_error&){caught=true;}
    assert(caught);
  }
  e=GravityExcitation{};
  std::ofstream f("results/forcing_schedule_1g.csv");
  f<<std::setprecision(17)<<"time_s,frequency_hz,phase_rad,excitation_active,excitation_acceleration_m_s2,total_acceleration_y_m_s2,pressure_pa\n";
  for(int i=0;i<=12000;++i) {
    double t=i*1e-5;
    f<<t<<','<<e.frequencyAt(t)<<','<<e.phase(t)<<','<<e.active(t)<<','<<e.at(t)<<','<<-9.81+e.at(t)<<",2500\n";
  }
  std::cout<<"PASS: 1,920,000 step schedule; 20 sets; 100 ms on / 20 ms rest; 10-100-10 Hz sweep; integrated phase; invalid options\n";
}
