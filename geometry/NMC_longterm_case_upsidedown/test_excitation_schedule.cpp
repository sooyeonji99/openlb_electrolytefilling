#include "excitation_schedule.h"
#include <cassert>
#include <cmath>
int main(){
 ExcitationSchedule s{0,.4,.01,.01};
 for(int n=0;n<20;++n){
  const double t=.02*n;
  assert(std::abs(s.phaseTime(t))<1e-12);
  assert(std::abs(s.phaseTime(t+.0025)-.0025)<1e-12);
  assert(std::abs(s.phaseTime(t+.0075)-.0075)<1e-12);
  assert(s.phaseTime(t+.01)<0);
  assert(s.phaseTime(t+.015)<0);
 }
 assert(s.phaseTime(.4)<0);assert(s.phaseTime(-.001)<0);
 ExcitationSchedule delayed{.2,.2,.01,.01};
 assert(delayed.phaseTime(.1)<0);assert(delayed.phaseTime(.21)<0);
 assert(std::abs(delayed.phaseTime(.22))<1e-12);
 ExcitationSchedule continuous{0,.4,0,0};
 assert(std::abs(continuous.phaseTime(.015)-.015)<1e-12);
 assert(continuous.phaseTime(.4)<0);
}
