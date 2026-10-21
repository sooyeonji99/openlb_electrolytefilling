#pragma once
#include <cmath>
#include <limits>
// The force window is [start,end); each burst resets the sine phase.
// on=off=0 preserves continuous excitation within the window.
struct ExcitationSchedule {
  double start,duration,on,off;
  double phaseTime(double time) const {
    double elapsed=time-start;
    const double eps=32*std::numeric_limits<double>::epsilon()*
      std::fmax(1.,std::fmax(std::abs(time),std::abs(start+duration)));
    if(elapsed < -eps || elapsed >= duration-eps)return -1.;
    if(elapsed<0)elapsed=0;
    if(on==0)return elapsed;
    const double period=on+off;
    double phase=std::fmod(elapsed,period);
    if(period-phase<=eps)phase=0;
    if(phase>=on-eps)return -1.;
    return phase;
  }
};
