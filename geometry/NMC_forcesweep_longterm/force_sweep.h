#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

struct SweepState {
  double frequency{}, frequencyRate{}, phase{}, inertial{};
  const char* segment="pre_rest";
  bool active=false;
};
struct ForceSweep {
  double pre=.03,up=.5,down=.5,post=.03,low=10,high=100;
  double duration()const{return pre+up+down+post;}
  void validate()const {
    for(double v:{pre,up,down,post,low,high})if(!std::isfinite(v))throw std::runtime_error("nonfinite sweep parameter");
    if(pre<0||post<0||up<=0||down<=0||low<=0||high<low||!std::isfinite(duration()))throw std::runtime_error("invalid sweep intervals or frequencies");
  }
  SweepState at(double time,double accelerationAmplitude)const {
    SweepState s;
    const double eps=32*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(time),duration()});
    if(time<pre-eps)return s;
    if(time>=pre+up+down-eps){s.segment="post_rest";return s;}
    double t=std::max(0.,time-pre),cycles;
    if(t<up-eps){
      s.segment="sweep_up";s.frequencyRate=(high-low)/up;
      s.frequency=low+s.frequencyRate*t;
      cycles=low*t+.5*s.frequencyRate*t*t;
    } else {
      s.segment="sweep_down";t=std::max(0.,t-up);
      s.frequencyRate=(low-high)/down;
      s.frequency=high+s.frequencyRate*t;
      cycles=.5*(low+high)*up+high*t+.5*s.frequencyRate*t*t;
    }
    const double pi=std::acos(-1.);
    s.active=true;s.phase=2*pi*cycles;
    // Acceleration-controlled chirp: fixed peak g level, variable frequency.
    // No displacement trajectory or frequency-rate correction is imposed.
    s.inertial=accelerationAmplitude*std::sin(s.phase);
    return s;
  }
};
