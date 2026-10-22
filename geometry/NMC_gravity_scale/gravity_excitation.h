#pragma once
#include <cmath>
#include <stdexcept>
struct GravityExcitation {
  static constexpr double referenceG=9.81;
  double frequency=100,scale=1;
  void validate()const {
    if(!std::isfinite(frequency)||frequency<=0||!std::isfinite(scale)||scale<0||!std::isfinite(peak()))
      throw std::runtime_error("frequency must be positive and gravity-scale must be nonnegative and finite");
  }
  double peak()const{return scale*referenceG;}
  double phase(double time)const{return 2*std::acos(-1.)*frequency*time;}
  double at(double time)const{return peak()*std::sin(phase(time));}
};
