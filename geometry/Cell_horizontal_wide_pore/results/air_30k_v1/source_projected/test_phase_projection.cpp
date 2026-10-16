#include "phase_projection.h"
#include <cassert>
#include <random>
int main(){
  auto sum=[](double v){return v;};
  for(auto input:{std::vector<double>{-.02,.2,.3,.8,1.01},std::vector<double>{-.01,.2,.7,1.08},std::vector<double>{0,.2,.8,1}}){
    auto expected=input;auto stats=projectPhase(input,sum);
    double a=0,b=0;for(double v:expected)a+=v;for(double v:input){assert(v>=0&&v<=1);b+=v;}
    assert(std::abs(a-b)<1e-10);assert(std::abs(stats.massBefore-stats.massAfter)<1e-10);
  }
  std::mt19937 gen(17);std::uniform_real_distribution<double> dist(-.1,1.1);
  for(int trial=0;trial<100;++trial){std::vector<double> v(1000);for(double& x:v)x=dist(gen);auto s=projectPhase(v,sum);for(double x:v)assert(x>=0&&x<=1);assert(std::abs(s.massAfter-s.massBefore)<1e-8);}
  bool rejected=false;try{std::vector<double> v{-.2,-.1};projectPhase(v,sum);}catch(const std::runtime_error&){rejected=true;}assert(rejected);
}
