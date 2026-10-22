#include "gravity_excitation.h"
#include <cassert>
#include <fstream>
#include <iomanip>
int main(){
 GravityExcitation e;e.validate();const double pi=std::acos(-1.);
 assert(std::abs(e.phase(.12)/(2*pi)-12)<1e-12);
 assert(std::abs(e.at(0))<1e-12);
 assert(std::abs(e.at(.0025)-9.81)<1e-12);
 assert(std::abs(e.at(.0075)+9.81)<1e-12);
 for(double scale:{0.,.5,1.,2.}){
  e.scale=scale;e.validate();
  for(int i=0;i<=12000;++i){double t=i*1e-5;assert(std::abs(e.at(t))<=scale*9.81+1e-12);}
  assert(std::abs(e.at(.0025)-scale*9.81)<1e-12);
 }
 bool bad=false;try{e.scale=-1;e.validate();}catch(...){bad=true;}assert(bad);
 e.scale=1;
 std::ofstream f("results/forcing_schedule_1g.csv");
 f<<std::setprecision(17)<<"time_s,frequency_hz,gravity_scale,excitation_acceleration_m_s2,total_acceleration_y_m_s2,pressure_pa\n";
 for(int i=0;i<=12000;++i){double t=i*1e-5;f<<t<<",100,1,"<<e.at(t)<<','<<-9.81+e.at(t)<<",2500\n";}
}
