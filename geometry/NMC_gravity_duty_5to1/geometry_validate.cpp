#include "geometry.h"
#include <iostream>
int main(int argc,char** argv){try{double dx=argc>1?std::stod(argv[1])*1e-6:5e-6;cell::Geometry g(dx,.4,50e-6,true,true);g.write(argc>2?argv[2]:"results/geometry_dx05");std::cout<<"Geometry PASS\n";}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
