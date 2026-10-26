#include "geometry.h"
#include <cassert>
int main(){
 cell::Geometry g(5e-6,.4,50e-6,true,true);
 assert(std::abs(cell::width-.0043)<1e-12);
 assert(std::abs(g.latticeHeight-2*.001835)<1e-12);
 const int nx=860,ny=734;
 for(int j=2;j<ny-1;++j){assert(g.material(5e-6,j*5e-6)==2);assert(g.material((nx-1)*5e-6,j*5e-6)==2);}
 for(int i=2;i<nx-1;++i){assert(g.material(i*5e-6,5e-6)==2);assert(g.material(i*5e-6,(ny-1)*5e-6)==7);}
 assert(g.material(.00215,.00025)==3);
}
