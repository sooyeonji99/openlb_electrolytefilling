#include "geometry.h"
#define cell original_cell
#include "../NMC_longterm_case/geometry.h"
#undef cell
#include <cassert>
#include <iostream>
int main(){
 cell::Geometry reflected(5e-6);original_cell::Geometry original(5e-6);
 int ny=std::lround(reflected.latticeHeight/reflected.dx),nx=860;
 for(int j=0;j<=ny;++j)for(int i=0;i<=nx;++i){
  int expected=original.material(i*5e-6,(ny-j)*5e-6);
  if(expected==7)expected=6;
  assert(reflected.material(i*5e-6,j*5e-6)==expected);
 }
 assert(reflected.material(cell::width/2,reflected.latticeHeight-20e-6)==3);
 assert(reflected.material(cell::width/2,reflected.latticeHeight-5e-6)==2);
 assert(reflected.material(cell::width/2,5e-6)==6);
 assert(reflected.material(5e-6,reflected.latticeHeight/2)==4);
 assert(reflected.material(cell::width-5e-6,reflected.latticeHeight/2)==5);
 reflected.write("results/geometry_validation");
 std::cout<<"All lattice cells match vertical reflection; top source and wall, side and bottom outlets PASS\n";
}
