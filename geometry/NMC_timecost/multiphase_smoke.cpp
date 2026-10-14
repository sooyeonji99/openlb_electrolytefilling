#include "olb2D.h"
#include "olb2D.hh"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::names;
using T = double;
using NSD = D2Q9<RHO,NABLARHO,FORCE,EXTERNAL_FORCE,TAU_EFF,STATISTIC,SCALAR>;
using CHD = D2Q9<CONV_POPS,FORCE,SOURCE,VELOCITY,OLD_PHIU,STATISTIC,THETA,BOUNDARY>;
using Op = LiangPostProcessor;
using Maps = meta::map<NavierStokes,descriptors::VALUED_DESCRIPTOR<T,NSD>,
                       Component1,descriptors::VALUED_DESCRIPTOR<T,CHD>>;
using Coupling = SuperLatticeCoupling<Op,Maps>;

struct Options {
  T radius=4.3e-3, electrolyteRadius=.8925e-3, dx=5e-6;
  T porosity=.30, pitch=50e-6, tau=.8, nu=1e-6, rhoL=1000, rhoG=100;
  T sigma=1e-5, theta=60, width=4;
  long steps=200, outputEvery=200, diagnosticEvery=10;
  std::string output="results/multiphase_smoke";
};

Options parse(int argc,char**argv) {
  Options o;
  for (int i=1;i<argc;++i) {
    if (i+1>=argc) throw std::runtime_error("missing option value");
    std::string key=argv[i],value=argv[++i];
    if(key=="--dx-um")o.dx=std::stod(value)*1e-6;
    else if(key=="--steps")o.steps=std::stol(value);
    else if(key=="--output-every")o.outputEvery=std::stol(value);
    else if(key=="--diagnostic-every")o.diagnosticEvery=std::stol(value);
    else if(key=="--output")o.output=value;
    else throw std::runtime_error("unknown option: "+key);
  }
  return o;
}

class RadialPhase : public AnalyticalF2D<T,T> {
  T _r,_w;
public:
  RadialPhase(T r,T w):AnalyticalF2D<T,T>(1),_r(r),_w(w){}
  bool operator()(T out[],const T p[]) override {
    out[0]=.5*(1-std::tanh(2*(std::hypot(p[0],p[1])-_r)/_w));
    return true;
  }
};

class Phi : public SuperF2D<T,T> {
  SuperLattice<T,CHD>& _lattice;
  SuperGeometry<T,2>& _geometry;
public:
  Phi(SuperLattice<T,CHD>& lattice, SuperGeometry<T,2>& geometry)
    :SuperF2D<T,T>(lattice,1),_lattice(lattice),_geometry(geometry){this->getName()="phi_fluid_only";}
  bool operator()(T out[],const int p[]) override {
    auto& load=_lattice.getLoadBalancer();
    if(!load.isLocal(p[0]))return false;
    if(_geometry.getBlockGeometry(load.loc(p[0])).get(p[1],p[2])!=1) {
      out[0]=std::numeric_limits<T>::quiet_NaN();
      return true;
    }
    out[0]=_lattice.getBlock(load.loc(p[0])).get(p[1],p[2]).template getFieldComponent<STATISTIC>(0);
    return true;
  }
};

struct Sample { T mean{},minimum{1e9},maximum{-1e9},cells{}; };
Sample sample(SuperGeometry<T,2>& geometry,SuperLattice<T,CHD>& ch) {
  Sample s;
  for (int iC=0; iC<ch.getLoadBalancer().size(); ++iC) {
    auto& g=geometry.getBlockGeometry(iC); auto& lattice=ch.getBlock(iC);
    g.forSpatialLocations([&](LatticeR<2> r){
      if(g.get(r)!=1||lattice.isPadding(r))return;
      T phi=lattice.get(r).template getFieldComponent<STATISTIC>(0);
      s.mean+=phi;s.minimum=std::min(s.minimum,phi);s.maximum=std::max(s.maximum,phi);++s.cells;
    });
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(s.mean,MPI_SUM);
  singleton::mpi().reduceAndBcast(s.cells,MPI_SUM);
  singleton::mpi().reduceAndBcast(s.minimum,MPI_MIN);
  singleton::mpi().reduceAndBcast(s.maximum,MPI_MAX);
#endif
  if(s.cells)s.mean/=s.cells;
  return s;
}

int main(int argc,char**argv) {
  initialize(&argc,&argv); const Options o=parse(argc,argv);
  std::filesystem::create_directories(o.output);
  singleton::directories().setOutputDir(o.output+"/");
  IndicatorCuboid2D<T> domain({o.radius,o.radius},{0,0});
  CuboidDecomposition2D<T> cuboids(domain, o.dx, singleton::mpi().getSize());
  HeuristicLoadBalancer<T> load(cuboids);
  SuperGeometry<T,2> geometry(cuboids,load,2); geometry.rename(0,1);
  const T particleRadius=o.pitch*std::sqrt((1-o.porosity)/M_PI);
  auto& block=geometry.getBlockGeometry(0);
  block.forSpatialLocations([&](LatticeR<2> r){
    const auto p=block.getPhysR(r); const T radius=std::hypot(p[0],p[1]);
    int material=1;
    if(radius>o.radius)material=0;
    else if(p[0]<1.51*o.dx||p[1]<1.51*o.dx||radius>o.radius-1.51*o.dx)material=2;
    else if(radius>o.electrolyteRadius) {
      const T cx=(std::floor(p[0]/o.pitch)+.5)*o.pitch;
      const T cy=(std::floor(p[1]/o.pitch)+.5)*o.pitch;
      if(std::hypot(p[0]-cx,p[1]-cy)<=particleRadius)material=2;
    }
    block.set(r,material);
  });
  geometry.getStatistics().getStatisticsStatus()=true;geometry.communicate();geometry.checkForErrors();geometry.print();

  SuperLattice<T,NSD> ns(geometry); SuperLattice<T,CHD> ch(geometry);
  Coupling coupling(Op{},NavierStokes{},ns,Component1{},ch);
  const int resolution=static_cast<int>(std::lround(o.radius/o.dx));
  ns.setUnitConverter<MultiPhaseUnitConverterFromRelaxationTime<T,NSD>>(
    resolution,o.tau,1.,o.radius,o.nu,o.rhoL);
  ch.setUnitConverter(ns.getUnitConverter()); const auto& converter=ns.getUnitConverter();
  const T sigmaLattice=o.sigma/converter.getConversionFactorSurfaceTension();
  const T rhoRatio=o.rhoG/o.rhoL;
  auto fluid=geometry.getMaterialIndicator(1),wall=geometry.getMaterialIndicator(2);
  auto all=geometry.getMaterialIndicator({0,1,2});
  dynamics::set<MultiPhaseIncompressibleBGKdynamics<T,NSD>>(ns,fluid);
  dynamics::set<ConservativePhaseFieldBGKdynamics<T,CHD>>(ch,fluid);
  boundary::set<boundary::BounceBackIncompressible>(ns,wall);
  boundary::set<boundary::PhaseFieldCurvedWall>(ch,wall);
  ch.addPostProcessor<stage::PreCoupling>(fluid,meta::id<RhoStatistics>());
  coupling.restrictTo(fluid); coupling.setParameter<Op::SIGMA>(sigmaLattice);
  coupling.setParameter<Op::W>(o.width); coupling.setParameter<Op::TAUS>({T(.53),o.tau});
  coupling.setParameter<Op::RHOS>({rhoRatio,T(1)}); coupling.setParameter<Op::SWITCH>(1.);
  ns.setParameter<OMEGA>(1/o.tau);ch.setParameter<OMEGA>(1/o.tau);
  ch.setParameter<INTERFACE_WIDTH>(o.width);ch.setParameter<THETA>(o.theta*M_PI/180);
  ch.setParameter<EPSILON>(1.5*o.width);
  auto& communicator=ch.getCommunicator(stage::PreCoupling());
  communicator.requestOverlap(1);communicator.requestField<STATISTIC>();communicator.exchangeRequests();

  RadialPhase phase(o.electrolyteRadius,o.width*o.dx);
  AnalyticalConst2D<T,T> zero(0),one(1),two(2),velocity(Vector<T,2>{0,0});
  AnalyticalConst2D<T,T> rhoGas(rhoRatio),rhoDifference(1-rhoRatio),tauGas(.53),tauDifference(o.tau-.53);
  AnalyticalIdentity2D<T,T> rho(rhoGas+rhoDifference*phase),tau(tauGas+tauDifference*phase);
  fields::set<RHO>(ns,all,rho);fields::set<TAU_EFF>(ns,all,tau);
  fields::set<SCALAR>(ns,all,one);fields::set<SCALAR>(ns,wall,two);
  fields::set<BOUNDARY>(ch,fluid,one);fields::set<BOUNDARY>(ch,wall,two);
  AnalyticalConst2D<T,T> wallTheta(o.theta*M_PI/180);fields::set<THETA>(ch,wall,wallTheta);
  fields::set<OLD_PHIU>(ch,all,velocity);
  momenta::setOrderParameter(ch,all,phase);momenta::setIncompressiblePressure(ns,all,zero);
  momenta::setVelocity(ns,all,velocity);ns.iniEquilibrium(all,zero,velocity);ch.iniEquilibrium(all,phase,velocity);
  ns.initialize();ch.initialize();ch.iniEquilibrium(all,phase,velocity);ns.iniEquilibrium(all,zero,velocity);
  communicator.communicate();ns.setProcessingContext(ProcessingContext::Simulation);ch.setProcessingContext(ProcessingContext::Simulation);
  ch.executePostProcessors(stage::PreCoupling());
  communicator.communicate();

  const bool isMain=singleton::mpi().isMainProcessor();
  std::ofstream csv;
  if(isMain) {
    csv.open(o.output+"/multiphase_diagnostics.csv");
    csv<<std::setprecision(17)<<"step,time_s,mean_phi,min_phi,max_phi,max_speed_lattice,finite\n";
  }
  Phi phi(ch,geometry);SuperLatticePhysVelocity2D<T,NSD> physVelocity(ns,converter);physVelocity.getName()="velocity_phys_m_s";
  SuperGeometryF<T,2> material(geometry);material.getName()="material";
  SuperVTMwriter2D<T> writer("nmc_timecost_multiphase");writer.addFunctor(material);writer.addFunctor(phi);writer.addFunctor(physVelocity);writer.createMasterFile();
  bool finite=true;
  for(long step=0;step<=o.steps;++step){
    if(step){ns.collideAndStream();ch.collideAndStream();communicator.communicate();ch.executePostProcessors(stage::PreCoupling());communicator.communicate();coupling.apply();}
    if(step%o.diagnosticEvery==0||step==o.steps){auto s=sample(geometry,ch);T u=ns.getStatistics().getMaxU();bool ok=std::isfinite(s.mean)&&std::isfinite(s.minimum)&&std::isfinite(s.maximum)&&std::isfinite(u);finite&=ok;if(isMain)csv<<step<<','<<converter.getPhysTime(step)<<','<<s.mean<<','<<s.minimum<<','<<s.maximum<<','<<u<<','<<ok<<'\n';}
    if(step%o.outputEvery==0||step==o.steps)writer.write(static_cast<int>(step));
  }
  auto final=sample(geometry,ch);
  const bool pass=finite&&final.minimum>=-.1&&final.maximum<=1.1&&ns.getStatistics().getMaxU()<.1;
  if(isMain) {
    std::ofstream result(o.output+"/multiphase_validation.csv");
    result<<"status,dx_m,dt_s,steps,physical_time_s,mean_phi,min_phi,max_phi,max_speed_lattice\n"
          <<(pass?"PASS":"FAIL")<<','<<o.dx<<','<<converter.getPhysDeltaT()<<','<<o.steps<<','<<converter.getPhysTime(o.steps)<<','<<final.mean<<','<<final.minimum<<','<<final.maximum<<','<<ns.getStatistics().getMaxU()<<'\n';
    std::cout<<"multiphase_status="<<(pass?"PASS":"FAIL")<<" mean_phi="<<final.mean<<" phi_range=["<<final.minimum<<','<<final.maximum<<"] max_u="<<ns.getStatistics().getMaxU()<<'\n';
  }
  return pass?0:2;
}
