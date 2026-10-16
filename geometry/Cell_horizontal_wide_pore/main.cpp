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
#include "geometry.h"
#include "physical_range_coupling.h"
#include "phase_projection.h"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::names;
using T = double;
using NSD = D2Q9<RHO,NABLARHO,LAPLACERHO,FORCE,EXTERNAL_FORCE,TAU_EFF,STATISTIC,SCALAR>;
using CHD = D2Q9<CONV_POPS,FORCE,SOURCE,VELOCITY,OLD_PHIU,STATISTIC,THETA,BOUNDARY>;
using Op = PhysicalRangeLiangPostProcessor;
using Maps = meta::map<NavierStokes,descriptors::VALUED_DESCRIPTOR<T,NSD>,
                       Component1,descriptors::VALUED_DESCRIPTOR<T,CHD>>;
using Coupling = SuperLatticeCoupling<Op,Maps>;

struct Options {
  bool openOuter=true, boundMaterial=true, projectPhi=true;
  std::string collision="crlb",sourceMode="legacy";
  T phaseTau=.53,magic=3./16;
  T radius=cell::width, dx=5e-6;
  T porosity=.40, tau=.53, nu=1e-6, rhoL=1000, rhoG=1.2;
  T nuGas=1.5e-5;
  T sigma=1e-5, theta=30, width=4;
  T reservoirPressure=2500, pressureRamp=.001;
  long steps=30000, outputEvery=100, diagnosticEvery=1;
  std::string output="results/air_theta30_p2500_30000";
};

Options parse(int argc,char**argv) {
  Options o;
  for (int i=1;i<argc;++i) {
    if (i+1>=argc) throw std::runtime_error("missing option value");
    std::string key=argv[i],value=argv[++i];
    if(key=="--dx-um")o.dx=std::stod(value)*1e-6;
    else if(key=="--reservoir-pressure-pa")o.reservoirPressure=std::stod(value);
    else if(key=="--outer-boundary") {if(value!="open"&&value!="closed")throw std::runtime_error("outer-boundary must be open or closed");o.openOuter=value=="open";}
    else if(key=="--pressure-ramp-s")o.pressureRamp=std::stod(value);
    else if(key=="--ns-collision")o.collision=value;
    else if(key=="--source-mode")o.sourceMode=value;
    else if(key=="--phase-tau")o.phaseTau=std::stod(value);
    else if(key=="--nu-liquid")o.nu=std::stod(value);
    else if(key=="--phase-projection") {if(value!="off"&&value!="conservative")throw std::runtime_error("phase-projection must be off or conservative");o.projectPhi=value=="conservative";}
    else if(key=="--material-interpolation") {if(value!="linear"&&value!="bounded")throw std::runtime_error("material-interpolation must be linear or bounded");o.boundMaterial=value=="bounded";}
    else if(key=="--rho-gas")o.rhoG=std::stod(value);
    else if(key=="--nu-gas")o.nuGas=std::stod(value);
    else if(key=="--contact-angle-deg")o.theta=std::stod(value);
    else if(key=="--sigma")o.sigma=std::stod(value);
    else if(key=="--tau")o.tau=std::stod(value);
    else if(key=="--steps")o.steps=std::stol(value);
    else if(key=="--output-every")o.outputEvery=std::stol(value);
    else if(key=="--diagnostic-every")o.diagnosticEvery=std::stol(value);
    else if(key=="--output")o.output=value;
    else throw std::runtime_error("unknown option: "+key);
  }
  if(!std::isfinite(o.rhoG)||o.rhoG<=0||o.rhoG>=o.rhoL)throw std::runtime_error("rho-gas must be finite and between zero and liquid density");
  if(!std::isfinite(o.theta)||o.theta<=0||o.theta>=180)throw std::runtime_error("contact-angle-deg must be finite and in (0,180)");
  if(!std::isfinite(o.nu)||o.nu<=0||!std::isfinite(o.nuGas)||o.nuGas<=0)throw std::runtime_error("viscosities must be finite and positive");
  if(!std::isfinite(o.tau)||o.tau<=.5)throw std::runtime_error("tau must exceed 0.5");
  if (!std::isfinite(o.dx) || o.dx<=0 || o.dx>5e-6 || !std::isfinite(o.reservoirPressure) || o.steps<0 || o.outputEvery<=0 || o.diagnosticEvery<=0) throw std::runtime_error("invalid options: dx must be in (0,5] um; intervals positive");
  if(!std::isfinite(o.pressureRamp)||o.pressureRamp<0)throw std::runtime_error("pressure-ramp-s must be finite and nonnegative");
  if(o.collision!="bgk"&&o.collision!="trt"&&o.collision!="itrt"&&o.collision!="rlb"&&o.collision!="crlb")throw std::runtime_error("ns-collision: bgk, trt, itrt, rlb or crlb");
  if(o.sourceMode!="legacy"&&o.sourceMode!="equilibrium"&&o.sourceMode!="finite"&&o.sourceMode!="decoupled")throw std::runtime_error("source-mode: legacy, equilibrium, finite or decoupled");
  if(o.phaseTau==0)o.phaseTau=o.tau;
  if(!std::isfinite(o.phaseTau)||o.phaseTau<=.5||!std::isfinite(o.sigma)||o.sigma<0)throw std::runtime_error("invalid phase tau or sigma");
  if(o.sourceMode=="finite"&&o.reservoirPressure!=0)throw std::runtime_error("finite reservoir comparison requires pressure-pa 0");
  return o;
}

class SourcePhase : public AnalyticalF2D<T,T> {
  T _w;
public:
  SourcePhase(T w):AnalyticalF2D<T,T>(1),_w(w){}
  bool operator()(T out[],const T p[]) override {
    const T distance=std::max(std::abs(p[0]-cell::width/2)-cell::sourceWidth/2,p[1]-cell::sourceHeight);
    out[0]=.5*(1-std::tanh(2*distance/_w)); return true;
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
    const int material=_geometry.getBlockGeometry(load.loc(p[0])).get(p[1],p[2]);
    if(material!=1 && material!=3 && !(material>=4&&material<=7)) {
      out[0]=std::numeric_limits<T>::quiet_NaN();
      return true;
    }
    out[0]=_lattice.getBlock(load.loc(p[0])).get(p[1],p[2]).template getFieldComponent<STATISTIC>(0);
    return true;
  }
};

struct Sample { T mean{},minimum{1e9},maximum{-1e9},cells{},rhoMin{1e9},tauMin{1e9},outside{},overshootL1{}; };
Sample sample(SuperGeometry<T,2>& geometry,SuperLattice<T,CHD>& ch,SuperLattice<T,NSD>& ns) {
  Sample s;
  for (int iC=0; iC<ch.getLoadBalancer().size(); ++iC) {
    auto& g=geometry.getBlockGeometry(iC); auto& lattice=ch.getBlock(iC);
    g.forSpatialLocations([&](LatticeR<2> r){
      int material=g.get(r);
      if(lattice.isPadding(r)||(material!=1&&material!=3&&!(material>=4&&material<=7)))return;
      T phi=lattice.get(r).template getFieldComponent<STATISTIC>(0);
      if(!std::isfinite(phi))s.mean=std::numeric_limits<T>::quiet_NaN();
      if(phi < -1e-12 || phi > 1.+1e-12)++s.outside;
      s.overshootL1+=std::abs(phi-std::clamp(phi,T(0),T(1)));
      T rho=ns.getBlock(iC).get(r).template getFieldComponent<RHO>(0);
      T tau=ns.getBlock(iC).get(r).template getFieldComponent<TAU_EFF>(0);
      if(!std::isfinite(rho)||!std::isfinite(tau))s.mean=std::numeric_limits<T>::quiet_NaN();
      s.rhoMin=std::min(s.rhoMin,rho);s.tauMin=std::min(s.tauMin,tau);
      s.minimum=std::min(s.minimum,phi);s.maximum=std::max(s.maximum,phi);
      if(material==1){s.mean+=phi;++s.cells;}
    });
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(s.outside,MPI_SUM);
  singleton::mpi().reduceAndBcast(s.overshootL1,MPI_SUM);
  singleton::mpi().reduceAndBcast(s.rhoMin,MPI_MIN);
  singleton::mpi().reduceAndBcast(s.tauMin,MPI_MIN);
  singleton::mpi().reduceAndBcast(s.mean,MPI_SUM);
  singleton::mpi().reduceAndBcast(s.cells,MPI_SUM);
  singleton::mpi().reduceAndBcast(s.minimum,MPI_MIN);
  singleton::mpi().reduceAndBcast(s.maximum,MPI_MAX);
#endif
  if(s.cells)s.mean/=s.cells;
  return s;
}

PhaseProjectionStats projectPhaseLattice(SuperGeometry<T,2>& geometry,SuperLattice<T,CHD>& ch) {
  struct Location {int block;LatticeR<2> r;};
  std::vector<Location> locations;std::vector<double> values;
  locations.reserve(125000);values.reserve(125000);
  for(int b=0;b<ch.getLoadBalancer().size();++b){auto& g=geometry.getBlockGeometry(b);auto& lattice=ch.getBlock(b);
    g.forSpatialLocations([&](LatticeR<2> r){int m=g.get(r);if(lattice.isPadding(r)||(m!=1&&m!=4&&m!=5&&m!=7))return;
      locations.push_back({b,r});values.push_back(lattice.get(r).computeRho());});
  }
  auto stats=projectPhase(values,[](double v){
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(v,MPI_SUM);
#endif
    return v;
  });
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(stats.rawMin,MPI_MIN);singleton::mpi().reduceAndBcast(stats.rawMax,MPI_MAX);
#endif
  for(size_t i=0;i<values.size();++i){auto cell=ch.getBlock(locations[i].block).get(locations[i].r);T delta=values[i]-cell.computeRho();
    // Preserve the first population moment: sum_i w_i*c_i = 0.
    for(int q=0;q<CHD::q;++q)cell[q]+=descriptors::t<T,CHD>(q)*delta;
  }
  return stats;
}

int main(int argc,char**argv) {
  initialize(&argc,&argv); try { const Options o=parse(argc,argv);
  std::filesystem::create_directories(o.output);
  singleton::directories().setOutputDir(o.output+"/");
  cell::Geometry layout(o.dx,o.porosity,50e-6,o.openOuter);
  if(singleton::mpi().isMainProcessor()) layout.write(o.output);
  IndicatorCuboid2D<T> domain({cell::width,layout.latticeHeight},{0,0});
  CuboidDecomposition2D<T> cuboids(domain, o.dx, singleton::mpi().getSize());
  HeuristicLoadBalancer<T> load(cuboids);
  SuperGeometry<T,2> geometry(cuboids,load,2); geometry.rename(0,1);
  for(int iC=0;iC<load.size();++iC) {
    auto& block=geometry.getBlockGeometry(iC);
    block.forSpatialLocations([&](LatticeR<2> r){
      const auto p=block.getPhysR(r);
      block.set(r,layout.material(p[0],p[1]));
    });
  }
  geometry.getStatistics().getStatisticsStatus()=true;geometry.communicate();geometry.checkForErrors();geometry.print();

  SuperLattice<T,NSD> ns(geometry); SuperLattice<T,CHD> ch(geometry);
  Coupling coupling(Op{},NavierStokes{},ns,Component1{},ch);
  const int resolution=static_cast<int>(std::lround(o.radius/o.dx));
  ns.setUnitConverter<MultiPhaseUnitConverterFromRelaxationTime<T,NSD>>(
    resolution,o.tau,1.,o.radius,o.nu,o.rhoL);
  ch.setUnitConverter(ns.getUnitConverter()); const auto& converter=ns.getUnitConverter();
  const T sigmaLattice=o.sigma/converter.getConversionFactorSurfaceTension();
  const T rhoRatio=o.rhoG/o.rhoL;
  const T tauGasValue=.5+(o.nuGas/o.nu)*(o.tau-.5);
  auto pore=geometry.getMaterialIndicator(1),wall=geometry.getMaterialIndicator(2);
  auto reservoir=geometry.getMaterialIndicator(3),fluid=geometry.getMaterialIndicator({1,3,4,5,6,7});
  auto all=geometry.getMaterialIndicator({0,1,2,3,4,5,6,7});
  if(o.collision=="itrt")dynamics::set<MultiPhaseIncompressibleInterfaceTRTdynamics<T,NSD>>(ns,fluid);
  else if(o.collision=="crlb")dynamics::set<CorrectedMultiPhaseIncompressibleRLBdynamics<T,NSD>>(ns,fluid);
  else if(o.collision=="rlb")dynamics::set<MultiPhaseIncompressibleRLBdynamics<T,NSD>>(ns,fluid);
  else if(o.collision=="trt")dynamics::set<MultiPhaseIncompressibleTRTdynamics<T,NSD>>(ns,fluid);
  else dynamics::set<MultiPhaseIncompressibleBGKdynamics<T,NSD>>(ns,fluid);
  dynamics::set<ConservativePhaseFieldBGKdynamics<T,CHD>>(ch,fluid);
  boundary::set<boundary::BounceBackIncompressible>(ns,wall);
  boundary::set<boundary::PhaseFieldCurvedWall>(ch,wall);
  if(o.openOuter) {
    auto outlets=geometry.getMaterialIndicator({4,5,6,7});
    boundary::set<boundary::IncompressibleZouHePressure>(ns,outlets);
    setConvectivePhaseFieldBoundary<T,CHD>(ch,outlets);
    ch.addPostProcessor<stage::InitOutlet>(geometry.getMaterialIndicator(4),meta::id<SetOutletCells2D<-1,0>>());
    ch.addPostProcessor<stage::InitOutlet>(geometry.getMaterialIndicator(5),meta::id<SetOutletCells2D<1,0>>());
    ch.addPostProcessor<stage::InitOutlet>(geometry.getMaterialIndicator(7),meta::id<SetOutletCells2D<0,1>>());
  }
  ch.addPostProcessor<stage::PreCoupling>(fluid,meta::id<RhoStatistics>());
  if(o.sourceMode=="decoupled")coupling.restrictTo(geometry.getMaterialIndicator({1,4,5,6,7}));
  else coupling.restrictTo(fluid); coupling.setParameter<Op::SIGMA>(sigmaLattice);
  coupling.setParameter<Op::CORRECTED>(o.collision=="crlb"?1.:0.);
  coupling.setParameter<Op::BOUND_MATERIAL>(o.boundMaterial?1.:0.);
  coupling.setParameter<Op::W>(o.width); coupling.setParameter<Op::TAUS>({tauGasValue,o.tau});
  coupling.setParameter<Op::RHOS>({rhoRatio,T(1)}); coupling.setParameter<Op::SWITCH>(1.);
  ns.setParameter<OMEGA>(1/o.tau);ch.setParameter<OMEGA>(1/o.phaseTau);
  if(o.collision=="trt")ns.setParameter<collision::TRT::MAGIC>(o.magic);
  if(o.collision=="itrt") {
    ns.setParameter<collision::ITRT::TAU_MINUS>(1.5);
    ns.setParameter<collision::ITRT::MAXNABLARHO>((1-rhoRatio)/o.width);
  }
  ch.setParameter<INTERFACE_WIDTH>(o.width);ch.setParameter<THETA>(o.theta*M_PI/180);
  ch.setParameter<EPSILON>(1.5*o.width);
  auto& communicator=ch.getCommunicator(stage::PreCoupling());
  communicator.requestOverlap(1);communicator.requestField<STATISTIC>();communicator.exchangeRequests();

  SourcePhase phase(o.width*o.dx);
  AnalyticalConst2D<T,T> zero(0),one(1),two(2),velocity(Vector<T,2>{0,0});
  auto pressureAt=[&](T time) {return o.reservoirPressure*(o.pressureRamp>0?std::clamp(time/o.pressureRamp,T(0),T(1)):T(1));};
  auto applyPressure=[&](T time) {
    AnalyticalConst2D<T,T> p(pressureAt(time)/converter.getConversionFactorPressure());
    ns.iniEquilibrium(reservoir,p,velocity);
  };
  AnalyticalConst2D<T,T> rhoGas(rhoRatio),rhoDifference(1-rhoRatio),tauGas(tauGasValue),tauDifference(o.tau-tauGasValue);
  AnalyticalIdentity2D<T,T> rho(rhoGas+rhoDifference*phase),tau(tauGas+tauDifference*phase);
  fields::set<RHO>(ns,all,rho);fields::set<TAU_EFF>(ns,all,tau);
  fields::set<SCALAR>(ns,all,one);fields::set<SCALAR>(ns,wall,two);
  fields::set<BOUNDARY>(ch,fluid,one);fields::set<BOUNDARY>(ch,wall,two);
  AnalyticalConst2D<T,T> wallTheta(o.theta*M_PI/180);fields::set<THETA>(ch,wall,wallTheta);
  fields::set<OLD_PHIU>(ch,all,velocity);
  momenta::setOrderParameter(ch,all,phase);momenta::setIncompressiblePressure(ns,all,zero);
  momenta::setVelocity(ns,all,velocity);ns.iniEquilibrium(all,zero,velocity);ch.iniEquilibrium(all,phase,velocity);
  ns.initialize();ch.initialize();ch.iniEquilibrium(all,phase,velocity);ns.iniEquilibrium(all,zero,velocity);
  auto applySource=[&](T time) {
    if(o.sourceMode=="finite")return;
    applyPressure(time);
    if(o.sourceMode=="equilibrium"||o.sourceMode=="decoupled") {
      ch.iniEquilibrium(reservoir,one,velocity);
      fields::set<OLD_PHIU>(ch,reservoir,velocity);
      fields::set<VELOCITY>(ch,reservoir,velocity);
      fields::set<FORCE>(ch,reservoir,velocity);
      if(o.sourceMode=="decoupled") {
        AnalyticalConst2D<T,T> liquidTau(o.tau);
        fields::set<RHO>(ns,reservoir,one);
        fields::set<TAU_EFF>(ns,reservoir,liquidTau);
        fields::set<FORCE>(ns,reservoir,velocity);
        fields::set<NABLARHO>(ns,reservoir,velocity);
        fields::set<LAPLACERHO>(ns,reservoir,zero);
      }
    } else momenta::setOrderParameter(ch,reservoir,one);
  };
  momenta::setOrderParameter(ch,reservoir,one); // identical initial liquid bath in all ablations
  if(o.sourceMode!="finite")applySource(0);
  if(o.openOuter) {
    momenta::setIncompressiblePressure(ns,geometry.getMaterialIndicator({4,5,6,7}),zero);
    ch.executePostProcessors(stage::InitOutlet());
  }
  communicator.communicate();ns.setProcessingContext(ProcessingContext::Simulation);ch.setProcessingContext(ProcessingContext::Simulation);
  ch.executePostProcessors(stage::PreCoupling());
  communicator.communicate();

  if(singleton::mpi().isMainProcessor()) {
    std::ofstream cfg(o.output+"/run_config.csv");
    cfg<<std::setprecision(17)<<"key,value\nphase_projection,"<<(o.projectPhi?"conservative":"off")<<"\nmaterial_interpolation,"<<(o.boundMaterial?"bounded":"linear")<<"\nns_collision,"<<o.collision<<"\nsource_mode,"<<o.sourceMode<<"\nphase_tau,"<<o.phaseTau<<"\ntrt_magic,"<<o.magic<<"\nitrt_tau_minus,1.5\nitrt_max_rho_gradient,"<<(1-rhoRatio)/o.width<<"\nrequested_steps,"<<o.steps<<"\nouter_boundary,"<<(o.openOuter?"open":"closed")<<"\nbottom_boundary,wall\noutlet_pressure_pa,0\nmax_speed_lattice_limit,disabled\nwidth_m,"<<cell::width<<"\nheight_m,"<<layout.latticeHeight<<"\nrequested_height_m,"<<cell::height<<"\nsource_width_m,"<<cell::sourceWidth<<"\nsource_height_m,"<<cell::sourceHeight<<"\ngeometry_version,random_mixed_disks_wide_pore_v1\nparticle_min_diameter_m,"<<cell::Geometry::minDiameter<<"\nparticle_max_diameter_m,"<<cell::Geometry::maxDiameter<<"\nparticle_min_gap_m,"<<cell::Geometry::minGap<<"\ntarget_porosity,"<<o.porosity<<"\nseed,"<<cell::seed<<"\ndx_m,"<<o.dx<<"\ndt_s,"<<converter.getPhysDeltaT()<<"\nrho_liquid_kg_m3,"<<o.rhoL<<"\nrho_gas_kg_m3,"<<o.rhoG<<"\nnu_liquid_m2_s,"<<o.nu<<"\ntau_liquid,"<<o.tau<<"\ntau_gas,"<<tauGasValue<<"\nnu_gas_m2_s,"<<o.nuGas<<"\nsigma_N_m,"<<o.sigma<<"\ncontact_angle_deg,"<<o.theta<<"\npressure_ramp_s,"<<o.pressureRamp<<"\npressure_ramp_shape,linear\npressure_pa,"<<o.reservoirPressure<<"\npressure_conversion_Pa,"<<converter.getConversionFactorPressure()<<"\npressure_lattice,"<<o.reservoirPressure/converter.getConversionFactorPressure()<<"\n";
  }
  const bool isMain=singleton::mpi().isMainProcessor();
  std::ofstream csv;
  if(isMain) {
    csv.open(o.output+"/multiphase_diagnostics.csv");
    csv<<std::setprecision(17)<<"step,time_s,reservoir_pressure_pa,pore_mean_phi,min_phi,max_phi,max_speed_lattice,finite,min_rho_lattice,min_tau_eff,out_of_range_phi_cells,phi_overshoot_l1,raw_min_phi,raw_max_phi,phase_projection_l1,phase_projection_mass_residual,phase_projection_cells\n";
  }
  SuperLatticePhysIncPressure2D<T,NSD> pressure(ns,converter);pressure.getName()="pressure_phys_Pa";
  Phi phi(ch,geometry);SuperLatticePhysVelocity2D<T,NSD> physVelocity(ns,converter);physVelocity.getName()="velocity_phys_m_s";
  SuperGeometryF<T,2> material(geometry);material.getName()="material";
  SuperVTMwriter2D<T> writer("cell_horizontal");writer.addFunctor(pressure);writer.addFunctor(material);writer.addFunctor(phi);writer.addFunctor(physVelocity);writer.createMasterFile();
  long completed=0;
  bool finite=true;
  PhaseProjectionStats projection;projection.rawMin=0;projection.rawMax=1;
  for(long step=0;step<=o.steps;++step){
    if(step){applySource(converter.getPhysTime(step-1));ns.collideAndStream();ch.collideAndStream();applySource(converter.getPhysTime(step));if(o.projectPhi)projection=projectPhaseLattice(geometry,ch);communicator.communicate();ch.executePostProcessors(stage::PreCoupling());communicator.communicate();coupling.apply();}
    completed=step;
    if(step%o.diagnosticEvery==0||step%o.outputEvery==0||step==o.steps){auto s=sample(geometry,ch,ns);T u=ns.getStatistics().getMaxU();bool ok=std::isfinite(s.mean)&&std::isfinite(s.minimum)&&std::isfinite(s.maximum)&&std::isfinite(u);ok=ok&&s.minimum>=-.1&&s.maximum<=1.1&&s.rhoMin>0&&s.tauMin>.5;finite&=ok;if(isMain)csv<<step<<','<<converter.getPhysTime(step)<<','<<pressureAt(converter.getPhysTime(step))<<','<<s.mean<<','<<s.minimum<<','<<s.maximum<<','<<u<<','<<ok<<','<<s.rhoMin<<','<<s.tauMin<<','<<s.outside<<','<<s.overshootL1<<','<<(o.projectPhi?projection.rawMin:s.minimum)<<','<<(o.projectPhi?projection.rawMax:s.maximum)<<','<<projection.l1<<','<<(projection.massAfter-projection.massBefore)<<','<<projection.changed<<'\n';csv.flush();}
    if(step%o.outputEvery==0||step==o.steps)writer.write(static_cast<int>(step));
    if(!finite) {
      if(step%o.outputEvery!=0 && step!=o.steps)writer.write(static_cast<int>(step));
      break;
    }
  }
  auto final=sample(geometry,ch,ns);
  const bool pass=finite&&final.minimum>=-.1&&final.maximum<=1.1&&std::isfinite(ns.getStatistics().getMaxU());
  if(isMain) {
    std::ofstream result(o.output+"/multiphase_validation.csv");
    result<<"status,dx_m,dt_s,steps,physical_time_s,reservoir_pressure_pa,pore_mean_phi,min_phi,max_phi,max_speed_lattice\n"
          <<(pass?"PASS":"FAIL")<<','<<o.dx<<','<<converter.getPhysDeltaT()<<','<<completed<<','<<converter.getPhysTime(completed)<<','<<pressureAt(converter.getPhysTime(completed))<<','<<final.mean<<','<<final.minimum<<','<<final.maximum<<','<<ns.getStatistics().getMaxU()<<'\n';
    std::cout<<"multiphase_status="<<(pass?"PASS":"FAIL")<<" mean_phi="<<final.mean<<" phi_range=["<<final.minimum<<','<<final.maximum<<"] max_u="<<ns.getStatistics().getMaxU()<<'\n';
  }
  return pass?0:2;
  } catch(const std::exception& e) {std::cerr<<"Cell_horizontal_wide_pore: "<<e.what()<<std::endl;return 1;}
}
