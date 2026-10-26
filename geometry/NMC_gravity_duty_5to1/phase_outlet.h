#pragma once
// Local open phase boundary. Match OpenLB's convective update for outward flow;
// prescribe ambient gas only on actual inward flow, sampled one cell inside.
namespace olb {
template<int nx,int ny> struct BackflowPhaseA {
  static constexpr OperatorScope scope=OperatorScope::PerCellWithParameters;
  using parameters=meta::list<descriptors::MAX_VELOCITY>;
  int getPriority()const{return 0;}
  template<class CELL,class PARAMETERS> void apply(CELL& cell,PARAMETERS& parameters) any_platform {
    using V=typename CELL::value_t;using D=typename CELL::descriptor_t;
    auto inner=cell.neighbor({-nx,-ny});
    auto innerU=inner.template getField<descriptors::VELOCITY>();
    auto saved=cell.template getField<descriptors::CONV_POPS>();
    const bool inflow=V(nx)*innerU[0]+V(ny)*innerU[1]<V(-1e-12);
    const V speed=parameters.template get<descriptors::MAX_VELOCITY>();
    V u[D::d]{};cell.computeU(u);
    for(int q=0;q<D::q;++q){
      saved[q]=inflow?equilibrium<D>::firstOrder(q,V(0),u):(saved[q]+speed*inner[q])/(V(1)+speed);
      cell[q]=saved[q];
    }
    cell.template setField<descriptors::CONV_POPS>(saved);
  }
};
template<int nx,int ny> struct BackflowPhaseB {
  static constexpr OperatorScope scope=OperatorScope::PerCellWithParameters;
  using parameters=meta::list<descriptors::MAX_VELOCITY>;
  int getPriority()const{return 0;}
  template<class CELL,class PARAMETERS> void apply(CELL& cell,PARAMETERS& parameters) any_platform {
    using V=typename CELL::value_t;using D=typename CELL::descriptor_t;
    auto saved=cell.template getField<descriptors::CONV_POPS>();
    for(int q=0;q<D::q;++q)cell[q]=saved[q];
    auto phi=cell.template getField<descriptors::STATISTIC>();phi[0]=cell.computeRho();
    cell.template setField<descriptors::STATISTIC>(phi);
    auto innerU=cell.neighbor({-nx,-ny}).template getField<descriptors::VELOCITY>();
    const bool inflow=V(nx)*innerU[0]+V(ny)*innerU[1]<V(-1e-12);
    const V speed=parameters.template get<descriptors::MAX_VELOCITY>();
    auto ghost=cell.neighbor({nx,ny}).template getField<descriptors::STATISTIC>();
    ghost[0]=inflow?V(0):(ghost[0]+speed*phi[0])/(V(1)+speed);
    // Bound only ghost reconstruction, as in native phase outlet. Preserve NaNs.
    if(ghost[0]>V(1))ghost[0]=V(1);if(ghost[0]<V(0))ghost[0]=V(0);
    cell.neighbor({nx,ny}).template setField<descriptors::STATISTIC>(ghost);
    Vector<int,D::d> tangent{ny,-nx};
    if(cell.neighbor(tangent).template getField<descriptors::BOUNDARY>()==2)
      cell.neighbor(Vector<int,D::d>{nx,ny}+tangent).template setField<descriptors::STATISTIC>(ghost);
    tangent*=-1;
    if(cell.neighbor(tangent).template getField<descriptors::BOUNDARY>()==2)
      cell.neighbor(Vector<int,D::d>{nx,ny}+tangent).template setField<descriptors::STATISTIC>(ghost);
  }
};
}
