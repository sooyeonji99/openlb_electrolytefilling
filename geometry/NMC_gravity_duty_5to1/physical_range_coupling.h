/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Tim Bingert, Michael Rennick
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#pragma once
// Local adaptation of OpenLB LiangPostProcessor / CorrectedLiangPostProcessor.
// Bounds material interpolation and its density gradient/laplacian. CORRECTED
// enables the native corrected density force and phase-flux correction used by
// risingBubble2d. This operator never clips phase populations or phase mass;
// the optional conservative projection is a separate operation in main.cpp.
namespace olb {
struct PhysicalRangeLiangPostProcessor {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  struct SIGMA       : public descriptors::FIELD_BASE<1> { };
  struct W           : public descriptors::FIELD_BASE<1> { };
  struct TAUS        : public descriptors::FIELD_BASE<2> { };
  struct RHOS        : public descriptors::FIELD_BASE<2> { };
  struct SWITCH      : public descriptors::FIELD_BASE<1> { };
  struct BOUND_MATERIAL : public descriptors::FIELD_BASE<1> { };
  struct CORRECTED : public descriptors::FIELD_BASE<1> { };
  using parameters = meta::list<SIGMA,W,TAUS,RHOS,SWITCH,BOUND_MATERIAL,CORRECTED>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR = typename CELL::template value_t<names::NavierStokes>::descriptor_t;

    auto& cellNS = cells.template get<names::NavierStokes>();
    auto& cellAC = cells.template get<names::Component1>();

    V phi = cellAC.template getFieldComponent<descriptors::STATISTIC>(0);
    const bool boundMaterial = parameters.template get<BOUND_MATERIAL>() > V(0.5);
    const V phiMaterial = boundMaterial ? util::min(V(1),util::max(V(0),phi)) : phi;
    Vector<V,DESCRIPTOR::d> gradMaterialPhi{};
    Vector<V,DESCRIPTOR::d> gradPhi{};
    V laplacePhi = 0, laplaceMaterialPhi = 0;
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      const V phi_i = cellAC.neighbor(descriptors::c<DESCRIPTOR>(iPop)).template getFieldComponent<descriptors::STATISTIC>(0);
      V tcs2 = descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
      gradPhi += phi_i * descriptors::c<DESCRIPTOR>(iPop) * tcs2;
      const V bounded_i = boundMaterial ? util::min(V(1),util::max(V(0),phi_i)) : phi_i;
      gradMaterialPhi += bounded_i * descriptors::c<DESCRIPTOR>(iPop) * tcs2;
      laplaceMaterialPhi += 2*(bounded_i-phiMaterial)*tcs2;
      laplacePhi += 2*(phi_i - phi) * tcs2;
    }

    auto scale = cellNS.template getField<descriptors::SCALAR>();
    auto sigma = parameters.template get<SIGMA>()*scale;
    auto w = parameters.template get<W>();
    auto rho_v = parameters.template get<RHOS>()[0];
    auto rho_l = parameters.template get<RHOS>()[1];
    V rho = rho_v + (rho_l-rho_v)*phiMaterial;
    auto tau_v = parameters.template get<TAUS>()[0];
    auto tau_l = parameters.template get<TAUS>()[1];
    V tau = (tau_v + (tau_l-tau_v)*phiMaterial)*scale + (1-scale)*5.;

    // Computation and storage of forces
    Vector<V,DESCRIPTOR::d> forceNS{};
    V u[DESCRIPTOR::d] {};
    V k = 1.5*sigma*w;
    V beta = 12*sigma/w;
    V mu = 4*beta*phi*(phi-1.)*(phi-0.5)-k*laplacePhi;
    forceNS += mu*gradPhi;
    auto externalForce = cellNS.template getField<descriptors::EXTERNAL_FORCE>();
    cellNS.template setField<descriptors::NABLARHO>((rho_l-rho_v)*gradMaterialPhi);
    cellNS.template setField<descriptors::LAPLACERHO>((rho_l-rho_v)*laplaceMaterialPhi);
    cellNS.template setField<descriptors::TAU_EFF>(tau);
    cellNS.template setField<descriptors::RHO>(rho);
    const bool corrected = parameters.template get<CORRECTED>() > V(0.5);
    const V forceDensity = corrected ? rho+(rho_l-rho_v)*laplaceMaterialPhi/V(4) : rho;
    cellNS.template setField<descriptors::FORCE>(externalForce + forceNS/forceDensity);
    cellNS.computeU(u);

    Vector<V,DESCRIPTOR::d> forceAC{};
    Vector<V,DESCRIPTOR::d> old_phiU{};
    Vector<V,DESCRIPTOR::d> n = normalizeWithEpsilon(gradPhi,V(1),V(1e-7));
    Vector<V,DESCRIPTOR::d> phiU{};
    auto velo_switch = parameters.template get<SWITCH>();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      u[iD] *= velo_switch;
      phiU[iD] = phi*u[iD];
    }
    old_phiU = cellAC.template getField<descriptors::OLD_PHIU>();
    cellAC.template setField<descriptors::OLD_PHIU>(phiU);
    V lambda = 4*phi*(1.-phi)/w;
    if(corrected) {
      // Same phase correction as OpenLB CorrectedLiangPostProcessor.
      V Y1 = util::min(util::max(V(0),util::min(V(1),(phi-V(0.8))/V(1e-4))),util::max(V(0),util::min(V(1),(V(1)-phi)/V(1e-4))));
      V Y2 = util::min(util::max(V(0),util::min(V(1),phi/V(1e-4))),util::max(V(0),util::min(V(1),(V(0.2)-phi)/V(1e-4))));
      V laplacePhiNorm = laplacePhi/(util::abs(laplacePhi)+V(1e-7));
      lambda += -V(22.5)*util::norm<DESCRIPTOR::d>(gradPhi)/w*((V(1)+laplacePhiNorm)*Y1+(V(1)-laplacePhiNorm)*Y2);
    }
    forceAC += (phiU - old_phiU) + lambda*n/descriptors::invCs2<V,DESCRIPTOR>();
    cellAC.template setField<descriptors::FORCE>(forceAC);
    cellAC.template setField<descriptors::VELOCITY>(u);
  }
};


} // namespace olb
