#ifndef NMC_SIMULATION_H
#define NMC_SIMULATION_H
#include "multiphase.h"
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <array>

namespace nmc {
inline void simulate(MultiphaseSystem& s, olb::SuperGeometry<T,2>& geometry,
                     int steps, int outputEvery) {
    PhaseOutput2D phi(s.cahnHilliard);
    olb::SuperGeometryMaterialOutput2D<T> material(geometry);
    olb::SuperLatticeVelocity2D<T,NSDescriptor> velocity(s.navierStokes);
    velocity.getName() = "velocity_lattice";
    olb::SuperLatticeExternalScalarField2D<T,CHDescriptor,olb::descriptors::CHEM_POTENTIAL> mu(s.cahnHilliard);
    mu.getName() = "chemical_potential";
    olb::SuperVTMwriter2D<T> writer("electrolyte_filling");
    writer.addFunctor(material); writer.addFunctor(phi);
    writer.addFunctor(velocity); writer.addFunctor(mu);
    writer.createMasterFile();
    std::ofstream csv("filling_history.csv");
    if (!csv) throw std::runtime_error("Cannot open filling_history.csv");
    csv << "step,time_s,electrode_saturation,phi_min,phi_max,max_speed_lattice\n" << std::setprecision(12);
    std::ofstream diagnostics("infiltration_diagnostics.csv"), profiles("height_phi_profile.csv");
    if (!diagnostics || !profiles) throw std::runtime_error("Cannot open diagnostics CSV files");
    diagnostics << "step,time_s,liquid_area_raw_m2,liquid_area_change_m2,max_phi_half_crossing_y_m,highest_liquid_node_y_m,inlet_volume_flux_m2_s,outlet_volume_flux_m2_s,inlet_liquid_advective_flux_m2_s,outlet_liquid_advective_flux_m2_s,cumulative_inlet_liquid_advective_m2,cumulative_outlet_liquid_advective_m2,advective_balance_residual_m2\n" << std::setprecision(12);
    profiles << "step,time_s,y_m,pore_nodes,mean_phi_raw,mean_phi_clamped\n" << std::setprecision(12);
    std::ofstream regions("regional_liquid_balance.csv");
    if (!regions) throw std::runtime_error("Cannot open regional_liquid_balance.csv");
    regions << "step,time_s,reservoir_area_raw_m2,electrode_area_raw_m2,gas_buffer_area_raw_m2,reservoir_change_m2,electrode_change_m2,gas_buffer_change_m2,total_change_m2,total_relative_change,reservoir_to_electrode_advective_flux_m2_s,electrode_to_buffer_advective_flux_m2_s,cumulative_reservoir_to_electrode_advective_m2,cumulative_electrode_to_buffer_advective_m2,reservoir_advective_residual_m2,electrode_advective_residual_m2,gas_buffer_advective_residual_m2\n" << std::setprecision(17);
    std::array<T,3> initialRegions{};
    std::array<T,2> previousInternal{}, cumulativeInternal{};
    T initialArea=0., previousIn=0., previousOut=0., cumulativeIn=0., cumulativeOut=0.;
    const T dt=s.navierStokes.getUnitConverter().getPhysDeltaT();
    const T velocityConversion=s.navierStokes.getUnitConverter().getConversionFactorVelocity();
    for (int step=0; step<=steps; ++step) {
        if (step>0) {
            s.navierStokes.collideAndStream();
            s.cahnHilliard.collideAndStream();
            auto& comm=s.cahnHilliard.getCommunicator(olb::stage::PreCoupling());
            comm.communicate();
            s.cahnHilliard.executePostProcessors(olb::stage::PreCoupling());
            comm.communicate();
            s.cahnHilliard.executePostProcessors(olb::stage::ChemPotCalc());
            comm.communicate();
            s.coupling.apply();
        }
        // Inspect every step; saturation excludes reservoir and gas buffer.
        T lo=std::numeric_limits<T>::max(), hi=-lo, maxSpeed=0., sum=0.;
        std::size_t count=0;
        std::array<T,3> regionalArea{};
        std::array<T,2> internalFlux{};
        T liquidArea=0., inletVolume=0., outletVolume=0., inletLiquid=0., outletLiquid=0.;
        T crossing=std::numeric_limits<T>::quiet_NaN(), liquidNode=crossing;
        struct HeightRow { std::size_t count=0; T raw=0., clamped=0.; };
        std::map<int,HeightRow> heightRows;
        auto updateHeight=[](T& maximum, T value) {
            if (!std::isfinite(maximum) || value>maximum) maximum=value;
        };
        for (int iC=0; iC<geometry.getLoadBalancer().size(); ++iC) {
            auto& g=geometry.getBlockGeometry(iC);
            auto& ch=s.cahnHilliard.getBlock(iC);
            auto& ns=s.navierStokes.getBlock(iC);
            g.forSpatialLocations([&](olb::LatticeR<2> r) {
                const int m=g.get(r);
                if (m!=param::openlbPoreMaterial && m!=param::liquidSupplyMaterial && m!=param::gasOutletMaterial) return;
                const T p=ch.get(r).template getFieldComponent<olb::descriptors::STATISTIC>(0);
                T u[2]{}; ns.get(r).computeU(u);
                const T speed=std::hypot(u[0],u[1]);
                const T chem=ch.get(r).template getField<olb::descriptors::CHEM_POTENTIAL>();
                if (!std::isfinite(p) || !std::isfinite(speed) || !std::isfinite(chem))
                    throw std::runtime_error("Non-finite multiphase field at step " + std::to_string(step));
                lo=std::min(lo,p); hi=std::max(hi,p); maxSpeed=std::max(maxSpeed,speed);
                const auto x=g.getPhysR(r);
                if (m==param::openlbPoreMaterial) {
                    const T area=p*param::dx*param::dx;
                    liquidArea+=area;
                    const T eps=param::dx*1e-7;
                    const int region=x[1]<-eps ? 0 : (x[1]<=param::ly+eps ? 1 : 2);
                    regionalArea[region]+=area;
                    // Sample liquid advection on the existing entrance/exit
                    // pore-node planes. These are not complete CH fluxes.
                    const T q=p*u[1]*velocityConversion*param::dx;
                    if (std::abs(x[1])<=eps) internalFlux[0]+=q;
                    if (std::abs(x[1]-param::ly)<=eps) internalFlux[1]+=q;
                }
                // Positive vertical flow: entering at bottom, leaving at top.
                if (m==param::liquidSupplyMaterial || m==param::gasOutletMaterial) {
                    const T q=u[1]*velocityConversion*param::dx;
                    if (m==param::liquidSupplyMaterial) { inletVolume+=q; inletLiquid+=p*q; }
                    else { outletVolume+=q; outletLiquid+=p*q; }
                }
                if (m==param::openlbPoreMaterial && x[1]>=-param::dx*1e-7 && x[1]<=param::ly+param::dx*1e-7) {
                    sum+=std::clamp(p,T(0),T(1)); ++count;
                    const int row=static_cast<int>(std::llround(x[1]/param::dx));
                    auto& h=heightRows[row]; ++h.count; h.raw+=p; h.clamped+=std::clamp(p,T(0),T(1));
                    if (p>=0.5) updateHeight(liquidNode,std::clamp(x[1],T(0),T(param::ly)));
                    // Interpolate vertical phi=0.5 crossings only across two
                    // adjacent pore nodes; never bridge a solid particle.
                    const olb::LatticeR<2> above=r+olb::Vector<int,2>{0,1};
                    if (g.isInside(above) && x[1]+param::dx<=param::ly+param::dx*1e-7
                        && g.get(above)==param::openlbPoreMaterial) {
                        const T upper=ch.get(above).template getFieldComponent<olb::descriptors::STATISTIC>(0);
                        if (std::isfinite(upper) && p!=upper && (p-0.5)*(upper-0.5)<=0.) {
                            const T y=x[1]+param::dx*(0.5-p)/(upper-p);
                            updateHeight(crossing,std::clamp(y,T(0),T(param::ly)));
                        }
                    }
                }
            });
        }
        if (lo < -0.1 || hi > 1.1 || maxSpeed > 0.1)
            throw std::runtime_error("Numerical limits exceeded at step " + std::to_string(step)
                + ": phi=[" + std::to_string(lo) + "," + std::to_string(hi)
                + "], speed=" + std::to_string(maxSpeed));
        const T time=s.navierStokes.getUnitConverter().getPhysTime(step);
        if (step==0) { initialArea=liquidArea; initialRegions=regionalArea; }
        else {
            for (int k=0; k<2; ++k)
                cumulativeInternal[k]+=0.5*(previousInternal[k]+internalFlux[k])*dt;
            cumulativeIn+=0.5*(previousIn+inletLiquid)*dt;
            cumulativeOut+=0.5*(previousOut+outletLiquid)*dt;
        }
        previousInternal=internalFlux;
        previousIn=inletLiquid; previousOut=outletLiquid;
        diagnostics << step << ',' << time << ',' << liquidArea << ',' << liquidArea-initialArea
            << ',' << crossing << ',' << liquidNode << ',' << inletVolume << ',' << outletVolume
            << ',' << inletLiquid << ',' << outletLiquid << ',' << cumulativeIn << ',' << cumulativeOut
            << ',' << liquidArea-initialArea-cumulativeIn+cumulativeOut << '\n';
        std::array<T,3> changes{};
        for (int k=0; k<3; ++k) changes[k]=regionalArea[k]-initialRegions[k];
        regions << step << ',' << time;
        for (T area : regionalArea) regions << ',' << area;
        for (T change : changes) regions << ',' << change;
        regions << ',' << liquidArea-initialArea << ','
            << (initialArea!=0. ? (liquidArea-initialArea)/initialArea : std::numeric_limits<T>::quiet_NaN())
            << ',' << internalFlux[0] << ',' << internalFlux[1]
            << ',' << cumulativeInternal[0] << ',' << cumulativeInternal[1]
            << ',' << changes[0]-cumulativeIn+cumulativeInternal[0]
            << ',' << changes[1]-cumulativeInternal[0]+cumulativeInternal[1]
            << ',' << changes[2]-cumulativeInternal[1]+cumulativeOut << '\n';
        csv << step << ',' << time << ',' << (count ? sum/count : 0.) << ',' << lo << ',' << hi << ',' << maxSpeed << '\n';
        if (step%outputEvery==0 || step==steps) {
            for (const auto& [row,h] : heightRows)
                profiles << step << ',' << time << ',' << row*param::dx << ',' << h.count
                    << ',' << h.raw/h.count << ',' << h.clamped/h.count << '\n';
            csv.flush(); diagnostics.flush(); profiles.flush(); regions.flush(); writer.write(step);
            std::cout << "Step " << step << ", time=" << time << " s, electrode saturation="
                      << (count ? sum/count : 0.) << ", phi=[" << lo << ',' << hi << "]\n" << std::flush;
        }
    }
}
}
#endif
