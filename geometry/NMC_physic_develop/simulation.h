#ifndef NMC_SIMULATION_H
#define NMC_SIMULATION_H
#include "multiphase.h"
#include <fstream>
#include <iomanip>
#include <limits>

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
                if (m==param::openlbPoreMaterial && x[1]>=-param::dx*1e-7 && x[1]<=param::ly+param::dx*1e-7) {
                    sum+=std::clamp(p,T(0),T(1)); ++count;
                }
            });
        }
        if (lo < -0.1 || hi > 1.1 || maxSpeed > 0.1)
            throw std::runtime_error("Numerical limits exceeded at step " + std::to_string(step)
                + ": phi=[" + std::to_string(lo) + "," + std::to_string(hi)
                + "], speed=" + std::to_string(maxSpeed));
        const T time=s.navierStokes.getUnitConverter().getPhysTime(step);
        csv << step << ',' << time << ',' << (count ? sum/count : 0.) << ',' << lo << ',' << hi << ',' << maxSpeed << '\n';
        if (step%outputEvery==0 || step==steps) {
            csv.flush(); writer.write(step);
            std::cout << "Step " << step << ", time=" << time << " s, electrode saturation="
                      << (count ? sum/count : 0.) << ", phi=[" << lo << ',' << hi << "]\n" << std::flush;
        }
    }
}
}
#endif
