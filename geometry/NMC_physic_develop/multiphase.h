#ifndef NMC_MULTIPHASE_H
#define NMC_MULTIPHASE_H

#include "olb2D.h"
#include "olb2D.hh"
#include "parameters.h"
#include "openlb_geometry.h"
#include <stdexcept>
#include <cmath>

namespace nmc {

using T = double;
using NSDescriptor = olb::descriptors::D2Q9<
    olb::descriptors::RHO, olb::descriptors::NABLARHO,
    olb::descriptors::FORCE, olb::descriptors::EXTERNAL_FORCE,
    olb::descriptors::TAU_EFF>;
using CHDescriptor = olb::descriptors::D2Q9<
    olb::descriptors::FORCE, olb::descriptors::SOURCE,
    olb::descriptors::SOURCE_OLD, olb::descriptors::PHIWETTING,
    olb::descriptors::VELOCITY, olb::descriptors::OLD_PHIU,
    olb::descriptors::STATISTIC, olb::descriptors::CHEM_POTENTIAL,
    olb::descriptors::BOUNDARY>;

using MixtureRules = olb::LinearTauViscosity;
using CouplingOperator = olb::WellBalancedCahnHilliardPostProcessor<MixtureRules>;
using CoupledLattices = olb::meta::map<
    olb::names::NavierStokes, olb::descriptors::VALUED_DESCRIPTOR<T, NSDescriptor>,
    olb::names::Component1, olb::descriptors::VALUED_DESCRIPTOR<T, CHDescriptor>>;
using Coupling = olb::SuperLatticeCoupling<CouplingOperator, CoupledLattices>;

// Reduced density ratio / surface tension for numerical verification only.
struct LatticeSetup {
    T tauLiquid = 0.5 + 3. * param::physLiquidNu * param::physTimeStep / (param::dx * param::dx);
    T tauGas = 0.5 + 3. * param::physGasNu * param::physTimeStep / (param::dx * param::dx);
    T tauPhase = param::phaseRelaxationTime;
    T rhoLiquid = 1.0;
    T rhoGas = param::physGasDensity / param::physLiquidDensity;
};

// Raster geometry uses half-way links consistently for NS and CH. This avoids
// inventing sub-grid circular surfaces different from the master material grid.
class HalfwayBoundary2D : public olb::IndicatorF2D<T> {
public:
    HalfwayBoundary2D() {
        this->_myMin = {-param::dx, -param::reservoirHeight-param::dx};
        this->_myMax = {param::lx+param::dx, param::ly+param::gasBufferHeight+param::dx};
    }
    bool operator()(bool out[], const T[]) override { out[0] = true; return true; }
    bool distance(T& dist, const olb::Vector<T,2>&,
                  const olb::Vector<T,2>& direction, int = -1) override {
        dist = 0.5 * std::sqrt(direction[0]*direction[0] + direction[1]*direction[1]);
        return true;
    }
};

// Geometry must outlive this object. Coupling is destroyed before its lattices.
// Initial populations may be set for visualization; initialize/apply/
// collideAndStream must wait until wetting and domain boundaries are ready.
class InitialPhase2D : public olb::AnalyticalF2D<T,T> {
    bool _statistic;
public:
    explicit InitialPhase2D(bool statistic = false)
        : olb::AnalyticalF2D<T,T>(statistic ? 2 : 1), _statistic(statistic) { }
    bool operator()(T output[], const T input[]) override {
        output[0] = 0.5 * (1. - std::tanh(
            2. * (input[1] - param::initialInterfaceY) / param::initialInterfaceWidth));
        if (_statistic) output[1] = 0.;
        return true;
    }
};

class MultiphaseSystem {
    LatticeSetup _setup;
public:
    olb::SuperLattice<T, NSDescriptor> navierStokes;
    olb::SuperLattice<T, CHDescriptor> cahnHilliard;
    Coupling coupling;

    explicit MultiphaseSystem(olb::SuperGeometry<T,2>& geometry,
                             const LatticeSetup& setup = {})
        : _setup(setup), navierStokes(geometry), cahnHilliard(geometry),
          coupling(CouplingOperator{},
                   olb::names::NavierStokes{}, navierStokes,
                   olb::names::Component1{}, cahnHilliard)
    {
        if (!(setup.tauLiquid > 0.5 && setup.tauGas > 0.5
              && setup.tauPhase > 0.5 && setup.rhoLiquid > 0.
              && setup.rhoGas > 0.)) {
            throw std::invalid_argument("Invalid multiphase lattice setup");
        }

        // OpenLB material 1 = pore; material 2 = solid circular particles.
        // Liquid/gas are represented by phi in the CH lattice, not material IDs.
        olb::dynamics::set<olb::NoDynamics>(navierStokes,
            geometry.getMaterialIndicator({0, param::openlbSolidMaterial, param::sideWallMaterial}));
        olb::dynamics::set<olb::NoDynamics>(cahnHilliard,
            geometry.getMaterialIndicator({0, param::openlbSolidMaterial, param::sideWallMaterial}));
        olb::dynamics::set<olb::MultiPhaseIncompressibleBGKdynamics<T, NSDescriptor>>(
            navierStokes, geometry.getMaterialIndicator({param::openlbPoreMaterial, param::liquidSupplyMaterial, param::gasOutletMaterial}));
        olb::dynamics::set<olb::WellBalancedCahnHilliardBGKdynamics<T, CHDescriptor>>(
            cahnHilliard, geometry.getMaterialIndicator({param::openlbPoreMaterial, param::liquidSupplyMaterial, param::gasOutletMaterial}));

        coupling.restrictTo(geometry.getMaterialIndicator({param::openlbPoreMaterial, param::liquidSupplyMaterial, param::gasOutletMaterial}));
        coupling.setParameter<MixtureRules::TAUS>({setup.tauGas, setup.tauLiquid});
        coupling.setParameter<MixtureRules::RHOS>({setup.rhoGas, setup.rhoLiquid});
        navierStokes.setParameter<olb::descriptors::OMEGA>(1. / setup.tauLiquid);
        cahnHilliard.setParameter<olb::descriptors::OMEGA>(1. / setup.tauPhase);

        cahnHilliard.addPostProcessor<olb::stage::PostStream>(
            geometry.getMaterialIndicator({param::openlbPoreMaterial, param::liquidSupplyMaterial, param::gasOutletMaterial}),
            olb::meta::id<olb::RhoWettingStatistics>());
        cahnHilliard.addPostProcessor<olb::stage::ChemPotCalc>(
            geometry.getMaterialIndicator(param::openlbPoreMaterial),
            olb::meta::id<olb::ChemPotentialPhaseFieldProcessor>());

        auto& communicator = cahnHilliard.getCommunicator(olb::stage::PreCoupling());
        communicator.requestOverlap(2);
        communicator.requestField<olb::descriptors::STATISTIC>();
        communicator.requestField<olb::descriptors::PHIWETTING>();
        communicator.requestField<olb::descriptors::CHEM_POTENTIAL>();
        communicator.exchangeRequests();
    }

    void setInitialPhases(olb::SuperGeometry<T,2>& geometry) {
        InitialPhase2D phi;
        InitialPhase2D phaseStatistic(true);
        olb::AnalyticalConst2D<T,T> zero(0.), two(2.);
        olb::AnalyticalConst2D<T,T> zeroVelocity(olb::Vector<T,2>(0., 0.));
        olb::AnalyticalConst2D<T,T> rhoGas(_setup.rhoGas);
        olb::AnalyticalConst2D<T,T> deltaRho(_setup.rhoLiquid - _setup.rhoGas);
        olb::AnalyticalConst2D<T,T> tauGas(_setup.tauGas);
        olb::AnalyticalConst2D<T,T> deltaTau(_setup.tauLiquid - _setup.tauGas);
        olb::AnalyticalIdentity2D<T,T> rho(rhoGas + deltaRho * phi);
        olb::AnalyticalIdentity2D<T,T> tau(tauGas + deltaTau * phi);

        auto pore = geometry.getMaterialIndicator({param::openlbPoreMaterial, param::liquidSupplyMaterial, param::gasOutletMaterial});
        auto solid = geometry.getMaterialIndicator({param::openlbSolidMaterial, param::sideWallMaterial});
        auto all = geometry.getMaterialIndicator({0, param::openlbPoreMaterial,
                                                 param::openlbSolidMaterial, param::liquidSupplyMaterial,
                                                 param::gasOutletMaterial, param::sideWallMaterial});
        // NS zeroth moment represents hydrodynamic pressure, NOT mixture density.
        navierStokes.defineRhoU(all, zero, zeroVelocity);
        navierStokes.iniEquilibrium(all, zero, zeroVelocity);
        cahnHilliard.defineRhoU(all, zero, zeroVelocity);
        cahnHilliard.iniEquilibrium(all, zero, zeroVelocity);
        cahnHilliard.defineRhoU(pore, phi, zeroVelocity);
        cahnHilliard.iniEquilibrium(pore, phi, zeroVelocity);

        navierStokes.defineField<olb::descriptors::RHO>(pore, rho);
        navierStokes.defineField<olb::descriptors::TAU_EFF>(pore, tau);
        navierStokes.defineField<olb::descriptors::FORCE>(all, zeroVelocity);
        navierStokes.defineField<olb::descriptors::EXTERNAL_FORCE>(all, zeroVelocity);
        navierStokes.defineField<olb::descriptors::BODY_FORCE>(all, zeroVelocity);
        navierStokes.defineField<olb::descriptors::NABLARHO>(all, zeroVelocity);
        cahnHilliard.defineField<olb::descriptors::STATISTIC>(all, zeroVelocity);
        cahnHilliard.defineField<olb::descriptors::STATISTIC>(pore, phaseStatistic);
        cahnHilliard.defineField<olb::descriptors::PHIWETTING>(pore, phi);
        // Exterior/solid stencil neighbors are boundaries; only fluid is bulk.
        cahnHilliard.defineField<olb::descriptors::BOUNDARY>(all, two);
        cahnHilliard.defineField<olb::descriptors::BOUNDARY>(pore, zero);
        cahnHilliard.defineField<olb::descriptors::BOUNDARY>(solid, two);
        cahnHilliard.defineField<olb::descriptors::CHEM_POTENTIAL>(all, zero);
        cahnHilliard.defineField<olb::descriptors::SOURCE>(all, zero);
        cahnHilliard.defineField<olb::descriptors::SOURCE_OLD>(all, zero);
        cahnHilliard.defineField<olb::descriptors::VELOCITY>(all, zeroVelocity);
        cahnHilliard.defineField<olb::descriptors::OLD_PHIU>(all, zeroVelocity);
        cahnHilliard.defineField<olb::descriptors::FORCE>(all, zeroVelocity);
        // Wall setup and capillary initialization are applied next.
    }

    void prepareCapillaryPhysics(olb::SuperGeometry<T,2>& geometry) {
        olb::OstreamManager out(std::cout, "CapillarySetup");
        out << "Setting converter and capillary parameters" << std::endl;
        // The ordinary UnitConverter leaves surface tension conversion unset.
        // Use the multiphase converter confirmed in bubbleChannel2d.
        navierStokes.setUnitConverter<olb::MultiPhaseUnitConverterFromRelaxationTime<T,NSDescriptor>>(
            param::nx, _setup.tauLiquid, _setup.rhoLiquid, param::lx,
            param::physLiquidNu, param::physLiquidDensity);
        cahnHilliard.setUnitConverter(navierStokes.getUnitConverter());
        const auto& converter = navierStokes.getUnitConverter();
        const T sigma = param::physSurfaceTension / converter.getConversionFactorSurfaceTension();
        const T width = param::initialInterfaceWidth / param::dx;

        HalfwayBoundary2D halfWay;
        out << "Setting half-way wall links" << std::endl;
        auto walls = geometry.getMaterialIndicator({param::openlbSolidMaterial, param::sideWallMaterial});
        auto fluid = geometry.getMaterialIndicator({param::openlbPoreMaterial, param::liquidSupplyMaterial, param::gasOutletMaterial});
        olb::setBouzidiBoundary(navierStokes, walls, fluid, halfWay);
        olb::setBouzidiBoundary(cahnHilliard, walls, fluid, halfWay);

        // The stock setBouzidiWellBalanced visits solid interiors and throws
        // for zero normals. Register its actual wetting processor only on skin
        // cells, choosing an existing fluid neighbor for ambiguous raster corners.
        std::size_t skinCount = 0, fallbackCount = 0;
        auto& load = geometry.getLoadBalancer();
        for (int iC=0; iC<load.size(); ++iC) {
            auto& g = geometry.getBlockGeometry(iC);
            auto& ch = cahnHilliard.getBlock(iC);
            g.forSpatialLocations([&](olb::LatticeR<2> r) {
                if (g.get(r) != param::openlbSolidMaterial && g.get(r) != param::sideWallMaterial) return;
                auto isFluid = [&](olb::LatticeR<2> n) {
                    if (!g.isInside(n)) return false;
                    const int m = g.get(n);
                    return m == param::openlbPoreMaterial || m == param::liquidSupplyMaterial || m == param::gasOutletMaterial;
                };
                olb::Vector<int,2> toward{};
                bool skin = false;
                for (int k=1; k<CHDescriptor::q; ++k) {
                    const auto c = olb::descriptors::c<CHDescriptor>(k);
                    if (isFluid(r + c)) { toward = c; skin = true; break; }
                }
                if (!skin || ch.isPadding(r)) return;
                const auto type = g.getStatistics().getType(r[0], r[1]);
                olb::Vector<int,2> normal{type[1], type[2]};
                if ((normal[0] == 0 && normal[1] == 0) || !isFluid(r - normal)) {
                    normal = toward * -1;
                    ++fallbackCount;
                }
                ch.addPostProcessor(typeid(olb::stage::PreCoupling), r,
                    olb::boundaryhelper::promisePostProcessorForNormal<T,CHDescriptor,olb::WellBalancedWallProcessor2D>(normal));
                ++skinCount;
            });
        }

        auto open = geometry.getMaterialIndicator({param::liquidSupplyMaterial, param::gasOutletMaterial});
        out << "Setting supply/outlet pressure and phase conditions" << std::endl;
        olb::boundary::set<olb::boundary::IncompressibleZouHePressure>(navierStokes, open);
        using PhaseDirichlet = olb::boundary::RegularizedTemperature<T,CHDescriptor,
            olb::WellBalancedCahnHilliardBGKdynamics<T,CHDescriptor>>;
        olb::boundary::set<PhaseDirichlet>(cahnHilliard, open);
        // setParameter updates currently registered operators only in this API.
        // Set all values AFTER registering wall and inlet/outlet operators.
        navierStokes.setParameter<olb::descriptors::OMEGA>(1. / _setup.tauLiquid);
        cahnHilliard.setParameter<olb::descriptors::OMEGA>(1. / _setup.tauPhase);
        cahnHilliard.setParameter<olb::descriptors::SCALAR>(sigma);
        cahnHilliard.setParameter<olb::descriptors::INTERFACE_WIDTH>(width);
        cahnHilliard.setParameter<olb::descriptors::THETA>(
            M_PI - param::contactAngleDegrees * M_PI / 180.);
        cahnHilliard.setParameter<olb::descriptors::ADDEND>(0.);
        olb::AnalyticalConst2D<T,T> zero(0.), liquid(1.), gas(0.);
        navierStokes.defineRho(open, zero); // equal hydrodynamic pressure baths
        cahnHilliard.defineRho(geometry.getMaterialIndicator(param::liquidSupplyMaterial), liquid);
        cahnHilliard.defineRho(geometry.getMaterialIndicator(param::gasOutletMaterial), gas);

        out << "Initializing lattices and capillary fields" << std::endl;
        navierStokes.initialize();
        cahnHilliard.initialize();
        // initialize() executes boundary postprocessors. Restore the prescribed
        // initial equilibria before evaluating forces (as contactAngle2d does),
        // so initial bounce-back does not create artificial phase overshoots.
        setInitialPhases(geometry);
        cahnHilliard.getCommunicator(olb::stage::PreCoupling()).communicate();
        cahnHilliard.executePostProcessors(olb::stage::PreCoupling());
        cahnHilliard.getCommunicator(olb::stage::PreCoupling()).communicate();
        cahnHilliard.executePostProcessors(olb::stage::ChemPotCalc());
        cahnHilliard.getCommunicator(olb::stage::PreCoupling()).communicate();
        coupling.apply();
        out << "Verification only: rho ratio=" << _setup.rhoLiquid / _setup.rhoGas
            << ", tau_l=" << _setup.tauLiquid << ", tau_g=" << _setup.tauGas
            << ", sigma_lattice=" << sigma << ", interface_width=" << width
            << ", target liquid contact angle=" << param::contactAngleDegrees
            << ", wall skin cells=" << skinCount << ", ambiguous normals=" << fallbackCount << std::endl;
    }
};

class PhaseOutput2D : public olb::SuperF2D<T,T> {
    olb::SuperLattice<T,CHDescriptor>& _lattice;
public:
    explicit PhaseOutput2D(olb::SuperLattice<T,CHDescriptor>& lattice)
        : olb::SuperF2D<T,T>(lattice, 1), _lattice(lattice) {
        this->getName() = "phi";
    }
    bool operator()(T output[], const int input[]) override {
        auto& load = _lattice.getLoadBalancer();
        if (!load.isLocal(input[0])) return false;
        output[0] = _lattice.getBlock(load.loc(input[0]))
            .get(input[1], input[2])
            .template getFieldComponent<olb::descriptors::STATISTIC>(0);
        return true;
    }
};

inline void writeInitialPhases(MultiphaseSystem& system,
                               olb::SuperGeometry<T,2>& geometry) {
    // NoDynamics has a constant density moment on solids. Output the phase
    // statistic as contactAngle2d does, rather than treating that moment as phi.
    // STATISTIC has two components; only component 0 is the phase variable.
    PhaseOutput2D phi(system.cahnHilliard);
    olb::SuperLatticeExternalScalarField2D<T,NSDescriptor,olb::descriptors::RHO>
        rho(system.navierStokes);
    rho.getName() = "rho_lattice";
    olb::SuperGeometryMaterialOutput2D<T> material(geometry);
    olb::SuperVTMwriter2D<T> writer("initial_phases");
    writer.addFunctor(material);
    writer.addFunctor(phi);
    writer.addFunctor(rho);
    writer.createMasterFile();
    writer.write(0);
}

inline void writeCapillarySetup(MultiphaseSystem& system, olb::SuperGeometry<T,2>& geometry) {
    PhaseOutput2D phi(system.cahnHilliard);
    olb::SuperGeometryMaterialOutput2D<T> material(geometry);
    olb::SuperLatticeExternalScalarField2D<T,CHDescriptor,olb::descriptors::CHEM_POTENTIAL> mu(system.cahnHilliard);
    mu.getName() = "chemical_potential";
    olb::SuperLatticeField2D<T,NSDescriptor,olb::descriptors::FORCE> force(system.navierStokes);
    force.getName() = "acceleration_lattice";
    olb::SuperVTMwriter2D<T> writer("capillary_setup");
    writer.addFunctor(material);
    writer.addFunctor(phi);
    writer.addFunctor(mu);
    writer.addFunctor(force);
    writer.createMasterFile();
    writer.write(0);
}

} // namespace nmc
#endif
