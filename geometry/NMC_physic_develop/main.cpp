
#include "olb2D.h"
#include "olb2D.hh"

#include "parameters.h"
#include "geometry.h"
#include "analysis.h"
#include "vtk_writer.h"
#include "openlb_geometry.h"
#include "multiphase.h"
#include "simulation.h"

#include <iostream>
#include <iomanip>

using namespace olb;

int main(int argc, char* argv[])
{
    // ========================================================
    // Initialize OpenLB
    // ========================================================

    initialize(&argc, &argv);

    int steps = 100, outputEvery = 10;
    for (int i=1; i<argc; ++i) {
        const std::string arg(argv[i]);
        if ((arg == "--steps" || arg == "--output-every") && i+1<argc) {
            const int value=std::stoi(argv[++i]);
            if (arg == "--steps") steps=value; else outputEvery=value;
        } else {
            std::cerr << "Usage: ./main [--steps N] [--output-every N]\n";
            return 1;
        }
    }
    if (steps<0 || outputEvery<1) return 1;

    using T = double;

    std::cout << std::setprecision(6);

    std::cout
        << "\n======================================\n"
        << "2D Cathode Random Porous Geometry\n"
        << "======================================\n";

    // ========================================================
    // Parameters
    // ========================================================

    std::cout
        << "Domain = "
        << param::lx * 1.0e6
        << " x "
        << param::ly * 1.0e6
        << " um\n"

        << "Resolution = "
        << param::dx * 1.0e6
        << " um\n"

        << "Grid = "
        << param::nx
        << " x "
        << param::ny
        << "\n"

        << "Target porosity = "
        << param::targetPorosity
        << "\n"

        << "Circle radius range = "
        << param::minParticleRadius * 1.0e6
        << " - "
        << param::maxParticleRadius * 1.0e6
        << " um\n";

    // ========================================================
    // 1. Generate solid circular particles
    // ========================================================

    auto particles = generateParticles();

    // ========================================================
    // 2. Circle radius statistics
    // ========================================================

    const double meanCircleRadius =
        calculateAverageParticleRadius(particles);

    const double minCircleRadius =
        calculateMinimumParticleRadius(particles);

    const double maxCircleRadius =
        calculateMaximumParticleRadius(particles);

    std::cout
        << "\n======================================\n"
        << "Solid Circle Statistics\n"
        << "======================================\n"

        << "Number of circles = "
        << particles.size()
        << "\n"

        << "Mean circle radius = "
        << meanCircleRadius * 1.0e6
        << " um\n"

        << "Minimum circle radius = "
        << minCircleRadius * 1.0e6
        << " um\n"

        << "Maximum circle radius = "
        << maxCircleRadius * 1.0e6
        << " um\n";

    // ========================================================
    // 3. Build custom material grid
    //
    // 0 = pore (outside circles)
    // 1 = solid (inside circles)
    // ========================================================

    auto material = buildMaterialGrid(particles);

    // ========================================================
    // 4. Analyze pore structure
    // ========================================================

    auto analysis = analyzePoreStructure(material);

    std::cout
        << "\n======================================\n"
        << "Pore Structure Analysis\n"
        << "======================================\n"

        << "Grid porosity = "
        << analysis.porosity
        << "\n"

        << "Vertical pore connectivity = "
        << (analysis.verticalConnectivity ? "YES" : "NO")
        << "\n"

        << "Bottom-connected pore fraction = "
        << analysis.connectedPoreFraction
        << "\n";

    // ========================================================
    // 5. VTK / CSV output
    // ========================================================

    writeGeometryVTK(material);

    writeCircleRadiusDistributionCSV(particles);

    // ========================================================
    // 6. Create OpenLB computational domain
    // ========================================================

    std::cout
        << "\n======================================\n"
        << "Creating OpenLB SuperGeometry\n"
        << "======================================\n";

    Vector<T,2> extend {
        static_cast<T>(param::lx + 2. * param::dx),
        static_cast<T>(param::ly + param::gasBufferHeight + param::reservoirHeight + 2. * param::dx)
    };

    Vector<T,2> origin {
        -static_cast<T>(param::dx),
        -static_cast<T>(param::reservoirHeight + param::dx)
    };

    IndicatorCuboid2D<T> cuboid(
        extend,
        origin
    );

    const int noOfCuboids = 1;

    CuboidDecomposition2D<T> cuboidDecomposition(
        cuboid,
        static_cast<T>(param::dx),
        noOfCuboids
    );

    HeuristicLoadBalancer<T> loadBalancer(
        cuboidDecomposition
    );

    SuperGeometry<T,2> superGeometry(
        cuboidDecomposition,
        loadBalancer
    );

    // ========================================================
    // 7. Map custom grid to OpenLB SuperGeometry
    //
    // OpenLB:
    // 0 = outside
    // 1 = pore / fluid
    // 2 = solid particle
    // ========================================================

    prepareOpenLBGeometry(
        superGeometry,
        material
    );

    // ========================================================
    // 8. OpenLB-native VTM output
    // ========================================================

    writeOpenLBGeometryVTK(
        superGeometry
    );

    // Step 2: construct NS/CH lattices and coupling on pore material 1.
    // Step 3 sets initial populations; boundary conditions still precede stepping.
    std::cout << "Constructing multiphase lattices" << std::endl;
    nmc::MultiphaseSystem multiphase(superGeometry);
    std::cout << "Setting initial phases" << std::endl;
    multiphase.setInitialPhases(superGeometry);
    std::cout << "Writing initial phases" << std::endl;
    nmc::writeInitialPhases(multiphase, superGeometry);
    multiphase.prepareCapillaryPhysics(superGeometry);
    nmc::writeCapillarySetup(multiphase, superGeometry);
    try {
        nmc::simulate(multiphase, superGeometry, steps, outputEvery);
    } catch (const std::exception& e) {
        std::cerr << "Simulation stopped: " << e.what() << std::endl;
        return 1;
    }

    std::cout
        << "\n======================================\n"
        << "Geometry, multiphase setup and initial phase output completed.\n"
        << "NS + CH bulk dynamics and coupling: material 1 (pore).\n"
        << "Material 2: solid; no fluid bulk dynamics.\n"
        << "Reservoir height = " << param::reservoirHeight * 1.e6 << " um\n"
        << "Initial interface y = " << param::initialInterfaceY * 1.e6 << " um\n"
        << "Step 5: time integration completed.\n"
        << "======================================\n";

    return 0;
}
