
#include "olb2D.h"
#include "olb2D.hh"

#include "parameters.h"
#include "geometry.h"
#include "analysis.h"
#include "vtk_writer.h"
#include "openlb_geometry.h"

#include <iostream>
#include <iomanip>

using namespace olb;

int main(int argc, char* argv[])
{
    // ========================================================
    // Initialize OpenLB
    // ========================================================

    initialize(&argc, &argv);

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
        static_cast<T>(param::lx),
        static_cast<T>(param::ly)
    };

    Vector<T,2> origin {
        T(0),
        T(0)
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

    std::cout
        << "\n======================================\n"
        << "Geometry preprocessing completed.\n"
        << "======================================\n";

    return 0;
}
