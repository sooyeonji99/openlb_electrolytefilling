#include "parameters.h"
#include "geometry.h"
#include "analysis.h"
#include "vtk_writer.h"
#include "openlb_geometry.h"

#include "olb2D.h"
// #include "olb2D.hh"

#include <iostream>
#include <iomanip>


using namespace olb;


int main(
    int argc,
    char* argv[])
{
    // ========================================================
    // Initialize OpenLB
    // ========================================================

    initialize(
        &argc,
        &argv
    );


    using T =
        double;


    std::cout
        << std::setprecision(6);


    std::cout
        << "\n======================================\n";

    std::cout
        << "2D Cathode Random Porous Geometry\n";

    std::cout
        << "======================================\n";


    // ========================================================
    // Parameters
    // ========================================================

    std::cout
        << "Domain = "
        << param::lx * 1.0e6
        << " x "
        << param::ly * 1.0e6
        << " um\n";


    std::cout
        << "Resolution = "
        << param::dx * 1.0e6
        << " um\n";


    std::cout
        << "Grid = "
        << param::nx
        << " x "
        << param::ny
        << "\n";


    std::cout
        << "Target porosity = "
        << param::targetPorosity
        << "\n";


    std::cout
        << "Particle radius = "
        << param::minParticleRadius
           * 1.0e6
        << " - "
        << param::maxParticleRadius
           * 1.0e6
        << " um\n";


    // ========================================================
    // 1. Generate particles
    // ========================================================

    auto particles =
        generateParticles();


    // ========================================================
    // 2. Particle statistics
    // ========================================================

    const double meanParticleRadius =
        calculateAverageParticleRadius(
            particles
        );


    const double minParticleRadius =
        calculateMinimumParticleRadius(
            particles
        );


    const double maxParticleRadius =
        calculateMaximumParticleRadius(
            particles
        );


    std::cout
        << "\n======================================\n";

    std::cout
        << "Particle Statistics\n";

    std::cout
        << "======================================\n";


    std::cout
        << "Number of particles = "
        << particles.size()
        << "\n";


    std::cout
        << "Mean particle radius = "
        << meanParticleRadius
           * 1.0e6
        << " um\n";


    std::cout
        << "Minimum particle radius = "
        << minParticleRadius
           * 1.0e6
        << " um\n";


    std::cout
        << "Maximum particle radius = "
        << maxParticleRadius
           * 1.0e6
        << " um\n";


    // ========================================================
    // 3. Build original material grid
    // ========================================================

    auto material =
        buildMaterialGrid(
            particles
        );


    // ========================================================
    // 4. Analyze original geometry
    // ========================================================

    auto analysis =
        analyzePoreStructure(
            material
        );


    std::cout
        << "\n======================================\n";

    std::cout
        << "Pore Structure Analysis\n";

    std::cout
        << "======================================\n";


    std::cout
        << "Grid porosity = "
        << analysis.porosity
        << "\n";


    std::cout
        << "Vertical pore connectivity = "
        << (
            analysis.verticalConnectivity
            ?
            "YES"
            :
            "NO"
        )
        << "\n";


    std::cout
        << "Connected pore fraction = "
        << analysis.connectedPoreFraction
        << "\n";


    std::cout
        << "Mean pore radius = "
        << analysis.meanPoreRadius
           * 1.0e6
        << " um\n";


    std::cout
        << "Maximum pore radius = "
        << analysis.maxPoreRadius
           * 1.0e6
        << " um\n";


    // ========================================================
    // 5. Existing VTK / CSV output
    // ========================================================

    writeGeometryVTK(
        material
    );


    writePoreRadiusDistributionCSV(
        analysis.poreRadiusValues
    );


    // ========================================================
    // 6. Create OpenLB computational domain
    // ========================================================

    std::cout
        << "\n======================================\n";

    std::cout
        << "Creating OpenLB SuperGeometry\n";

    std::cout
        << "======================================\n";


    Vector<T,2> extend {
        static_cast<T>(
            param::lx
        ),

        static_cast<T>(
            param::ly
        )
    };


    Vector<T,2> origin {
        T(0),
        T(0)
    };


    IndicatorCuboid2D<T> cuboid(
        extend,
        origin
    );


    // Single CPU / single cuboid for now
    const int noOfCuboids =
        1;


    CuboidDecomposition2D<T>
        cuboidDecomposition(
            cuboid,
            static_cast<T>(
                param::dx
            ),
            noOfCuboids
        );


    HeuristicLoadBalancer<T>
        loadBalancer(
            cuboidDecomposition
        );


    SuperGeometry<T,2>
        superGeometry(
            cuboidDecomposition,
            loadBalancer
        );


    // ========================================================
    // 7. Particle2D -> OpenLB material geometry
    // ========================================================

    prepareOpenLBGeometry(
        superGeometry,
        material
    );


    // ========================================================
    // 8. Print OpenLB geometry information
    // ========================================================

    // std::cout
    //     << "\n======================================\n";

    // std::cout
    //     << "OpenLB Geometry Statistics\n";

    // std::cout
    //     << "======================================\n";


    // superGeometry.print();


    std::cout
        << "\n======================================\n";

    std::cout
        << "Geometry preprocessing completed.\n";

    std::cout
        << "======================================\n";


    return 0;
}