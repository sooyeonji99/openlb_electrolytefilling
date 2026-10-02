#include "parameters.h"
#include "geometry.h"
#include "analysis.h"
#include "vtk_writer.h"

#include <iostream>
#include <iomanip>

int main()
{
    std::cout << std::setprecision(6);

    std::cout << "\n======================================\n";
    std::cout << "2D Cathode Random Porous Geometry\n";
    std::cout << "======================================\n";

    auto particles = generateParticles();

    std::cout << "Number of particles = "
              << particles.size()
              << "\n";

    auto material = buildMaterialGrid(particles);

    auto analysis = analyzePoreStructure(material);

    std::cout << "Porosity = "
              << analysis.porosity
              << "\n";

    std::cout << "Vertical pore connectivity = "
              << (analysis.verticalConnectivity ? "YES" : "NO")
              << "\n";

    std::cout << "Connected pore fraction = "
              << analysis.connectedPoreFraction
              << "\n";

    std::cout << "Mean pore radius = "
              << analysis.meanPoreRadius * 1.0e6
              << " um\n";

    std::cout << "Maximum pore radius = "
              << analysis.maxPoreRadius * 1.0e6
              << " um\n";

    writeGeometryVTK(material);

    writePoreRadiusDistributionCSV(
        analysis.poreRadiusValues
    );

    return 0;
}