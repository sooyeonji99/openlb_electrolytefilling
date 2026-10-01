#include "parameters.h"
#include "geometry.h"
#include "vtk_writer.h"
#include "analysis.h"

#include <iostream>


int main()
{
    // =================================================
    // 1. Generate random cathode particles
    // =================================================

    auto particles = generateParticles();


    // =================================================
    // 2. Calculate analytical geometry properties
    // =================================================

    double analyticalPorosity =
        calculatePorosity(particles);

    double averageRadius =
        calculateAverageRadius(particles);


    // =================================================
    // 3. Define computational grid
    // =================================================

    int nx = static_cast<int>(
        Parameters::domainLengthX /
        Parameters::dx
    );

    int ny = static_cast<int>(
        Parameters::domainLengthY /
        Parameters::dx
    );


    // =================================================
    // 4. Convert particles to binary material grid
    // =================================================

    auto materialGrid =
        createMaterialGrid(
            particles,
            nx,
            ny,
            Parameters::dx
        );


    // =================================================
    // 5. Calculate grid-based porosity
    // =================================================

    double gridPorosity =
        calculateGridPorosity(materialGrid);

    bool verticallyConnected =
    checkVerticalPoreConnectivity(
        materialGrid,
        nx,
        ny
    );

    double connectedPoreFraction =
        calculateConnectedPoreFraction(
            materialGrid,
            nx,
            ny
        );

    
    auto poreDistanceMap =
        calculatePoreDistanceMap(
            materialGrid,
            nx,
            ny,
            Parameters::dx
        );


    double meanPoreDiameter =
        calculateMeanPoreDiameter(
            poreDistanceMap,
            materialGrid
        );


    double maxPoreDiameter =
        calculateMaxPoreDiameter(
            poreDistanceMap,
            materialGrid
        );
    // =================================================
    // 6. Print geometry information
    // =================================================

    std::cout << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "2D Cathode Random Geometry" << std::endl;
    std::cout << "======================================" << std::endl;

    std::cout
        << "Domain size = "
        << Parameters::domainLengthX * 1e6
        << " x "
        << Parameters::domainLengthY * 1e6
        << " um"
        << std::endl;

    std::cout
        << "Grid resolution = "
        << Parameters::dx * 1e6
        << " um"
        << std::endl;

    std::cout
        << "Grid size = "
        << nx
        << " x "
        << ny
        << std::endl;

    std::cout
        << "Number of particles = "
        << particles.size()
        << std::endl;

    std::cout
        << "Average particle radius = "
        << averageRadius * 1e6
        << " um"
        << std::endl;

    std::cout
        << "Analytical porosity = "
        << analyticalPorosity
        << std::endl;

    std::cout
        << "Grid porosity = "
        << gridPorosity
        << std::endl;

    std::cout << "======================================" << std::endl;


    std::cout
    << "Vertical pore connectivity = "
    << (verticallyConnected ? "YES" : "NO")
    << std::endl;

    std::cout
        << "Connected pore fraction = "
        << connectedPoreFraction
        << std::endl;

    // =================================================
    // 7. Write VTK file
    // =================================================
    std::cout
    << "Mean local pore diameter = "
    << meanPoreDiameter * 1e6
    << " um"
    << std::endl;

    std::cout
        << "Maximum local pore diameter = "
        << maxPoreDiameter * 1e6
        << " um"
        << std::endl;

    
    writeVTK(
        materialGrid,
        poreDistanceMap,
        nx,
        ny,
        Parameters::dx,
        "cathode_geometry.vtk"
    );

    writePoreSizeDistributionCSV(
        poreDistanceMap,
        materialGrid,
        0.5e-6,
        "pore_size_distribution.csv"
    );

    return 0;
}