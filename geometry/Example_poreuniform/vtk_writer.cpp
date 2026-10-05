#include "vtk_writer.h"

#include <cmath>
#include <fstream>
#include <iostream>


// =====================================================
// Create binary material grid
// =====================================================

std::vector<int> createMaterialGrid(
    const std::vector<Particle2D>& particles,
    int nx,
    int ny,
    double dx
)
{
    // Default: all cells are pore (= 0)
    std::vector<int> materialGrid(nx * ny, 0);

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            // Physical position of cell center
            double x = (i + 0.5) * dx;
            double y = (j + 0.5) * dx;

            // Check whether this cell lies inside a particle
            for (const auto& particle : particles)
            {
                double deltaX = x - particle.x;
                double deltaY = y - particle.y;

                double distanceSquared =
                    deltaX * deltaX +
                    deltaY * deltaY;

                if (distanceSquared <=
                    particle.radius * particle.radius)
                {
                    materialGrid[j * nx + i] = 1;
                    break;
                }
            }
        }
    }

    return materialGrid;
}


// =====================================================
// Calculate porosity from discretized geometry
// =====================================================

double calculateGridPorosity(
    const std::vector<int>& materialGrid
)
{
    int poreCells = 0;

    for (int material : materialGrid)
    {
        if (material == 0)
        {
            ++poreCells;
        }
    }

    return static_cast<double>(poreCells)
           / static_cast<double>(materialGrid.size());
}


// =====================================================
// Write VTK file
// =====================================================

void writeVTK(
    const std::vector<int>& materialGrid,
    const std::vector<double>& poreDistanceMap,
    int nx,
    int ny,
    double dx,
    const std::string& filename
)
{
    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cerr
            << "Error: could not open VTK file."
            << std::endl;

        return;
    }


    // =================================================
    // VTK header
    // =================================================

    file << "# vtk DataFile Version 3.0\n";
    file << "2D Cathode Random Porous Geometry\n";
    file << "ASCII\n";

    file << "DATASET STRUCTURED_POINTS\n";

    file << "DIMENSIONS "
         << nx + 1 << " "
         << ny + 1 << " "
         << 1 << "\n";

    file << "ORIGIN 0 0 0\n";

    file << "SPACING "
         << dx << " "
         << dx << " "
         << dx << "\n";


    // =================================================
    // Cell data
    // =================================================

    file << "CELL_DATA "
         << nx * ny
         << "\n";


    // =================================================
    // Material
    //
    // 0 = pore
    // 1 = NMC
    // =================================================

    file << "SCALARS material int 1\n";
    file << "LOOKUP_TABLE default\n";

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            int index =
                j * nx + i;

            file
                << materialGrid[index]
                << "\n";
        }
    }


    // =================================================
    // Local pore diameter [um]
    // =================================================

    file << "SCALARS poreDiameter_um double 1\n";
    file << "LOOKUP_TABLE default\n";

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            int index =
                j * nx + i;


            double poreDiameter = 0.0;


            if (materialGrid[index] == 0)
            {
                poreDiameter =
                    2.0 *
                    poreDistanceMap[index] *
                    1e6;
            }


            file
                << poreDiameter
                << "\n";
        }
    }


    file.close();


    std::cout
        << "VTK file written: "
        << filename
        << std::endl;
}