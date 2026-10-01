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
// Write legacy VTK file
// =====================================================

void writeVTK(
    const std::vector<int>& materialGrid,
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

    file << "CELL_DATA "
         << nx * ny << "\n";

    file << "SCALARS material int 1\n";
    file << "LOOKUP_TABLE default\n";

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            file << materialGrid[j * nx + i]
                 << "\n";
        }
    }

    file.close();

    std::cout
        << "VTK file written: "
        << filename
        << std::endl;
}