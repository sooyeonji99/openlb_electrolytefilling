#ifndef VTK_WRITER_H
#define VTK_WRITER_H

#include "geometry.h"

#include <string>
#include <vector>

// Convert particle geometry to a binary grid
// 0 = pore
// 1 = NMC solid
std::vector<int> createMaterialGrid(
    const std::vector<Particle2D>& particles,
    int nx,
    int ny,
    double dx
);

// Calculate porosity from the discretized grid
double calculateGridPorosity(
    const std::vector<int>& materialGrid
);

// Write the binary geometry as a legacy VTK file
void writeVTK(
    const std::vector<int>& materialGrid,
    int nx,
    int ny,
    double dx,
    const std::string& filename
);

#endif