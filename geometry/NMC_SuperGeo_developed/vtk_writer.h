#ifndef VTK_WRITER_H
#define VTK_WRITER_H

#include "geometry.h"

#include <vector>


void writeGeometryVTK(
    const std::vector<int>& material
);

void writeCircleRadiusDistributionCSV(
    const std::vector<Particle2D>& particles
);

#endif