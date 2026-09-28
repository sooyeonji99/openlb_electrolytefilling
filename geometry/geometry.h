#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "olb2D.h"
#include "olb2D.hh"

#include <vector>  // 이건 원래 c++ 표준 라이브러리

using namespace olb;

struct Particle2D
{
    double x;
    double y;
    double radius;
};


// Random particle generation
std::vector<Particle2D> generateParticles();


// Calculate analytical particle porosity
double calculatePorosity(
    const std::vector<Particle2D>& particles
);


// Calculate average particle radius
double calculateAverageRadius(
    const std::vector<Particle2D>& particles
);


// Build OpenLB geometry
void buildGeometry(
    SuperGeometry<double,2>& superGeometry,
    const std::vector<Particle2D>& particles
);


// Write geometry for ParaView
void writeGeometryVTK(
    SuperGeometry<double,2>& superGeometry
);

#endif