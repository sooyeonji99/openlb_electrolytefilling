#ifndef GEOMETRY_H
#define GEOMETRY_H

// #include "olb2D.h"
// #include "olb2D.hh"

#include <vector>

// using namespace olb;


// =====================================================
// Particle structure
// =====================================================

struct Particle2D
{
    double x;
    double y;
    double radius;
};


// =====================================================
// Function declarations
// =====================================================

// Generate random cathode particles
std::vector<Particle2D> generateParticles();


// Calculate analytical porosity
double calculatePorosity(
    const std::vector<Particle2D>& particles
);


// Calculate average particle radius
double calculateAverageRadius(
    const std::vector<Particle2D>& particles
);


// Build OpenLB material geometry
// void buildGeometry(
//     SuperGeometry<double, 2>& superGeometry,
//     const std::vector<Particle2D>& particles
// );


// // Write geometry for ParaView
// void writeGeometryVTK(
//     SuperGeometry<double, 2>& superGeometry
// );


#endif