
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>

// ============================================================
// Solid circular particle
// ============================================================

struct Particle2D
{
    double x;       // [m]
    double y;       // [m]
    double radius;  // [m]
};

// Generate random solid circular particles
std::vector<Particle2D> generateParticles();

// ============================================================
// Custom material grid
//
// 0 = pore (outside circles)
// 1 = solid (inside circles)
// ============================================================

std::vector<int> buildMaterialGrid(
    const std::vector<Particle2D>& particles
);

// ============================================================
// Circle radius statistics
// ============================================================

double calculateAverageParticleRadius(
    const std::vector<Particle2D>& particles
);

double calculateMinimumParticleRadius(
    const std::vector<Particle2D>& particles
);

double calculateMaximumParticleRadius(
    const std::vector<Particle2D>& particles
);

#endif
