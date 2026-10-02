#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>

struct Particle2D
{
    double x;
    double y;
    double radius;
};


// ------------------------------------------------------------
// Particle generation
// ------------------------------------------------------------

std::vector<Particle2D> generateParticles();


// ------------------------------------------------------------
// Geometry construction
// ------------------------------------------------------------

// 0 = pore
// 1 = solid
std::vector<int> buildMaterialGrid(
    const std::vector<Particle2D>& particles
);


// ------------------------------------------------------------
// Particle statistics
// ------------------------------------------------------------

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