#ifndef PARAMETERS_H
#define PARAMETERS_H

namespace Parameters
{

// =====================================================
// Computational domain
// =====================================================

// Physical domain size [m]
constexpr double domainLengthX = 100e-6;
constexpr double domainLengthY = 100e-6;

// Lattice resolution [m]
constexpr double dx = 0.5e-6;
constexpr double poreSizeBinWidth = 0.5e-6;

// =====================================================
// Cathode particle parameters
// =====================================================

// Particle radius range [m]
constexpr double minParticleRadius = 2.0e-6;
constexpr double maxParticleRadius = 5.0e-6;

// Target pore volume fraction
// In 2D this corresponds to pore area fraction.
constexpr double targetPorosity = 0.30;


// =====================================================
// Random generation
// =====================================================

constexpr unsigned int randomSeed = 21;

// Maximum attempts for random particle placement
constexpr int maxPlacementAttempts = 100000;


// =====================================================
// Material numbers for OpenLB
// =====================================================

constexpr int poreMaterial = 1;
constexpr int solidMaterial = 2;

}

#endif