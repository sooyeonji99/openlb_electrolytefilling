#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <cstdint>

namespace param
{

// ============================================================
// Domain
// ============================================================

// Physical domain size
constexpr double lx = 100.0e-6;   // [m]
constexpr double ly = 100.0e-6;   // [m]

// Grid resolution
constexpr double dx = 0.1e-6;     // [m]

// Number of grid cells
constexpr int nx = static_cast<int>(lx / dx);
constexpr int ny = static_cast<int>(ly / dx);


// ============================================================
// Particle size distribution
// ============================================================

// Particle radius range
constexpr double minParticleRadius = 1.0e-7;   // [m] = 0.1 um
constexpr double maxParticleRadius = 10.0e-6;  // [m] = 10 um

// Target porosity
constexpr double targetPorosity = 0.30;

// Tolerance for terminating particle generation
constexpr double porosityTolerance = 0.002;


// ============================================================
// Random particle generation
// ============================================================

constexpr std::uint32_t randomSeed = 12345;

// Maximum number of placement trials
constexpr int maxPlacementAttempts = 2000000;

// Minimum gap between particles
// 0.0 = particles may touch but cannot overlap
constexpr double minimumParticleGap = 0.0;


// ============================================================
// Material IDs
// ============================================================

constexpr int poreMaterial  = 0;
constexpr int solidMaterial = 1;


// ============================================================
// Connectivity
// ============================================================

// 4-neighbor connectivity in 2D
constexpr bool useDiagonalConnectivity = false;


// ============================================================
// Pore-size analysis
// ============================================================

// Number of bins for pore-radius histogram
constexpr int poreRadiusBins = 100;


// ============================================================
// Output
// ============================================================

const std::string vtkFileName =
    "cathode_geometry.vtk";

const std::string poreRadiusCsvFileName =
    "pore_radius_distribution.csv";

}

#endif