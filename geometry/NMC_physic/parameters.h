#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <cstdint>

namespace param
{

// ============================================================
// Physical domain
// ============================================================

// Electrode dimensions
constexpr double lx = 100.0e-6;   // [m]
constexpr double ly = 100.0e-6;   // [m]

// Spatial resolution
//
// 0.1 um resolution:
// 100 um / 0.1 um = 1000 cells
//
constexpr double dx = 0.1e-6;     // [m]

// Number of grid cells
constexpr int nx =
    static_cast<int>(lx / dx);

constexpr int ny =
    static_cast<int>(ly / dx);


// ============================================================
// Particle size distribution
// ============================================================

// Particle radius range
constexpr double minParticleRadius =
    1.0e-7;                       // [m] = 0.1 um

constexpr double maxParticleRadius =
    10.0e-6;                      // [m] = 10 um

// Target pore fraction
constexpr double targetPorosity =
    0.30;

// Allowed deviation during particle generation
constexpr double porosityTolerance =
    0.002;


// ============================================================
// Random generation
// ============================================================

constexpr std::uint32_t randomSeed =
    12345;

// Maximum consecutive failed placement attempts
constexpr int maxPlacementAttempts =
    2000000;

// Minimum gap between particles
//
// 0 = particles may touch
//     but cannot overlap
//
constexpr double minimumParticleGap =
    0.0;


// ============================================================
// Custom grid material IDs
//
// NOTE:
// These are NOT OpenLB material numbers.
// ============================================================

constexpr int poreMaterial =
    0;

constexpr int solidMaterial =
    1;


// ============================================================
// OpenLB material IDs
// ============================================================

constexpr int openlbSolidMaterial =
    1;

constexpr int openlbPoreMaterial =
    2;


// ============================================================
// Connectivity
// ============================================================

// Current connectivity analysis:
// 4-neighbor connectivity
constexpr bool useDiagonalConnectivity =
    false;


// ============================================================
// Output
// ============================================================

const std::string vtkFileName =
    "cathode_geometry.vtk";

const std::string poreRadiusCsvFileName =
    "pore_radius_distribution.csv";

}

#endif