// electrolyte의 경우 surface tension, contact angle, viscosity 값 필요


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
constexpr double lx = 100.0e-6;  // [m]
constexpr double ly = 100.0e-6;  // [m]

// Spatial resolution
constexpr double dx = 0.1e-6;    // [m]

// Number of grid cells
constexpr int nx =
    static_cast<int>(lx / dx);

constexpr int ny =
    static_cast<int>(ly / dx);


// ============================================================
// Particle / circle size distribution
// ============================================================

// Circle = solid electrode particle
constexpr double minParticleRadius =
    1.0e-7;                      // [m]

constexpr double maxParticleRadius =
    10.0e-6;                     // [m]

// Porosity = area outside circles / total area
constexpr double targetPorosity =
    0.30;

constexpr double porosityTolerance =
    0.002;


// ============================================================
// Random generation
// ============================================================

constexpr std::uint32_t randomSeed =
    12345;

constexpr int maxPlacementAttempts =
    2000000;

constexpr double minimumParticleGap =
    0.0;


// ============================================================
// Custom grid material IDs
//
// Circle interior = solid
// Circle exterior = pore
//
// These are NOT OpenLB material numbers.
// ============================================================

constexpr int poreMaterial =
    0;

constexpr int solidMaterial =
    1;


// ============================================================
// OpenLB material IDs
//
// 0 = outside / unassigned
// 1 = pore / fluid
// 2 = solid electrode particle
// ============================================================

constexpr int openlbPoreMaterial =
    1;

constexpr int openlbSolidMaterial =
    2;


// ============================================================
// Connectivity
// ============================================================

// 4-neighbor connectivity
constexpr bool useDiagonalConnectivity =
    false;


// ============================================================
// Output
// ============================================================

const std::string vtkFileName =
    "cathode_geometry.vtk";

// Distribution of generated solid-circle radii
const std::string circleRadiusCsvFileName =
    "circle_radius_distribution.csv";

}

#endif
