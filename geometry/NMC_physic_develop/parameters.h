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

// Step 3: extend only the OpenLB domain below the unchanged electrode (y=0).
constexpr double reservoirHeight = 20.0e-6; // [m], reservoir y in [-height, 0)
constexpr double initialInterfaceY = 0.0; // [m], phi=0.5 at electrode entrance
constexpr double initialInterfaceWidth = 4.0 * dx; // [m], initial smoothing only
static_assert(reservoirHeight > 0.0);
static_assert(initialInterfaceWidth > 0.0);
static_assert(initialInterfaceY > -reservoirHeight && initialInterfaceY <= 0.0);

// Step 4 numerical verification properties; NOT measured electrolyte/air data.
constexpr double gasBufferHeight = 10.e-6;
constexpr double physTimeStep = 1.e-9; // [s]
constexpr double physLiquidDensity = 1000.; // [kg/m^3]
constexpr double physGasDensity = 100.; // reduced density ratio = 10
constexpr double physLiquidNu = 1.e-6; // [m^2/s]
constexpr double physGasNu = 1.e-6;
constexpr double physSurfaceTension = 1.e-3; // [N/m], verification only
constexpr double contactAngleDegrees = 60.; // liquid-side target
constexpr double phaseRelaxationTime = 1.;
constexpr int liquidSupplyMaterial = 3; // pressure boundary, fixed phi=1
constexpr int gasOutletMaterial = 4; // pressure boundary, fixed phi=0 gas bath
constexpr int sideWallMaterial = 5; // stationary wall, same wetting as particles

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
