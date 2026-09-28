#ifndef PARAMETERS_H
#define PARAMETERS_H

namespace param {

// ============================================================
// Domain
// ============================================================

// Physical domain size [m], 계산 영역 정의
constexpr double lx = 100.0e-6;   // 100 um, constexpr: compile-time constant를 정의할 때 사용
constexpr double ly = 100.0e-6;   // 100 um

// Lattice resolution
constexpr double dx = 1.0e-6;     // 1 um


// ============================================================
// Particle generation
// ============================================================

// Number of particles
constexpr int numParticles = 100;

// Particle radius range [m]
constexpr double radiusMin = 2.0e-6;
constexpr double radiusMax = 5.0e-6;

// Random seed
constexpr unsigned int randomSeed = 12345; // random seed는 숫자 무조건 고정 필요


// ============================================================
// Material numbers
// ============================================================

// pore / fluid region
constexpr int poreMaterial = 1;

// solid electrode particle
constexpr int solidMaterial = 2;


// ============================================================
// Output
// ============================================================

constexpr const char* vtkName = "porousGeometry";

}

#endif