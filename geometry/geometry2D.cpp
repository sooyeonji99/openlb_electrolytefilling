#include "geometry.h"
#include "parameters.h"

#include <random>
#include <cmath>
#include <iostream>

using namespace olb;


// ============================================================
// Generate random particles
// ============================================================

std::vector<Particle2D> generateParticles()
{
    std::vector<Particle2D> particles;

    std::mt19937 generator(param::randomSeed);

    std::uniform_real_distribution<double> radiusDistribution(
        param::radiusMin,
        param::radiusMax
    );

    std::uniform_real_distribution<double> xDistribution(
        param::radiusMax,
        param::lx - param::radiusMax
    );

    std::uniform_real_distribution<double> yDistribution(
        param::radiusMax,
        param::ly - param::radiusMax
    );


    for (int i = 0; i < param::numParticles; ++i)
    {
        Particle2D p;

        p.radius = radiusDistribution(generator);
        p.x = xDistribution(generator);
        p.y = yDistribution(generator);

        particles.push_back(p);
    }

    return particles;
}


// ============================================================
// Calculate average particle radius
// ============================================================

double calculateAverageRadius(
    const std::vector<Particle2D>& particles
)
{
    double sum = 0.0;

    for (const auto& p : particles)
    {
        sum += p.radius;
    }

    return sum / particles.size();
}


// ============================================================
// Calculate analytical porosity
// ============================================================

double calculatePorosity(
    const std::vector<Particle2D>& particles
)
{
    double solidArea = 0.0;

    for (const auto& p : particles)
    {
        solidArea += M_PI * p.radius * p.radius;
    }

    const double totalArea =
        param::lx * param::ly;

    const double porosity =
        1.0 - solidArea / totalArea;

    return porosity;
}


// ============================================================
// Build OpenLB geometry
// ============================================================

void buildGeometry(
    SuperGeometry<double,2>& superGeometry,
    const std::vector<Particle2D>& particles
)
{
    OstreamManager clout(
        std::cout,
        "buildGeometry"
    );

    clout << "Creating porous geometry..." << std::endl;


    // Entire domain initially = pore
    superGeometry.rename(
        0,
        param::poreMaterial
    );


    // Add particles
    for (const auto& p : particles)
    {
        Vector<double,2> center(
            p.x,
            p.y
        );

        IndicatorCircle2D<double> particle(
            center,
            p.radius
        );

        superGeometry.rename(
            param::poreMaterial,
            param::solidMaterial,
            particle
        );
    }


    superGeometry.clean();
    superGeometry.checkForErrors();

    superGeometry.getStatistics().print();


    clout << "Geometry generation complete."
          << std::endl;
}