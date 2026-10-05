#include "geometry.h"
#include "parameters.h"

#include <cmath>
#include <iostream>
#include <random>


// =====================================================
// Generate random particles
// =====================================================

std::vector<Particle2D> generateParticles()
{
    std::vector<Particle2D> particles;

    std::mt19937 generator(Parameters::randomSeed);

    std::uniform_real_distribution<double> radiusDistribution(
        Parameters::minParticleRadius,
        Parameters::maxParticleRadius
    );

    std::uniform_real_distribution<double> xDistribution(
        0.0,
        Parameters::domainLengthX
    );

    std::uniform_real_distribution<double> yDistribution(
        0.0,
        Parameters::domainLengthY
    );


    const double domainArea =
        Parameters::domainLengthX *
        Parameters::domainLengthY;

    const double targetSolidArea =
        domainArea *
        (1.0 - Parameters::targetPorosity);


    double currentSolidArea = 0.0;
    int attempts = 0;


    while (
        currentSolidArea < targetSolidArea &&
        attempts < Parameters::maxPlacementAttempts
    )
    {
        attempts++;

        double radius = radiusDistribution(generator);

        double x = xDistribution(generator);
        double y = yDistribution(generator);


        // -------------------------------------------------
        // Prevent particles from crossing domain boundaries
        // -------------------------------------------------

        if (
            x - radius < 0.0 ||
            x + radius > Parameters::domainLengthX ||
            y - radius < 0.0 ||
            y + radius > Parameters::domainLengthY
        )
        {
            continue;
        }


        // -------------------------------------------------
        // Check particle overlap
        // -------------------------------------------------

        bool overlap = false;

        for (const auto& particle : particles)
        {
            double dx = x - particle.x;
            double dy = y - particle.y;

            double distance =
                std::sqrt(dx * dx + dy * dy);

            if (distance < radius + particle.radius)
            {
                overlap = true;
                break;
            }
        }


        if (overlap)
        {
            continue;
        }


        // -------------------------------------------------
        // Add particle
        // -------------------------------------------------

        Particle2D newParticle;

        newParticle.x = x;
        newParticle.y = y;
        newParticle.radius = radius;

        particles.push_back(newParticle);


        currentSolidArea +=
            M_PI * radius * radius;
    }


    std::cout
        << "Number of particles: "
        << particles.size()
        << std::endl;

    std::cout
        << "Placement attempts: "
        << attempts
        << std::endl;


    return particles;
}


// =====================================================
// Calculate porosity
// =====================================================

double calculatePorosity(
    const std::vector<Particle2D>& particles
)
{
    double solidArea = 0.0;

    for (const auto& particle : particles)
    {
        solidArea +=
            M_PI *
            particle.radius *
            particle.radius;
    }


    double totalArea =
        Parameters::domainLengthX *
        Parameters::domainLengthY;


    double porosity =
        1.0 - solidArea / totalArea;


    return porosity;
}


// =====================================================
// Calculate average radius
// =====================================================

double calculateAverageRadius(
    const std::vector<Particle2D>& particles
)
{
    if (particles.empty())
    {
        return 0.0;
    }


    double radiusSum = 0.0;

    for (const auto& particle : particles)
    {
        radiusSum += particle.radius;
    }


    return radiusSum /
           static_cast<double>(particles.size());
}