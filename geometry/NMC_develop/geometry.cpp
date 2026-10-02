#include "geometry.h"
#include "parameters.h"

#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>


namespace
{

bool particlesOverlap(
    const Particle2D& a,
    const Particle2D& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;

    const double distanceSquared =
        dx * dx + dy * dy;

    const double minimumDistance =
        a.radius +
        b.radius +
        param::minimumParticleGap;

    return distanceSquared <
           minimumDistance * minimumDistance;
}


bool insideDomain(const Particle2D& particle)
{
    return
        particle.x - particle.radius >= 0.0 &&
        particle.x + particle.radius <= param::lx &&
        particle.y - particle.radius >= 0.0 &&
        particle.y + particle.radius <= param::ly;
}


double calculateAnalyticalPorosity(
    const std::vector<Particle2D>& particles)
{
    double solidArea = 0.0;

    for (const auto& p : particles)
    {
        solidArea +=
            M_PI * p.radius * p.radius;
    }

    const double totalArea =
        param::lx * param::ly;

    return 1.0 - solidArea / totalArea;
}

}


// ============================================================
// Generate particles
// ============================================================

std::vector<Particle2D> generateParticles()
{
    std::vector<Particle2D> particles;

    std::mt19937 generator(param::randomSeed);

    std::uniform_real_distribution<double>
        radiusDistribution(
            param::minParticleRadius,
            param::maxParticleRadius);

    std::uniform_real_distribution<double>
        xDistribution(0.0, param::lx);

    std::uniform_real_distribution<double>
        yDistribution(0.0, param::ly);


    int attempts = 0;

    double currentPorosity = 1.0;


    while (
        currentPorosity >
        param::targetPorosity +
        param::porosityTolerance
    )
    {
        if (attempts >= param::maxPlacementAttempts)
        {
            std::cout
                << "Maximum placement attempts reached."
                << std::endl;

            break;
        }


        Particle2D candidate;

        candidate.radius =
            radiusDistribution(generator);

        candidate.x =
            xDistribution(generator);

        candidate.y =
            yDistribution(generator);


        if (!insideDomain(candidate))
        {
            ++attempts;
            continue;
        }


        bool overlap = false;

        for (const auto& particle : particles)
        {
            if (particlesOverlap(
                    candidate,
                    particle))
            {
                overlap = true;
                break;
            }
        }


        if (!overlap)
        {
            particles.push_back(candidate);

            currentPorosity =
                calculateAnalyticalPorosity(
                    particles);

            // Successful placement:
            // restart consecutive failure counter
            attempts = 0;
        }
        else
        {
            ++attempts;
        }
    }


    std::cout
        << "Particle generation finished."
        << std::endl;

    std::cout
        << "Analytical porosity = "
        << currentPorosity
        << std::endl;

    std::cout
        << "Number of particles = "
        << particles.size()
        << std::endl;


    return particles;
}


// ============================================================
// Build material grid
// ============================================================

std::vector<int> buildMaterialGrid(
    const std::vector<Particle2D>& particles)
{
    std::vector<int> material(
        param::nx * param::ny,
        param::poreMaterial);


    for (int j = 0; j < param::ny; ++j)
    {
        for (int i = 0; i < param::nx; ++i)
        {
            const double x =
                (i + 0.5) * param::dx;

            const double y =
                (j + 0.5) * param::dx;


            for (const auto& p : particles)
            {
                const double rx = x - p.x;
                const double ry = y - p.y;

                if (
                    rx * rx + ry * ry
                    <=
                    p.radius * p.radius
                )
                {
                    material[
                        j * param::nx + i
                    ] =
                        param::solidMaterial;

                    break;
                }
            }
        }
    }


    return material;
}


// ============================================================
// Particle statistics
// ============================================================

double calculateAverageParticleRadius(
    const std::vector<Particle2D>& particles)
{
    if (particles.empty())
        return 0.0;


    double sum = 0.0;

    for (const auto& p : particles)
        sum += p.radius;


    return sum /
           static_cast<double>(
               particles.size());
}


double calculateMinimumParticleRadius(
    const std::vector<Particle2D>& particles)
{
    if (particles.empty())
        return 0.0;


    double value =
        particles.front().radius;

    for (const auto& p : particles)
        value =
            std::min(value, p.radius);


    return value;
}


double calculateMaximumParticleRadius(
    const std::vector<Particle2D>& particles)
{
    if (particles.empty())
        return 0.0;


    double value =
        particles.front().radius;

    for (const auto& p : particles)
        value =
            std::max(value, p.radius);


    return value;
}