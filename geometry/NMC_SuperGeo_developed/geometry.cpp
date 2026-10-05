#include "geometry.h"
#include "parameters.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <vector>


// ============================================================
// Generate random non-overlapping particles
// ============================================================

std::vector<Particle2D> generateParticles()
{
    std::vector<Particle2D> particles;

    std::mt19937 generator(
        param::randomSeed
    );

    std::uniform_real_distribution<double> radiusDistribution(
        param::minParticleRadius,
        param::maxParticleRadius
    );

    std::uniform_real_distribution<double> xDistribution(
        0.0,
        param::lx
    );

    std::uniform_real_distribution<double> yDistribution(
        0.0,
        param::ly
    );


    // --------------------------------------------------------
    // Total domain area
    // --------------------------------------------------------

    const double domainArea =
        param::lx * param::ly;


    // --------------------------------------------------------
    // Target solid area
    //
    // porosity = pore area / total area
    //
    // Therefore:
    //
    // solid fraction = 1 - porosity
    // --------------------------------------------------------

    const double targetSolidArea =
        domainArea
        *
        (1.0 - param::targetPorosity);


    double currentSolidArea =
        0.0;


    std::size_t placementAttempts =
        0;


    // ========================================================
    // Particle placement loop
    // ========================================================

    while (
        currentSolidArea < targetSolidArea
        &&
        placementAttempts < param::maxPlacementAttempts
    )
    {
        ++placementAttempts;


        // ----------------------------------------------------
        // Generate particle radius
        // ----------------------------------------------------

        const double radius =
            radiusDistribution(generator);


        // ----------------------------------------------------
        // Particle center
        //
        // Keep the entire particle inside the domain.
        // ----------------------------------------------------

        if (
            2.0 * radius > param::lx
            ||
            2.0 * radius > param::ly
        )
        {
            continue;
        }


        std::uniform_real_distribution<double>
            particleXDistribution(
                radius,
                param::lx - radius
            );

        std::uniform_real_distribution<double>
            particleYDistribution(
                radius,
                param::ly - radius
            );


        const double x =
            particleXDistribution(generator);

        const double y =
            particleYDistribution(generator);


        // ----------------------------------------------------
        // Check overlap with existing particles
        // ----------------------------------------------------

        bool overlap =
            false;


        for (const auto& existing : particles)
        {
            const double dx =
                x - existing.x;

            const double dy =
                y - existing.y;


            const double distanceSquared =
                dx * dx
                +
                dy * dy;


            const double minimumDistance =
                radius
                +
                existing.radius;


            if (
                distanceSquared
                <
                minimumDistance * minimumDistance
            )
            {
                overlap =
                    true;

                break;
            }
        }


        if (overlap)
        {
            continue;
        }


        // ----------------------------------------------------
        // Accept particle
        // ----------------------------------------------------

        Particle2D particle;

        particle.x =
            x;

        particle.y =
            y;

        particle.radius =
            radius;


        particles.push_back(
            particle
        );


        currentSolidArea +=
            M_PI
            *
            radius
            *
            radius;
    }


    std::cout
        << "Particle generation completed.\n";


    std::cout
        << "Placement attempts = "
        << placementAttempts
        << "\n";


    std::cout
        << "Analytical solid fraction = "
        << currentSolidArea / domainArea
        << "\n";


    std::cout
        << "Analytical porosity = "
        << 1.0 - currentSolidArea / domainArea
        << "\n";


    if (
        currentSolidArea
        <
        targetSolidArea
    )
    {
        std::cout
            << "WARNING: Target porosity was not reached "
            << "before maximum placement attempts.\n";
    }


    return particles;
}


// ============================================================
// Calculate average particle radius
// ============================================================

double calculateAverageParticleRadius(
    const std::vector<Particle2D>& particles)
{
    if (particles.empty())
    {
        return 0.0;
    }


    double sum =
        0.0;


    for (const auto& particle : particles)
    {
        sum +=
            particle.radius;
    }


    return sum
        /
        static_cast<double>(
            particles.size()
        );
}


// ============================================================
// Calculate minimum particle radius
// ============================================================

double calculateMinimumParticleRadius(
    const std::vector<Particle2D>& particles)
{
    if (particles.empty())
    {
        return 0.0;
    }


    double minimumRadius =
        std::numeric_limits<double>::max();


    for (const auto& particle : particles)
    {
        minimumRadius =
            std::min(
                minimumRadius,
                particle.radius
            );
    }


    return minimumRadius;
}


// ============================================================
// Calculate maximum particle radius
// ============================================================

double calculateMaximumParticleRadius(
    const std::vector<Particle2D>& particles)
{
    if (particles.empty())
    {
        return 0.0;
    }


    double maximumRadius =
        0.0;


    for (const auto& particle : particles)
    {
        maximumRadius =
            std::max(
                maximumRadius,
                particle.radius
            );
    }


    return maximumRadius;
}


// ============================================================
// Build master material grid
//
// Material convention:
//
//     0 = pore
//     1 = solid particle
//
// Grid:
//
//     nx × ny
//
// Index:
//
//     index = j * nx + i
//
// Cell-center coordinate:
//
//     x = (i + 0.5) * dx
//     y = (j + 0.5) * dx
//
// This material grid is the MASTER GEOMETRY.
// The same grid will later be mapped to OpenLB.
// ============================================================

std::vector<int> buildMaterialGrid(
    const std::vector<Particle2D>& particles)
{
    std::cout
        << "\nBuilding material grid...\n";


    // --------------------------------------------------------
    // Initially every cell is pore
    // --------------------------------------------------------

    std::vector<int> material(
        static_cast<std::size_t>(
            param::nx
        )
        *
        static_cast<std::size_t>(
            param::ny
        ),
        param::poreMaterial
    );


    // ========================================================
    // Rasterize particles
    // ========================================================

    for (const auto& particle : particles)
    {
        // ----------------------------------------------------
        // Determine particle bounding box in grid coordinates
        //
        // This avoids checking every particle against every
        // grid cell.
        // ----------------------------------------------------

        int iMin =
            static_cast<int>(
                std::floor(
                    (particle.x - particle.radius)
                    /
                    param::dx
                )
            );


        int iMax =
            static_cast<int>(
                std::floor(
                    (particle.x + particle.radius)
                    /
                    param::dx
                )
            );


        int jMin =
            static_cast<int>(
                std::floor(
                    (particle.y - particle.radius)
                    /
                    param::dx
                )
            );


        int jMax =
            static_cast<int>(
                std::floor(
                    (particle.y + particle.radius)
                    /
                    param::dx
                )
            );


        // ----------------------------------------------------
        // Clamp bounding box to domain
        // ----------------------------------------------------

        iMin =
            std::max(
                0,
                iMin
            );

        jMin =
            std::max(
                0,
                jMin
            );


        iMax =
            std::min(
                param::nx - 1,
                iMax
            );

        jMax =
            std::min(
                param::ny - 1,
                jMax
            );


        const double radiusSquared =
            particle.radius
            *
            particle.radius;


        // ----------------------------------------------------
        // Check cells within particle bounding box
        // ----------------------------------------------------

        for (
            int j = jMin;
            j <= jMax;
            ++j
        )
        {
            const double y =
                (
                    static_cast<double>(j)
                    +
                    0.5
                )
                *
                param::dx;


            for (
                int i = iMin;
                i <= iMax;
                ++i
            )
            {
                const double x =
                    (
                        static_cast<double>(i)
                        +
                        0.5
                    )
                    *
                    param::dx;


                const double deltaX =
                    x
                    -
                    particle.x;


                const double deltaY =
                    y
                    -
                    particle.y;


                const double distanceSquared =
                    deltaX * deltaX
                    +
                    deltaY * deltaY;


                // --------------------------------------------
                // Cell center lies inside particle
                // --------------------------------------------

                if (
                    distanceSquared
                    <=
                    radiusSquared
                )
                {
                    const std::size_t index =
                        static_cast<std::size_t>(j)
                        *
                        static_cast<std::size_t>(
                            param::nx
                        )
                        +
                        static_cast<std::size_t>(i);


                    material[index] =
                        param::solidMaterial;
                }
            }
        }
    }


    std::cout
        << "Material grid completed.\n";


    return material;
}