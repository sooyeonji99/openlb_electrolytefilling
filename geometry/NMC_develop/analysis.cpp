#include "analysis.h"
#include "parameters.h"

#include <queue>
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>


namespace
{

inline int index(int i, int j)
{
    return j * param::nx + i;
}


bool isInside(int i, int j)
{
    return
        i >= 0 &&
        i < param::nx &&
        j >= 0 &&
        j < param::ny;
}

}


// ============================================================
// Main analysis
// ============================================================

PoreAnalysisResult analyzePoreStructure(
    const std::vector<int>& material)
{
    PoreAnalysisResult result{};


    // --------------------------------------------------------
    // Porosity
    // --------------------------------------------------------

    long long poreCount = 0;

    for (int value : material)
    {
        if (value == param::poreMaterial)
            ++poreCount;
    }


    const long long totalCells =
        static_cast<long long>(
            param::nx) *
        param::ny;


    result.porosity =
        static_cast<double>(poreCount) /
        static_cast<double>(totalCells);


    // --------------------------------------------------------
    // Connectivity analysis
    // --------------------------------------------------------

    std::vector<bool> visited(
        totalCells,
        false);

    std::queue<int> queue;


    // Start from bottom boundary
    for (int i = 0; i < param::nx; ++i)
    {
        const int id = index(i, 0);

        if (
            material[id] ==
            param::poreMaterial
        )
        {
            visited[id] = true;
            queue.push(id);
        }
    }


    const int di4[4] =
        {1, -1, 0, 0};

    const int dj4[4] =
        {0, 0, 1, -1};


    while (!queue.empty())
    {
        const int current =
            queue.front();

        queue.pop();


        const int i =
            current % param::nx;

        const int j =
            current / param::nx;


        for (int n = 0; n < 4; ++n)
        {
            const int ni =
                i + di4[n];

            const int nj =
                j + dj4[n];


            if (!isInside(ni, nj))
                continue;


            const int nid =
                index(ni, nj);


            if (
                !visited[nid] &&
                material[nid] ==
                param::poreMaterial
            )
            {
                visited[nid] = true;
                queue.push(nid);
            }
        }
    }


    long long connectedPoreCount = 0;

    result.verticalConnectivity = false;


    for (int j = 0; j < param::ny; ++j)
    {
        for (int i = 0; i < param::nx; ++i)
        {
            const int id =
                index(i, j);


            if (
                material[id] ==
                    param::poreMaterial &&
                visited[id]
            )
            {
                ++connectedPoreCount;


                if (j == param::ny - 1)
                {
                    result.verticalConnectivity =
                        true;
                }
            }
        }
    }


    if (poreCount > 0)
    {
        result.connectedPoreFraction =
            static_cast<double>(
                connectedPoreCount) /
            static_cast<double>(
                poreCount);
    }
    else
    {
        result.connectedPoreFraction = 0.0;
    }


    // --------------------------------------------------------
    // Local pore radius
    //
    // Distance from each pore cell to nearest solid cell.
    //
    // Simple direct search version.
    // --------------------------------------------------------

    std::vector<std::pair<int,int>>
        solidCells;


    for (int j = 0; j < param::ny; ++j)
    {
        for (int i = 0; i < param::nx; ++i)
        {
            if (
                material[index(i,j)] ==
                param::solidMaterial
            )
            {
                solidCells.emplace_back(i,j);
            }
        }
    }


    result.poreRadiusValues.clear();

    double sumRadius = 0.0;
    double maximumRadius = 0.0;


    for (int j = 0; j < param::ny; ++j)
    {
        for (int i = 0; i < param::nx; ++i)
        {
            if (
                material[index(i,j)] !=
                param::poreMaterial
            )
                continue;


            double minimumDistanceSquared =
                1.0e100;


            for (const auto& solid :
                 solidCells)
            {
                const double di =
                    static_cast<double>(
                        i - solid.first);

                const double dj =
                    static_cast<double>(
                        j - solid.second);


                const double distanceSquared =
                    di * di +
                    dj * dj;


                if (
                    distanceSquared <
                    minimumDistanceSquared
                )
                {
                    minimumDistanceSquared =
                        distanceSquared;
                }
            }


            const double radius =
                std::sqrt(
                    minimumDistanceSquared
                ) *
                param::dx;


            result.poreRadiusValues
                .push_back(radius);


            sumRadius += radius;

            maximumRadius =
                std::max(
                    maximumRadius,
                    radius);
        }
    }


    if (
        !result.poreRadiusValues.empty()
    )
    {
        result.meanPoreRadius =
            sumRadius /
            static_cast<double>(
                result
                .poreRadiusValues
                .size());
    }
    else
    {
        result.meanPoreRadius = 0.0;
    }


    result.maxPoreRadius =
        maximumRadius;


    return result;
}


// ============================================================
// Write pore-radius distribution
// ============================================================

void writePoreRadiusDistributionCSV(
    const std::vector<double>& values)
{
    std::ofstream file(
        param::poreRadiusCsvFileName);


    file << "pore_radius_m,"
         << "pore_radius_um\n";


    for (double radius : values)
    {
        file
            << radius
            << ","
            << radius * 1.0e6
            << "\n";
    }


    file.close();


    std::cout
        << "Pore radius CSV written: "
        << param::poreRadiusCsvFileName
        << std::endl;
}