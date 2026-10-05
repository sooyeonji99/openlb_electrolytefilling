#include "analysis.h"
#include "parameters.h"

#include <queue>
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <utility>


namespace
{

inline int gridIndex(
    int i,
    int j)
{
    return
        j * param::nx + i;
}


bool isInside(
    int i,
    int j)
{
    return
        i >= 0
        &&
        i < param::nx
        &&
        j >= 0
        &&
        j < param::ny;
}

}


// ============================================================
// Main analysis
// ============================================================

PoreAnalysisResult
analyzePoreStructure(
    const std::vector<int>& material)
{
    PoreAnalysisResult result{};


    std::cout
        << "\nAnalyzing pore structure..."
        << std::endl;


    // ========================================================
    // 1. Porosity
    // ========================================================

    long long poreCount =
        0;


    for (
        int value :
        material
    )
    {
        if (
            value
            ==
            param::poreMaterial
        )
        {
            ++poreCount;
        }
    }


    const long long totalCells =
        static_cast<long long>(
            param::nx
        )
        *
        static_cast<long long>(
            param::ny
        );


    result.porosity =
        static_cast<double>(
            poreCount
        )
        /
        static_cast<double>(
            totalCells
        );


    // ========================================================
    // 2. Vertical pore connectivity
    // ========================================================

    std::vector<bool> visited(
        totalCells,
        false
    );


    std::queue<int> queue;


    // Start BFS from bottom pore cells
    for (
        int i = 0;
        i < param::nx;
        ++i
    )
    {
        const int id =
            gridIndex(
                i,
                0
            );


        if (
            material[id]
            ==
            param::poreMaterial
        )
        {
            visited[id] =
                true;

            queue.push(
                id
            );
        }
    }


    const int di[4] =
        {1, -1, 0, 0};

    const int dj[4] =
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


        for (
            int n = 0;
            n < 4;
            ++n
        )
        {
            const int ni =
                i + di[n];

            const int nj =
                j + dj[n];


            if (
                !isInside(
                    ni,
                    nj
                )
            )
            {
                continue;
            }


            const int nid =
                gridIndex(
                    ni,
                    nj
                );


            if (
                !visited[nid]
                &&
                material[nid]
                ==
                param::poreMaterial
            )
            {
                visited[nid] =
                    true;

                queue.push(
                    nid
                );
            }
        }
    }


    long long connectedPoreCount =
        0;


    result.verticalConnectivity =
        false;


    for (
        int j = 0;
        j < param::ny;
        ++j
    )
    {
        for (
            int i = 0;
            i < param::nx;
            ++i
        )
        {
            const int id =
                gridIndex(
                    i,
                    j
                );


            if (
                material[id]
                ==
                param::poreMaterial
                &&
                visited[id]
            )
            {
                ++connectedPoreCount;


                if (
                    j
                    ==
                    param::ny - 1
                )
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
                connectedPoreCount
            )
            /
            static_cast<double>(
                poreCount
            );
    }
    else
    {
        result.connectedPoreFraction =
            0.0;
    }


    // ========================================================
    // 3. Collect solid cells
    // ========================================================

    std::vector<
        std::pair<int,int>
    > solidCells;


    solidCells.reserve(
        totalCells
    );


    for (
        int j = 0;
        j < param::ny;
        ++j
    )
    {
        for (
            int i = 0;
            i < param::nx;
            ++i
        )
        {
            if (
                material[
                    gridIndex(i,j)
                ]
                ==
                param::solidMaterial
            )
            {
                solidCells.emplace_back(
                    i,
                    j
                );
            }
        }
    }


    // ========================================================
    // 4. Local pore radius
    //
    // Distance from pore cell to nearest solid cell.
    //
    // NOTE:
    // This is a brute-force implementation.
    // ========================================================

    result.poreRadiusValues.clear();


    double sumRadius =
        0.0;

    double maximumRadius =
        0.0;


    std::cout
        << "Calculating pore-radius field..."
        << std::endl;


    for (
        int j = 0;
        j < param::ny;
        ++j
    )
    {
        for (
            int i = 0;
            i < param::nx;
            ++i
        )
        {
            if (
                material[
                    gridIndex(i,j)
                ]
                !=
                param::poreMaterial
            )
            {
                continue;
            }


            double minimumDistanceSquared =
                1.0e100;


            for (
                const auto& solid :
                solidCells
            )
            {
                const double deltaI =
                    static_cast<double>(
                        i - solid.first
                    );

                const double deltaJ =
                    static_cast<double>(
                        j - solid.second
                    );


                const double distanceSquared =
                    deltaI * deltaI
                    +
                    deltaJ * deltaJ;


                if (
                    distanceSquared
                    <
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
                )
                *
                param::dx;


            result.poreRadiusValues.push_back(
                radius
            );


            sumRadius +=
                radius;


            maximumRadius =
                std::max(
                    maximumRadius,
                    radius
                );
        }
    }


    if (
        !result.poreRadiusValues.empty()
    )
    {
        result.meanPoreRadius =
            sumRadius
            /
            static_cast<double>(
                result
                .poreRadiusValues
                .size()
            );
    }
    else
    {
        result.meanPoreRadius =
            0.0;
    }


    result.maxPoreRadius =
        maximumRadius;


    std::cout
        << "Pore analysis completed."
        << std::endl;


    return result;
}


// ============================================================
// CSV
// ============================================================

void writePoreRadiusDistributionCSV(
    const std::vector<double>& values)
{
    std::ofstream file(
        param::poreRadiusCsvFileName
    );


    if (!file.is_open())
    {
        std::cerr
            << "Error: Cannot open pore-radius CSV."
            << std::endl;

        return;
    }


    file
        << "pore_radius_m,"
        << "pore_radius_um\n";


    for (
        double radius :
        values
    )
    {
        file
            << radius
            << ","
            << radius * 1.0e6
            << "\n";
    }


    file.close();


    std::cout
        << "Pore-radius CSV written: "
        << param::poreRadiusCsvFileName
        << std::endl;
}