#include "analysis.h"

#include <queue>
#include <vector>
#include <cmath>
#include <limits>
#include <utility>

#include <fstream>
#include <algorithm>
#include <string>

// =====================================================
// Check vertical pore connectivity
// =====================================================

bool checkVerticalPoreConnectivity(
    const std::vector<int>& materialGrid,
    int nx,
    int ny
)
{
    std::vector<bool> visited(nx * ny, false);

    std::queue<int> q;


    // -------------------------------------------------
    // Start from pore cells at the bottom boundary
    // -------------------------------------------------

    for (int i = 0; i < nx; ++i)
    {
        int index = i;

        if (materialGrid[index] == 0)
        {
            visited[index] = true;
            q.push(index);
        }
    }


    // Four-neighbor connectivity
    const int dx[4] = {1, -1, 0, 0};
    const int dy[4] = {0, 0, 1, -1};


    while (!q.empty())
    {
        int index = q.front();
        q.pop();

        int x = index % nx;
        int y = index / nx;


        // Reached the top boundary
        if (y == ny - 1)
        {
            return true;
        }


        for (int direction = 0; direction < 4; ++direction)
        {
            int newX = x + dx[direction];
            int newY = y + dy[direction];


            if (
                newX < 0 ||
                newX >= nx ||
                newY < 0 ||
                newY >= ny
            )
            {
                continue;
            }


            int newIndex =
                newY * nx + newX;


            if (
                materialGrid[newIndex] == 0 &&
                !visited[newIndex]
            )
            {
                visited[newIndex] = true;
                q.push(newIndex);
            }
        }
    }


    return false;
}


// =====================================================
// Connected pore fraction
// =====================================================

double calculateConnectedPoreFraction(
    const std::vector<int>& materialGrid,
    int nx,
    int ny
)
{
    std::vector<bool> visited(nx * ny, false);

    std::queue<int> q;


    // Start from bottom pore cells
    for (int i = 0; i < nx; ++i)
    {
        int index = i;

        if (materialGrid[index] == 0)
        {
            visited[index] = true;
            q.push(index);
        }
    }


    const int dx[4] = {1, -1, 0, 0};
    const int dy[4] = {0, 0, 1, -1};


    int connectedPoreCells = 0;


    while (!q.empty())
    {
        int index = q.front();
        q.pop();

        ++connectedPoreCells;

        int x = index % nx;
        int y = index / nx;


        for (int direction = 0; direction < 4; ++direction)
        {
            int newX = x + dx[direction];
            int newY = y + dy[direction];


            if (
                newX < 0 ||
                newX >= nx ||
                newY < 0 ||
                newY >= ny
            )
            {
                continue;
            }


            int newIndex =
                newY * nx + newX;


            if (
                materialGrid[newIndex] == 0 &&
                !visited[newIndex]
            )
            {
                visited[newIndex] = true;
                q.push(newIndex);
            }
        }
    }


    // Count total pore cells
    int totalPoreCells = 0;

    for (int material : materialGrid)
    {
        if (material == 0)
        {
            ++totalPoreCells;
        }
    }


    if (totalPoreCells == 0)
    {
        return 0.0;
    }


    return static_cast<double>(connectedPoreCells)
         / static_cast<double>(totalPoreCells);


    
}

// =====================================================
// Calculate pore distance map
// =====================================================

std::vector<double> calculatePoreDistanceMap(
    const std::vector<int>& materialGrid,
    int nx,
    int ny,
    double dx
)
{
    std::vector<double> distanceMap(
        nx * ny,
        0.0
    );


    // Store coordinates of solid cells
    std::vector<std::pair<int, int>> solidCells;

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            int index = j * nx + i;

            if (materialGrid[index] == 1)
            {
                solidCells.push_back({i, j});
            }
        }
    }


    // For every pore cell, find nearest solid cell
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            int index = j * nx + i;


            // Solid cell
            if (materialGrid[index] == 1)
            {
                distanceMap[index] = 0.0;
                continue;
            }


            double minDistanceSquared =
                std::numeric_limits<double>::max();


            for (const auto& solid : solidCells)
            {
                double deltaX =
                    static_cast<double>(i - solid.first);

                double deltaY =
                    static_cast<double>(j - solid.second);


                double distanceSquared =
                    deltaX * deltaX +
                    deltaY * deltaY;


                if (distanceSquared < minDistanceSquared)
                {
                    minDistanceSquared =
                        distanceSquared;
                }
            }


            distanceMap[index] =
                std::sqrt(minDistanceSquared) * dx;
        }
    }


    return distanceMap;
}

// =====================================================
// Mean local pore diameter
// =====================================================

double calculateMeanPoreDiameter(
    const std::vector<double>& poreDistanceMap,
    const std::vector<int>& materialGrid
)
{
    double sumDiameter = 0.0;
    int poreCells = 0;


    for (std::size_t i = 0;
         i < materialGrid.size();
         ++i)
    {
        if (materialGrid[i] == 0)
        {
            sumDiameter +=
                2.0 * poreDistanceMap[i];

            ++poreCells;
        }
    }


    if (poreCells == 0)
    {
        return 0.0;
    }


    return sumDiameter /
           static_cast<double>(poreCells);
}


// =====================================================
// Maximum local pore diameter
// =====================================================

double calculateMaxPoreDiameter(
    const std::vector<double>& poreDistanceMap,
    const std::vector<int>& materialGrid
)
{
    double maxRadius = 0.0;


    for (std::size_t i = 0;
         i < materialGrid.size();
         ++i)
    {
        if (
            materialGrid[i] == 0 &&
            poreDistanceMap[i] > maxRadius
        )
        {
            maxRadius =
                poreDistanceMap[i];
        }
    }


    return 2.0 * maxRadius;
}

// =====================================================
// Write local pore diameter distribution to CSV
// =====================================================

void writePoreSizeDistributionCSV(
    const std::vector<double>& poreDistanceMap,
    const std::vector<int>& materialGrid,
    double binWidth,
    const char* filename
)
{
    // -------------------------------------------------
    // Find maximum pore diameter
    // -------------------------------------------------

    double maxDiameter = 0.0;

    for (std::size_t i = 0;
         i < materialGrid.size();
         ++i)
    {
        if (materialGrid[i] == 0)
        {
            double diameter =
                2.0 * poreDistanceMap[i];

            if (diameter > maxDiameter)
            {
                maxDiameter = diameter;
            }
        }
    }


    // -------------------------------------------------
    // Number of histogram bins
    // -------------------------------------------------

    int numberOfBins =
        static_cast<int>(
            std::ceil(maxDiameter / binWidth)
        );

    std::vector<int> histogram(
        numberOfBins,
        0
    );


    // -------------------------------------------------
    // Fill histogram
    // -------------------------------------------------

    int totalPoreCells = 0;

    for (std::size_t i = 0;
         i < materialGrid.size();
         ++i)
    {
        if (materialGrid[i] == 0)
        {
            double diameter =
                2.0 * poreDistanceMap[i];

            int bin =
                static_cast<int>(
                    diameter / binWidth
                );

            if (bin >= numberOfBins)
            {
                bin = numberOfBins - 1;
            }

            histogram[bin]++;

            totalPoreCells++;
        }
    }


    // -------------------------------------------------
    // Write CSV
    // -------------------------------------------------

    std::ofstream file(filename);

    if (!file.is_open())
    {
        return;
    }


    file
        << "bin_min_um,"
        << "bin_max_um,"
        << "bin_center_um,"
        << "count,"
        << "fraction\n";


    for (int bin = 0;
         bin < numberOfBins;
         ++bin)
    {
        double binMin =
            bin * binWidth;

        double binMax =
            (bin + 1) * binWidth;

        double binCenter =
            0.5 * (binMin + binMax);


        double fraction =
            static_cast<double>(
                histogram[bin]
            )
            /
            static_cast<double>(
                totalPoreCells
            );


        file
            << binMin * 1e6 << ","
            << binMax * 1e6 << ","
            << binCenter * 1e6 << ","
            << histogram[bin] << ","
            << fraction
            << "\n";
    }


    file.close();
}