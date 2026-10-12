
#include "analysis.h"
#include "parameters.h"

#include <queue>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{

inline int gridIndex(int i, int j)
{
    return j * param::nx + i;
}

inline bool isInside(int i, int j)
{
    return i >= 0 && i < param::nx
        && j >= 0 && j < param::ny;
}

}

// ============================================================
// Analyze pore space outside solid circles
// ============================================================

PoreAnalysisResult analyzePoreStructure(
    const std::vector<int>& material)
{
    PoreAnalysisResult result{};

    const int totalCells = param::nx * param::ny;

    if (material.size() !=
        static_cast<std::size_t>(totalCells))
    {
        throw std::runtime_error(
            "Invalid material grid size"
        );
    }

    std::cout
        << "\nAnalyzing pore structure...\n";

    // --------------------------------------------------------
    // 1. Porosity
    // --------------------------------------------------------

    long long poreCount = 0;

    for (int value : material)
    {
        if (value == param::poreMaterial)
        {
            ++poreCount;
        }
    }

    result.porosity =
        static_cast<double>(poreCount)
        / static_cast<double>(totalCells);

    // --------------------------------------------------------
    // 2. Bottom-connected pore network (4-neighbor BFS)
    // --------------------------------------------------------

    std::vector<bool> visited(
        static_cast<std::size_t>(totalCells),
        false
    );

    std::queue<int> q;

    // Initialize from bottom pore cells (j = 0)
    for (int i = 0; i < param::nx; ++i)
    {
        const int id = gridIndex(i, 0);

        if (material[id] == param::poreMaterial)
        {
            visited[id] = true;
            q.push(id);
        }
    }

    const int di[4] = {1, -1, 0, 0};
    const int dj[4] = {0, 0, 1, -1};

    while (!q.empty())
    {
        const int current = q.front();
        q.pop();

        const int i = current % param::nx;
        const int j = current / param::nx;

        for (int n = 0; n < 4; ++n)
        {
            const int ni = i + di[n];
            const int nj = j + dj[n];

            if (!isInside(ni, nj))
            {
                continue;
            }

            const int nid = gridIndex(ni, nj);

            if (!visited[nid] &&
                material[nid] == param::poreMaterial)
            {
                visited[nid] = true;
                q.push(nid);
            }
        }
    }

    // --------------------------------------------------------
    // 3. Connectivity statistics
    // --------------------------------------------------------

    long long connectedPoreCount = 0;

    for (int j = 0; j < param::ny; ++j)
    {
        for (int i = 0; i < param::nx; ++i)
        {
            const int id = gridIndex(i, j);

            if (material[id] == param::poreMaterial
                && visited[id])
            {
                ++connectedPoreCount;

                if (j == param::ny - 1)
                {
                    result.verticalConnectivity = true;
                }
            }
        }
    }

    result.connectedPoreFraction =
        poreCount > 0
        ? static_cast<double>(connectedPoreCount)
          / static_cast<double>(poreCount)
        : 0.0;

    std::cout
        << "Pore analysis completed.\n";

    return result;
}
