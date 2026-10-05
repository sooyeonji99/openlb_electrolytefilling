#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <vector>


// Check whether a continuous pore pathway exists
// from the bottom boundary to the top boundary.
bool checkVerticalPoreConnectivity(
    const std::vector<int>& materialGrid,
    int nx,
    int ny
);


// Calculate the fraction of pore cells connected
// to the bottom boundary.
double calculateConnectedPoreFraction(
    const std::vector<int>& materialGrid,
    int nx,
    int ny
);


// Calculate local pore radius using distance to
// the nearest solid cell.
std::vector<double> calculatePoreDistanceMap(
    const std::vector<int>& materialGrid,
    int nx,
    int ny,
    double dx
);


// Calculate mean local pore diameter.
double calculateMeanPoreDiameter(
    const std::vector<double>& poreDistanceMap,
    const std::vector<int>& materialGrid
);


// Calculate maximum local pore diameter.
double calculateMaxPoreDiameter(
    const std::vector<double>& poreDistanceMap,
    const std::vector<int>& materialGrid
);

// Write local pore diameter distribution to CSV
void writePoreSizeDistributionCSV(
    const std::vector<double>& poreDistanceMap,
    const std::vector<int>& materialGrid,
    double binWidth,
    const char* filename
);

#endif