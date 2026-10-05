#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <vector>


struct PoreAnalysisResult
{
    double porosity;

    bool verticalConnectivity;

    double connectedPoreFraction;

    double meanPoreRadius;

    double maxPoreRadius;

    std::vector<double>
        poreRadiusValues;
};


// ============================================================
// Main pore analysis
// ============================================================

PoreAnalysisResult
analyzePoreStructure(
    const std::vector<int>& material
);


// ============================================================
// CSV output
// ============================================================

void writePoreRadiusDistributionCSV(
    const std::vector<double>& poreRadiusValues
);


#endif