
#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <vector>

// ============================================================
// Pore structure analysis results
//
// Circle interior = solid electrode particle
// Circle exterior = pore / fluid region
// ============================================================

struct PoreAnalysisResult
{
    // Fraction of pore cells in total domain
    double porosity = 0.0;

    // Whether pore space connects bottom to top
    bool verticalConnectivity = false;

    // Fraction of pore cells accessible from bottom
    double connectedPoreFraction = 0.0;
};

// ============================================================
// Main pore analysis
// ============================================================

PoreAnalysisResult analyzePoreStructure(
    const std::vector<int>& material
);

#endif
