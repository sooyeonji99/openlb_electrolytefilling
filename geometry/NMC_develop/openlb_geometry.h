#ifndef OPENLB_GEOMETRY_H
#define OPENLB_GEOMETRY_H

#include "olb2D.h"
#include "olb2D.hh"

#include "parameters.h"

#include <vector>
#include <iostream>
#include <algorithm>

namespace olb
{

// ============================================================
// Indicator generated directly from the master material grid
//
// Custom grid:
//   0 = pore
//   1 = solid
//
// Indicator:
//   true  = solid
//   false = pore
//
// IMPORTANT:
// The custom grid is defined by CELL CENTERS:
//
//   x_i = (i + 0.5) * dx
//   y_j = (j + 0.5) * dx
//
// OpenLB queries physical coordinates. We map each coordinate
// to the corresponding custom-grid cell.
// ============================================================

template<typename T>
class MaterialGridIndicator2D
    : public IndicatorF2D<T>
{
private:

    const std::vector<int>& _material;

public:

    explicit MaterialGridIndicator2D(
        const std::vector<int>& material)
        : IndicatorF2D<T>(),
          _material(material)
    {
        // Physical extent of the master geometry
        this->_myMin =
            Vector<T,2>(
                T(0),
                T(0)
            );

        this->_myMax =
            Vector<T,2>(
                static_cast<T>(param::lx),
                static_cast<T>(param::ly)
            );
    }


    bool operator()(
        bool output[1],
        const T input[2]) override
    {
        // ----------------------------------------------------
        // Physical coordinate
        // ----------------------------------------------------

        const T x =
            input[0];

        const T y =
            input[1];


        // ----------------------------------------------------
        // Reject coordinates outside the domain
        // ----------------------------------------------------

        if (
            x < T(0) ||
            x > static_cast<T>(param::lx) ||
            y < T(0) ||
            y > static_cast<T>(param::ly)
        )
        {
            output[0] = false;

            return true;
        }


        // ----------------------------------------------------
        // Physical coordinate -> master-grid index
        //
        // OpenLB contains nodes at:
        //
        // x = 0, dx, 2dx, ..., lx
        //
        // whereas our master geometry contains 1000 cells.
        //
        // floor(x/dx) maps the OpenLB coordinate to the
        // corresponding master-grid cell.
        // ----------------------------------------------------

        int i =
            static_cast<int>(
                x / static_cast<T>(param::dx)
            );

        int j =
            static_cast<int>(
                y / static_cast<T>(param::dx)
            );


        // ----------------------------------------------------
        // Special treatment of the upper/right boundary
        //
        // x = lx would give i = nx.
        // Clamp it to the final master-grid cell.
        // ----------------------------------------------------

        i = std::clamp(
            i,
            0,
            param::nx - 1
        );

        j = std::clamp(
            j,
            0,
            param::ny - 1
        );


        const std::size_t index =
            static_cast<std::size_t>(j)
            *
            static_cast<std::size_t>(param::nx)
            +
            static_cast<std::size_t>(i);


        // ----------------------------------------------------
        // Master-grid material:
        //
        // 0 = pore  -> false
        // 1 = solid -> true
        // ----------------------------------------------------

        output[0] =
            (
                _material[index]
                ==
                param::solidMaterial
            );


        return true;
    }
};


// ============================================================
// Prepare OpenLB geometry directly from master material grid
//
// OpenLB:
//   material 1 = solid
//   material 2 = pore / fluid
// ============================================================

template<typename T>
void prepareOpenLBGeometry(
    SuperGeometry<T,2>& geometry,
    const std::vector<int>& material)
{
    OstreamManager clout(
        std::cout,
        "OpenLBGeometry"
    );


    clout
        << "Preparing OpenLB geometry from master grid..."
        << std::endl;


    // ========================================================
    // 1. Entire OpenLB domain -> fluid
    //
    // 0 -> 2
    // ========================================================

    geometry.rename(
        0,
        param::openlbPoreMaterial
    );


    // ========================================================
    // 2. Create solid indicator directly from master grid
    // ========================================================

    MaterialGridIndicator2D<T>
        solidIndicator(material);


    // ========================================================
    // 3. Fluid -> solid wherever master grid says solid
    //
    // 2 -> 1
    // ========================================================

    geometry.rename(
        param::openlbPoreMaterial,
        param::openlbSolidMaterial,
        solidIndicator
    );


    // ========================================================
    // IMPORTANT:
    //
    // Do NOT call:
    //
    // geometry.clean();
    // geometry.innerClean();
    //
    // because the master porous structure should be preserved.
    // ========================================================


    // ========================================================
    // 4. Geometry check
    // ========================================================

    geometry.checkForErrors();


    clout
        << "OpenLB geometry mapping completed."
        << std::endl;


    clout
        << "Material "
        << param::openlbSolidMaterial
        << " = solid"
        << std::endl;


    clout
        << "Material "
        << param::openlbPoreMaterial
        << " = pore/fluid"
        << std::endl;


    geometry.print();
}

}

#endif