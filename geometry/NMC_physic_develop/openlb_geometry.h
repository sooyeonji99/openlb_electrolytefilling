#ifndef OPENLB_GEOMETRY_H
#define OPENLB_GEOMETRY_H

#include "olb2D.h"
// #include "olb2D.hh"

#include "parameters.h"

#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <cmath>


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

        // Extending the lattice below y=0 introduces round-off at electrode
        // grid lines. Snap only within a tiny tolerance (1e-7 grid spacing).
        const T tolerance = static_cast<T>(param::dx * 1.e-7);
        if (
            x < -tolerance ||
            x > static_cast<T>(param::lx) + tolerance ||
            y < -tolerance ||
            y > static_cast<T>(param::ly) + tolerance
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
                std::floor(x / static_cast<T>(param::dx) + T(1.e-7))
            );

        int j =
            static_cast<int>(
                std::floor(y / static_cast<T>(param::dx) + T(1.e-7))
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
//   material 1 = pore / fluid (outside circles)
//   material 2 = solid (inside circles)
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
    // 1. Entire OpenLB domain (electrode + lower reservoir) -> fluid
    //
    // 0 -> 1 (pore / fluid)
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
    // 1 -> 2 (solid)
    // ========================================================

    geometry.rename(
        param::openlbPoreMaterial,
        param::openlbSolidMaterial,
        solidIndicator
    );

    // One additional lattice layer surrounds the unchanged electrode/reservoir.
    // Side-wall corners take priority over supply/outlet pressure boundaries.
    auto& load = geometry.getLoadBalancer();
    for (int iC = 0; iC < load.size(); ++iC) {
        auto& block = geometry.getBlockGeometry(iC);
        block.forSpatialLocations([&](LatticeR<2> r) {
            const auto p = block.getPhysR(r);
            if (block.get(r) == param::openlbPoreMaterial) {
                if (p[0] < -param::dx * 0.5 || p[0] > param::lx + param::dx * 0.5)
                    block.set(r, param::sideWallMaterial);
                else if (p[1] < -param::reservoirHeight - param::dx * 0.5)
                    block.set(r, param::liquidSupplyMaterial);
                else if (p[1] > param::ly + param::gasBufferHeight + param::dx * 0.5)
                    block.set(r, param::gasOutletMaterial);
            }
        });
    }
    geometry.getStatistics().getStatisticsStatus() = true;
    geometry.communicate();


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

// ============================================================
// OpenLB SuperGeometry material-number output functor
//
// ParaView output:
//   0 = outside
//   1 = pore / fluid (outside circles)
//   2 = solid (inside circles)
// ============================================================

template<typename T>
class SuperGeometryMaterialOutput2D
    : public olb::SuperF2D<T,T>
{
private:

    olb::SuperGeometry<T,2>& _geometry;

public:

    SuperGeometryMaterialOutput2D(
        olb::SuperGeometry<T,2>& geometry)
        :
        olb::SuperF2D<T,T>(geometry, 1),
        _geometry(geometry)
    {
        this->getName() = "material";
    }


    bool operator()(
        T output[],
        const int input[]) override
    {
        // input[0] = global cuboid ID
        // input[1] = lattice x
        // input[2] = lattice y

        auto& load =
            _geometry.getLoadBalancer();

        // The VTK writer may ask for a cuboid
        // that is not local to this process.
        if (!load.isLocal(input[0])) {
            return false;
        }

        // Convert global cuboid ID -> local cuboid ID
        const int iCloc =
            load.loc(input[0]);

        const int iX = input[1];
        const int iY = input[2];

        const int material =
            _geometry
                .getBlockGeometry(iCloc)
                .get(iX, iY);

        output[0] =
            static_cast<T>(material);

        return true;
    }
};

// ============================================================
// Write actual OpenLB SuperGeometry material field
// ============================================================

template<typename T>
void writeOpenLBGeometryVTK(
    olb::SuperGeometry<T,2>& geometry)
{
    olb::OstreamManager clout(
        std::cout,
        "OpenLBGeometryVTK"
    );

    clout
        << "Writing OpenLB material geometry..."
        << std::endl;


    SuperGeometryMaterialOutput2D<T>
        materialField(geometry);


    olb::SuperVTMwriter2D<T>
        vtkWriter("openlb_geometry");


    vtkWriter.addFunctor(
        materialField,
        "material"
    );


    vtkWriter.createMasterFile();
    vtkWriter.write(0);


    clout
        << "OpenLB material geometry written."
        << std::endl;
}

}

#endif
