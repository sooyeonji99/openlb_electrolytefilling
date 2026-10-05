#include "vtk_writer.h"
#include "parameters.h"

#include <fstream>
#include <iostream>


void writeGeometryVTK(
    const std::vector<int>& material)
{
    std::ofstream file(
        param::vtkFileName
    );


    if (!file.is_open())
    {
        std::cerr
            << "Error: Could not open VTK file."
            << std::endl;

        return;
    }


    file
        << "# vtk DataFile Version 3.0\n";

    file
        << "Cathode porous geometry\n";

    file
        << "ASCII\n";

    file
        << "DATASET STRUCTURED_POINTS\n";


    file
        << "DIMENSIONS "
        << param::nx
        << " "
        << param::ny
        << " "
        << 1
        << "\n";


    file
        << "ORIGIN 0 0 0\n";


    file
        << "SPACING "
        << param::dx
        << " "
        << param::dx
        << " "
        << param::dx
        << "\n";


    file
        << "POINT_DATA "
        << param::nx * param::ny
        << "\n";


    file
        << "SCALARS material int 1\n";

    file
        << "LOOKUP_TABLE default\n";


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
            file
                << material[
                    j * param::nx + i
                ]
                << "\n";
        }
    }


    file.close();


    std::cout
        << "VTK file written: "
        << param::vtkFileName
        << std::endl;
}