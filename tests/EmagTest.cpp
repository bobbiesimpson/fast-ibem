#include <iostream>
#include "GalerkinAssembly.h"
#include "Functor.h"
#include "Kernel.h"
#include "Point3D.h"

using namespace nurbs;

int main(int argc, char* argv[])
{
    fastibem::EmagPlaneWave functor(Point3D(0.0, 1.0, 0.0),
                                Point3D(0.0, 1.0, 0.0));
    
    std::cout << functor(Point3D(0.0, 0.0, 0.0)) << "\n"
                << functor(Point3D(0.0, 1.0, 0.0)) <<  "\n"
                << functor(Point3D(0.0, 0.0, 1.0)) << "\n";
    
    const double k = 3.14;
    fastibem::EmagKernel ekernel(std::make_pair("wavenumber", k));
    fastibem::GalerkinEmagAssembly assembly(ekernel);
    std::cout << assembly.evalForceVec(functor) << "\n";
    return EXIT_SUCCESS;
}