#include <iostream>

#include "Forest.h"
#include "Geometry.h"

#include "CollocationAssembly.h"

int main(const int argc, const char* argv[])
{
    try {
        std::cout << "Running fast-ibem integrator test...\n";
        std::ifstream ifs(argv[1]);
        if(!ifs)
            nurbs::error("Cannot open file for reading\n");
        nurbs::Geometry g;
        if(!g.loadHBSFile(ifs))
            nurbs::error("Failed to load geometry from hbs data");
        nurbs::Forest forest(g);

        fastibem::HelmholtzKernel hkernel(std::make_pair("k", 10.0));
        fastibem::CollocationAssembly<fastibem::HelmholtzKernel> assembly(&forest,
                                                                          hkernel,
                                                                          true);
        
        return EXIT_SUCCESS;
    }
    catch(const std::exception& e) {
        std::cerr << e.what() << "\n";
        return EXIT_FAILURE;
    }
    catch(...) {
        std::cerr << "Unknown exception occured\n";
        return EXIT_FAILURE;
    }
}