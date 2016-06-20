#include "CollocationAssembly.h"

namespace fastibem {
    
    

    std::vector<long int> calculateThreadBounds(long int parts, long int mem)
    {
        std::vector<long int>bnd;
        long int delta = mem / parts;
        long int reminder = mem % parts;
        long int N1 = 0, N2 = 0;
        bnd.push_back(N1);
        for (long int i = 0; i < parts; ++i) {
            N2 = N1 + delta;
            if (i == parts - 1)
                N2 += reminder;
            bnd.push_back(N2);
            N1 = N2;
        }
        return bnd;
    }
}