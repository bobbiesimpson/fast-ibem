#ifndef FASTIBEM_KERNEL_INTEGRATOR_H
#define FASTIBEM_KERNEL_INTEGRATOR_H

#include "Forest.h"

namespace fastibem {
    
    /// Template class for integrating terms over a NURBS 'forest' for
    /// the kernel K.
    
    template<typename K>
    class KernelIntegrator {
        
    public:
        
        /// Handy typedef for kernel type
        typedef K KernelType;
        
        /// typedef for kernel return type
        typedef T K::ReturnType;
        
        /// Default constructor
        KernelIntegrator()
        : KernelIntegrator(nullptr) {}
        
        /// Construct with a forest
        KernelIntegrator(const nurbs::Forest* f)
        : mpForest(f) {}
        
        /// Evaluate the 
        T eval(const uint icolloc,
               const uint ielem) const;
        
    private:
        
        /// Pointer to the forest (mesh)
        const nurbs::Forest* mpForest;
    };
    
}
#endif