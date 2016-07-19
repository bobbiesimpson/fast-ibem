#ifndef FASTIBEM_GALERKIN_ASSEMBLY_H
#define FASTIBEM_GALERKIN_ASSEMBLY_H

#include <stdexcept>
#include "Functor.h"
#include "Kernel.h"

namespace fastibem {
    
    
    
    ///
    /// A class used to represent assembly implementation for Galerkin BE
    /// computations. It is intended to be used with the HLibPro library
    /// which requests submatrices of the dense linear operator.
    ///
    /// This class is also reponsible for computing terms of the force vector.
    template<typename K>
    class GalerkinAssembly {
        
    protected:
        
        /// Typedef of the kernel type
        typedef K KernelType;
        
        /// Typedef of the kernel return type
        typedef typename K::ReturnType DataType;
        
        
    public:
        
        /// Evaluate the submatrix of the LHS matrix
        virtual std::vector<std::vector<DataType>> evalSubmatrixLHS(const std::vector<uint>& rows,
                                                                    const std::vector<uint>& cols) const
        {
            throw std::runtime_error("evalSubmatrixLHS() is not implemented for this class");
        }
        
        /// Evaluate the submatrix of the RHS matrix
        virtual std::vector<std::vector<DataType>> evalSubmatrixRHS(const std::vector<uint>& rows,
                                                                    const std::vector<uint>& cols) const
        {
            throw std::runtime_error("evalSubmatrixRHS() is not implemented for this class");
        }
        
        /// Evaluat the global force vector
        virtual std::vector<DataType> evalForceVec(const Functor<DataType>& func) const
        {
            throw std::runtime_error("evalForceVec() is not implemented for this class");
        }
        
        /// Kernel accessor
        const KernelType& kernel() const { return mKernel; }
        
    protected:
        
        /// Abstract class therefore constructor is protected.
        GalerkinAssembly(const KernelType& k)
        :
        mKernel(k) {}
        
    private:
        
        /// The kernel instance
        const KernelType& mKernel;
    
    };
    
    /// Specialisation for electromagnetic galerkin assembly
    class GalerkinEmagAssembly : public GalerkinAssembly<EmagKernel> {
      

    public:
        
        /// Constructor
        GalerkinEmagAssembly(const EmagKernel& kernel)
        :
        GalerkinAssembly(kernel) {}
        
        virtual std::vector<std::vector<DataType>> evalSubmatrixLHS(const std::vector<uint>& rows,
                                                                    const std::vector<uint>& cols) const
        {
            // TODO: implement equivalent routine as in collocation assembly whereby
            // we integrate over all connected elements and simply discard any terms we
            // don't need.
            
            // We need to determine when to apply each of the singular quadrature
            // transformations (vertex,edge,identical elements).
            
            return std::vector<std::vector<DataType>>{};
        }
        
        virtual std::vector<DataType> evalForceVec(const Functor<DataType>& func) const
        { return std::vector<DataType>{}; }
        
    };
    
}

#endif