#ifndef FASTIBEM_FUNCTOR_H
#define FASTIBEM_FUNCTOR_H

#include "Point3D.h"
#include <cmath>

namespace fastibem {
    
    template<typename T>
    class Functor {
        
    public:
        
        /// Return typedef
        typedef T ReturnType;
        
        /// Number of dimensions
        virtual unsigned dimN() const = 0;
        
        /// The function operator
        virtual std::vector<ReturnType> operator()(const nurbs::Point3D& x) const = 0;
    };
    
    /// Some typedefs
    typedef Functor<double> DFunctor;
    typedef Functor<std::complex<double>> ComplexDFunctor;
    
    /// Specialisation for an electromagnetic plane wave
    class SinusoidalFunctor : public DFunctor {
        
    public:
        
        /// Constructor
        SinusoidalFunctor() = default;
        
        /// 3 components
        unsigned dimN() const
        { return 3; }
        
        /// override function operator
        std::vector<ReturnType> operator()(const nurbs::Point3D& x) const
        {
            //return {x[0], 0.0, 0.0};
            return { std::cos(x[0]), std::cos(x[1]), std::cos(x[2])};
        }
        
    private:
        
    };
    
    /// Specialisation for an electromagnetic plane wave
    class EmagPlaneWave : public ComplexDFunctor {
        
    public:
        
        /// Constructor
        EmagPlaneWave(const nurbs::Point3D& kvec,
                      const nurbs::Point3D& pvec)
        :
        mWaveVector(kvec),
        mPolarisationVector(pvec) {}
        
        
        /// 3 components
        unsigned dimN() const
        { return 3; }
        
        /// override function operator
        std::vector<ReturnType> operator()(const nurbs::Point3D& x) const
        {
            const auto& p = polarVec();
            const auto& k = waveVec();
            
            std::vector<std::complex<double>> result{p[0], p[1], p[2]};
            const auto wave = std::exp(-std::complex<double>(0.0, dot(k, x)));
            for(auto& r : result)
                r *= wave;
            return result;
        }
        
    private:
        
        const nurbs::Point3D& waveVec() const { return mWaveVector; }
        
        const nurbs::Point3D& polarVec() const { return mPolarisationVector; }
        
        const nurbs::Point3D mWaveVector;
        
        const nurbs::Point3D mPolarisationVector;
    };
}
#endif
