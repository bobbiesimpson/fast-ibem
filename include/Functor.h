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
        std::vector<ReturnType> operator()(const nurbs::Point3D& p) const
        {
            const double x = p[0];
            const double y = p[1];
            const double z = p[2];
            //return {x[0], 0.0, 0.0};
            return
            {
                x*z / std::sqrt(x*x + y*y),
                y*z / std::sqrt(x*x + y*y),
                -std::sqrt(x * x + y * y)
            };
            
//            return { std::cos(x[0]), std::cos(x[1]), std::cos(x[2])};
        }
        
    private:
        
    };
    
    /// Specialisation for an electromagnetic plane wave
    class EmagPlaneWave : public ComplexDFunctor {
        
    public:
        
        /// Constructor
        EmagPlaneWave(const nurbs::Point3D& kvec,
                      const std::vector<std::complex<double>>& pvec)
        :
        mWaveVector(kvec),
        mPolarisationVector(pvec) {}
        
        
        /// 3 components
        unsigned dimN() const
        { return 3; }
        
        /// override function operator
        std::vector<ReturnType> operator()(const nurbs::Point3D& x) const
        {
            const std::complex<double> iconst(0.0, 1.0);
            
            const auto& p = polarVec();
            const auto& k = waveVec();
            
            std::vector<std::complex<double>> result(3);
            const std::complex<double> wave = std::exp(-iconst * dot(k, x));
            for(size_t i = 0; i < 3; ++i)
                result[i] =  p[i] * wave;
            return result;
        }
        
    private:
        
        const nurbs::Point3D& waveVec() const { return mWaveVector; }
        
        const std::vector<std::complex<double>>& polarVec() const { return mPolarisationVector; }
        
        const nurbs::Point3D mWaveVector;
        
        const std::vector<std::complex<double>> mPolarisationVector;
    };
}
#endif
