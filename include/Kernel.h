#ifndef FASTIBEM_KERNEL_H
#define FASTIBEM_KERNEL_H

#include <complex>
#include <algorithm>
#include <complex>

#include "Point3D.h"
#include "base.h"

namespace fastibem {
    
    template<typename T>
    class Kernel {
        
    public:
        /// Typedef of kernel return type
        typedef T ReturnType;
        
        /// Evaluation function which must be overloaded
        virtual T evalSLP(const nurbs::Point3D& xs,
                       const nurbs::Point3D& xf) const = 0;
        
        virtual T evalDLP(const nurbs::Point3D& xs,
                          const nurbs::Point3D& xf,
                          const nurbs::Point3D& n) const = 0;
        
        void setMaterialProp(const std::string& s,
                             const double val)
        {
            mMaterialProps[s] = val;
        }
        
        /// material property getter
        std::pair<bool, double> materialPropVal(const std::string& s) const
        {
            auto search = mMaterialProps.find(s);
            if(search != mMaterialProps.end())
                return std::make_pair(true, search->second);
            else {
                std::cout << "Couldn't find material parameter: " << s << "\n";
                return std::make_pair(false, 0.0);
            }
        }
        
        
    protected:
        
        /// Construct with one material property
        Kernel(const std::pair<std::string, double>& m)
        {
            mMaterialProps.insert(m);
        }
        
        /// Construct with material properties map
        Kernel(const std::map<std::string, double>& m)
        : mMaterialProps(m) {}
        
        
    private:
        
        /// Map of material properties (if required)
        std::map<std::string, double> mMaterialProps;
        
    };
    
    /// Helmholtz specialisation
    class HelmholtzKernel : public Kernel<std::complex<double>> {
        
    public:
        
        /// Default constructor
        HelmholtzKernel()
        : HelmholtzKernel(std::make_pair("k", 1.0))
        {
            std::cout << "Constructing default Helmholtz kernel with wavenumber of 1.0\n";
        }
        
        HelmholtzKernel(const std::pair<std::string, double>& m)
        : Kernel<std::complex<double>>(m)
        {
            std::string s = m.first;
            std::transform(s.begin(),
                           s.end(),
                           s.begin(),
                           ::tolower);
            if(s != "wavenumber" && s != "k") {
                throw std::runtime_error("Must specify wavenumber during Helmholtz kernel construction with "
                                         "string specifier 'wavenumber' or 'k'. Setting default value of 1.0");
                mWavenumber = 1.0;
            }
            else
                mWavenumber = m.second;
        }
        
        ReturnType evalSLP(const nurbs::Point3D& xs,
                           const nurbs::Point3D& xf) const
        {
            const double r = dist(xs, xf);
            return std::exp(std::complex<double>(0.0, 1.0 * wavenumber() * r)) / (4.0 * nurbs::PI * r);
        }
        
        ReturnType evalDLP(const nurbs::Point3D& xs,
                           const nurbs::Point3D& xf,
                           const nurbs::Point3D& n) const
        {
            const double r = dist(xs, xf);
            const std::complex<double> const1(0.0, wavenumber() * r);
            const double drdn = nurbs::dot((xf-xs) / r, n);
            return std::exp(const1) / (4.0 * nurbs::PI * r * r ) * (const1 - 1.0) * drdn;
        }
        
    private:
        
        double wavenumber() const { return mWavenumber; }
        
        double mWavenumber;
        
    };
    
}


#endif