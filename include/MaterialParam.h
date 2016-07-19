#ifndef FASTIBEM_MATERIAL_PARAM_H
#define FASTIBEM_MATERIAL_PARAM_H

#include <map>
#include <utility>
#include <algorithm>
#include <string>
#include "Point3D.h"

namespace fastibem
{
    ///
    /// A singleton class for holding material parameters used by kernels.
    /// Examples include wavenumber, Young's modulus, Poisson's ratio.
    ///
    
    class MaterialParam
    {
    public:
        
        /// Get the one and only instance
        static MaterialParam& Instance()
        {
            static MaterialParam theInstance;
            return theInstance;
        }
        
        /// Add a material parameter with a string identifier
        bool addDoubleParam(const std::string& name,
                            const double val)
        {
            auto result = data().insert(std::make_pair(name, val));
            return result.second;
        }
        
        /// Get the parameter. Throws a runtime error if not found
        double getDoubleParam(const std::string& name) const
        {
            auto search = data().find(name);
            if(search == data().end())
                throw std::runtime_error("Cannot find parameter with name: " + name);
            return search->second;
        }
        
        /// Add point data
        bool addPointParam(const std::string& name,
                           const nurbs::Point3D& p)
        {
            auto result = pointdata().insert(std::make_pair(name, p));
            return result.second;
        }
        
        /// Get the parameter. Throws a runtime error if not found
        const nurbs::Point3D& getPointParam(const std::string& name) const
        {
            auto search = pointdata().find(name);
            if(search == pointdata().end())
                throw std::runtime_error("Cannot find parameter with name: " + name);
            return search->second;
        }
        
    private:
        
        /// Typedef of data type
        typedef std::map<std::string, double> DMapType;
        
        /// Typedef of data type
        typedef std::map<std::string, nurbs::Point3D> PointMapType;
        
        /// Constructor (made private)
        MaterialParam() = default;
        
        /// Destructor (made private)
        ~MaterialParam()  = default;
        
        /// Copy constructor (made private)
        MaterialParam(const MaterialParam& m)  = default;
        
        /// Move constructor (made private)
        MaterialParam(MaterialParam&& m) = default;
        
        /// Assignment copy constructor (made private)
        MaterialParam& operator=(const MaterialParam& m) = default;
        
        /// Move assignment (made private)
        MaterialParam& operator=(MaterialParam&& m) = default;
        
        /// Data accessor
        DMapType& data() { return mData; }
        
        /// Const data accessor
        const DMapType& data() const { return mData; }
        
        /// Point data accessor
        PointMapType& pointdata() { return mPointData; }
        
        /// Const point data accessor
        const PointMapType& pointdata() const { return mPointData; }
        
        /// Map of material parameters
        std::map<std::string, double> mData;
        
        /// Map of 3D point data
        std::map<std::string, nurbs::Point3D> mPointData;
        
        /// Overload output operator
        friend std::ostream& operator<<(std::ostream& ost, const MaterialParam& p);
    };
}

#endif