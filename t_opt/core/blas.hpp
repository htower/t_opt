#pragma once

#include "types.hpp"

namespace t_opt
{

namespace blas
{

inline void
scal(double a, DVector& x)
{
    x *= a;
}

inline void
copy(const DVector& src, DVector& dest)
{
    dest = src;
//     std::memcpy(dest.data(), src.data(), src.size() * sizeof(double));
}

inline void
scal_copy(double a, const DVector& src, DVector& dest)
{
    dest = a * src;
}

inline void
set_zero(DVector& x)
{
    x.setZero();
//     std::memset(x.data(), 0, x.size() * sizeof(double));
}

inline void
axpy(double a, const DVector& x, DVector& y)
{
    y += a * x;
//     for (uint32_t i = 0; i < x.size(); ++i)
//     {
//         y[i] += a * x[i];
//     }
}

inline void
axpyz(double a, const DVector& x, const DVector& y, DVector& z)
{
    z = a * x + y;
//     for (uint32_t i = 0; i < x.size(); ++i)
//     {
//         z[i] = a * x[i] + y[i];
//     }
}

inline void
axpbyz(double a, const DVector& x, double b, const DVector& y, DVector& z)
{
    z = a * x + b * y;
//     for (uint32_t i = 0; i < x.size(); ++i)
//     {
//         z[i] = a * x[i] + b * y[i];
//     }
}

inline void
xmy(const DVector& x, DVector& y)
{
    y = x - y;
}

inline void
xmyz(const DVector& x, const DVector& y, DVector& z)
{
    z = x - y;
//     for (uint32_t i = 0; i < x.size(); ++i)
//     {
//         z[i] = x[i] - y[i];
//     }
}

inline double
dot(const DVector& x, const DVector& y)
{
    return x.dot(y);

//     // FIXME add vectorization/loop unrolling
//     double d = 0.0;
//     for (uint32_t i = 0; i < x.size(); ++i)
//     {
//         d += x[i] * y[i];
//     }
//     return d;
}

inline double
nrm2(const DVector& x)
{
    return x.stableNorm();
//     return std::sqrt(dot(x, x));
}

inline void
mv(const DMatrix& m, const DVector& x, DVector& y)
{
    y.noalias() = m * x;
}

}

}
