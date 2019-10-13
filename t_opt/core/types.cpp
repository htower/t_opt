#include "types.hpp"

#include "blas.hpp"
#include "problem.hpp"

namespace t_opt
{

String
to_string(ProblemProperty p)
{
    switch (p)
    {
        case ProblemProperty::Gradient:
        {
            static const String s("gradient");
            return s;
        }

        case ProblemProperty::LipschitzConstant:
        {
            static const String s("lipschitz constant");
            return s;
        }
    };

    // -Wreturn-type warning fix
    return {};
}

Point::Point()
{
    reset();
}

Point::Point(const Problem& problem)
{
    resize(problem);
}

void
Point::reset()
{
    f = g_nrm2 = g_nrm2_2 = limits<double>::quiet_NaN();
}

void
Point::resize(const Problem& problem)
{
    x.resize(problem.size());

    if (problem.has(ProblemProperty::Gradient))
    {
        g.resize(problem.size());
    }
    else
    {
        g.resize(0);
    }

    blas::set_zero(x);
    blas::set_zero(g);

    reset();
}

void
Point::swap(Point& other)
{
    x.swap(other.x);
    g.swap(other.g);

    std::swap(f, other.f);
    std::swap(g_nrm2, other.g_nrm2);
    std::swap(g_nrm_1, other.g_nrm_1);
    std::swap(g_nrm2_2, other.g_nrm2_2);
    std::swap(g_nrm_inf, other.g_nrm_inf);
}

}
