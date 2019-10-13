#pragma once

#include "types.hpp"

namespace t_opt
{

class Point;
class Problem;

struct LineSearchProbe
{
    double step;
    double f;
};

class LineSearchMethod
{
public:

    virtual void
    setup(Problem& problem) = 0;

    virtual LineSearchProbe
    search(Problem& problem, const Point& point, const DVector& dir, bool dir_is_gradient, double start_step) = 0;

protected:
    inline double
    fix_step(double step, bool dir_is_gradient)
    {
        return dir_is_gradient ? -std::abs(step) : std::abs(step);
    }
};

}
