#pragma once

#include "core/method.hpp"

namespace t_opt
{

namespace local
{

class AGMsDR : public Method
{
public:
    explicit
    AGMsDR(LineSearchMethod& ls, double epsilon);

protected:
    void
    before(Problem& problem, Point& point) override;

    bool
    iteration(Problem& problem, Point& point, size_t iter) override;

private:
    LineSearchMethod& ls;
    double ls_start_step = 1.0;
    double ls_step = 1.0;

    double epsilon;

    double A;
    Point y;
    DVector v;
    DVector vx;
};

}

}
