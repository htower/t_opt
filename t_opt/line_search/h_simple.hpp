#pragma once

#include "core/line_search.hpp"

namespace t_opt
{

namespace line_search
{

struct HSimple : public LineSearchMethod
{
    HSimple();
    HSimple(double step_minus_k, double step_plus_k);
    HSimple(double fixed_step);

    void
    setup(Problem& problem) override;

    LineSearchProbe
    search(Problem& problem, const Point& point, const DVector& dir, bool dir_is_gradient, double start_step) override;

    double step_minus_k;
    double step_plus_k;
    double step_min;

private:
    LineSearchProbe
    probe(Problem& problem, const Point& point, const DVector& dir, double step);

    Point probe_point;
    double fixed_step;
};

}

}
