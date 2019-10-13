#pragma once

#include "core/line_search.hpp"

namespace t_opt
{

namespace line_search
{

struct Parabolic : public LineSearchMethod
{
    Parabolic(uint8_t probes, bool use_gradient);

    void
    setup(Problem& problem) override;

    LineSearchProbe
    search(Problem& problem, const Point& point, const DVector& dir, bool dir_is_gradient, double start_step) override;

private:
    void
    probe(Problem& problem, const Point& point, const DVector& dir, uint8_t i);

    void
    sort_probes(uint8_t count);

    uint8_t probes;
    uint8_t max_probes;
    bool use_gradient;

    Point probe_point;
    LineSearchProbe p[4]; // FIXME rename
};

}

}
