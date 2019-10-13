#pragma once

#include "core/problem.hpp"

namespace problem
{

using namespace t_opt;

struct Quadratic : public Problem
{
    Quadratic(size_t n);

    void
    f(Point& p) override;

    void
    df(Point& p) override;

private:
    DVector m_k;
    DVector m_w;
};

}
