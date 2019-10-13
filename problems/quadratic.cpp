#include "quadratic.hpp"

#include "core/problem.hpp"

namespace problem
{

Quadratic::Quadratic(size_t n)
    : Problem("Test quadratic problem", n, ProblemProperty::Gradient)
    , m_k(n)
    , m_w(n)
{
    for (size_t i = 0; i < n; ++i)
    {
        m_k[i] = i + 1;
    }
}

void
Quadratic::f(Point& p)
{
    double k = 1.0;

    p.f = 0.0;
    for (size_t i = 0; i < m_size; ++i)
    {
        p.f += k * (p.x[i] - 1.0) * (p.x[i] - 1.0);
        k += 1.0;
    }
}

void
Quadratic::df(Point& p)
{
    double k = 2.0;
    for (size_t i = 0; i < m_size; ++i)
    {
        p.g[i] = k * (p.x[i] - 1.0);
        k += 2.0;
    }
}

}
