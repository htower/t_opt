#pragma once

#include "types.hpp"

namespace t_opt
{

class Problem
{
public:
    inline const String&
    name() const
    {
        return m_name;
    }

    inline size_t
    size() const
    {
        return m_size;
    }

    inline size_t
    dual_size() const
    {
        return m_dual_size;
    }

    inline bool
    has(ProblemProperty p) const
    {
        return bool(m_properties & p);
    }

    inline double
    L() const
    {
        return m_l;
    }

    virtual void
    f(Point& p) = 0;

    virtual void
    df(Point& p) = 0;

    virtual void // FIXME remove
    emoe(Point& point) {};

    virtual void
    dual_x(Point& p, Point& dual_p) {}

    virtual void
    dual_f(Point& dual_p) {}

protected:
    Problem(const Problem& problem) = default;
    Problem(const String& name, size_t size, ProblemPropertyFlags properties = ProblemPropertyFlags());

    String m_name;
    size_t m_size;
    size_t m_dual_size;
    double m_l;
    ProblemPropertyFlags m_properties;
};

}
