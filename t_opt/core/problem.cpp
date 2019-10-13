#include "problem.hpp"

namespace t_opt
{

Problem::Problem(const String& name, size_t size, ProblemPropertyFlags properties)
    : m_name(name)
    , m_size(size)
    , m_l(limits<double>::quiet_NaN())
    , m_properties(properties)
{
}

}
