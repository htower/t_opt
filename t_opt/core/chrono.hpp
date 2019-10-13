#pragma once

#include <chrono>

namespace t_opt
{

namespace chrono
{

using Point = std::chrono::time_point<std::chrono::system_clock>;

inline Point
now()
{
    return std::chrono::system_clock::now();
}

inline double
ms(const Point& start, const Point& end)
{
    return std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - start).count();
}

inline double
ms(const Point& start)
{
    return ms(start, now());
}

inline double
s(const Point& start, const Point& end)
{
    return std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
}

inline double
s(const Point& start)
{
    return s(start, now());
}

}

}
