#pragma once

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <vector>

#include <Eigen/Dense>
#include <flags/flags.hpp>

namespace t_opt
{

template<typename T>
using limits = std::numeric_limits<T>;

template<typename T>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>;

// FIXME use RowMajor here ??? we need do some tests to check speed changes
using DVector = Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor>;

using DMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

using String = std::string;

enum class ProblemProperty : uint8_t
{
    Gradient = 1 << 0,
    LipschitzConstant = 1 << 1,
};
using ProblemPropertyFlags = flags::flags<ProblemProperty>;

String
to_string(ProblemProperty p);

class Problem;

class Point
{
public:
    Point();

    Point(const Point& other) = default;

    explicit
    Point(const Problem& problem);

    void
    resize(const Problem& problem);

    void
    swap(Point& other);

    DVector x;
    DVector g;

    double f;
    double g_nrm_1;
    double g_nrm2; // FIXME rename to g_nrm_2
    double g_nrm2_2;
    double g_nrm_inf;

private:
    void
    reset();
};

}

ALLOW_FLAGS_FOR_ENUM(t_opt::ProblemProperty)
