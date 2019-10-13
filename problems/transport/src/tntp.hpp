#pragma once

#include "core/types.hpp"

namespace tntp
{

using t_opt::String;

struct Edge
{
    uint32_t source;
    uint32_t target;
    double   capacity;
    double   free_flow_time;
    double   b;
    double   p;

    double   travel_time;

    inline double
    calc_travel_time(double flow)
    {
        return free_flow_time * (1.0 + b * std::pow(flow / capacity, p));
    }

    inline double
    update_travel_time(double flow)
    {
        return travel_time = calc_travel_time(flow);
    }
};

struct Trip
{
    uint32_t source;
    uint32_t target;
    double   flow;
};

struct Data
{
    std::vector<Edge> edges;
    std::vector<Trip> trips;
    std::vector<uint32_t> sources;

    t_opt::DVector flow;
    double total_flow;

    uint32_t max_node_index;
};

void
load_tntp_data(const String& path, const String& name, Data& data);

}

