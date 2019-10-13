#include "tntp.hpp"

#include <fstream>
#include <map>
#include <set>

#include <fmt/color.h>

namespace tntp
{

using File = std::fstream;
using NodesSet = std::set<uint32_t>;
using EdgesMap = std::map<uint64_t, uint32_t>;

inline uint64_t
edge_key(uint32_t source, uint32_t target)
{
    return ((uint64_t)source << 32) + target;
}

// FIXME move into core ?
class StringView
{
public:
    static const size_t
    npos = String::npos;

    StringView(const String& str)
        : StringView(str.data(), str.size())
    {
    }

    StringView(const char* data, size_t size)
        : m_data(data)
        , m_size(size)
        , m_trimmed(false)
    {
    }

    inline const char*
    data() const
    {
        return m_data;
    }

    inline size_t
    size() const
    {
        return m_size;
    }

    bool
    empty() const
    {
        return m_size == 0;
    }

    bool
    starts_with(char c)
    {
        return !empty() && (m_data[0] == c);
    }

    inline StringView&
    trim()
    {
        if (m_trimmed == false)
        {
            while (m_size > 0 && std::isspace(m_data[0]))
            {
                ++m_data;
                --m_size;
            }

            while (m_size > 0 && std::isspace(m_data[m_size - 1]))
            {
                --m_size;
            }

            m_trimmed = true;
        }

        return *this;
    }

    size_t
    find_first_of(char c, size_t start) const
    {
        for (size_t i = start; i < m_size; ++i)
        {
            if (m_data[i] == c)
            {
                return i;
            }
        }

        return npos;
    }

    void
    split(char delimeter, std::vector<StringView>& result, int max_count = -1) const
    {
        result.clear();
        for (size_t start = 0, end = 0; end < m_size; start = end + 1)
        {
            end = find_first_of(delimeter, start);
            if (end == String::npos)
            {
                end = m_size;
            }

            result.emplace_back(&m_data[start], end - start);

            if ((max_count > 0) && result.size() == (size_t)max_count)
            {
                return;
            }
        }
    }

    inline String
    to_string() const
    {
        return String(m_data, m_size);
    }

private:
  const char* m_data;
  size_t m_size;
  bool m_trimmed;
};

class ParserError : public std::exception
{
public:
    ParserError(const String& text)
    {
        static auto bright_red = fmt::fg(fmt::terminal_color::bright_red);

        m_message = fmt::format("{}: {}", fmt::format(bright_red, "parser error"), text);
    }

    ParserError(const String& text, const String& file_name, size_t line_num)
        : ParserError(text)
    {
        m_message += fmt::format("\n  --> {}:{}", file_name, line_num);
    }

    const char*
    what() const noexcept override
    {
        return m_message.c_str();
    }

private:
    String m_message;
};

long
to_long(StringView& str, const char* value_name, const String file_name, size_t line_num)
{
    if (str.size() > 0)
    {
        str.trim();

        char* end = nullptr;
        long value = strtol(str.data(), &end, 10);
        if (end - str.data() == (long)str.size())
        {
            return value;
        }
    }

    auto message = fmt::format("cant parse '{}' as {}", str.to_string(), value_name);
    throw ParserError(message, file_name, line_num);
}

double
to_double(StringView& str, const char* value_name, const String file_name, size_t line_num)
{
    if (str.size() > 0)
    {
        str.trim();

        char* end = nullptr;
        double value = strtod(str.data(), &end);
        if (end - str.data() == (long)str.size())
        {
            return value;
        }
    }

    auto message = fmt::format("cant parse '{}' as {}", str.to_string(), value_name);
    throw ParserError(message, file_name, line_num);
}

void
load_net_data(
    File& file,
    const String& file_name,
    NodesSet& nodes_set,
    EdgesMap& edges_map,
    Data& data)
{
    Edge edge;
    size_t line_num = 0;
    String line_str;
    std::vector<StringView> items;

    data.max_node_index = 0;
    while (std::getline(file, line_str))
    {
        ++line_num;

        auto line = StringView(line_str).trim();

        // skip empty, comment and metadata lines
        if (line.empty() || line.starts_with('~') || line.starts_with('<'))
        {
            continue;
        }

        line.split('\t', items);
        if (items.size() < 7)
        {
            throw ParserError("line contains less than 7 columns", file_name, line_num);
        }

        edge.source = to_long(items[0], "edge source", file_name, line_num);
        edge.target = to_long(items[1], "edge target", file_name, line_num);
        edge.capacity = to_double(items[2], "edge capacity", file_name, line_num);
        edge.free_flow_time = to_double(items[4], "edge free flow", file_name, line_num);
        edge.b = to_double(items[5], "edge B parameter", file_name, line_num);
        edge.p = to_double(items[6], "edge P parameter", file_name, line_num);

        auto key = edge_key(edge.source, edge.target);
        auto it = edges_map.find(key);
        if (it == edges_map.end())
        {
            edges_map.insert(std::make_pair(key, (uint32_t)data.edges.size()));
            data.edges.push_back(edge);
        }
        else
        {
            printf("warning: edge '%u -> %u' already exists, replace with new one\n", edge.source, edge.target);
            data.edges[it->second] = edge;
        }

        nodes_set.insert(edge.source);
        nodes_set.insert(edge.target);

        data.max_node_index = std::max(data.max_node_index, edge.source);
        data.max_node_index = std::max(data.max_node_index, edge.target);
    }

    printf("  nodes: %zu", nodes_set.size());
    if (nodes_set.size() != data.max_node_index)
    {
        printf(" (max node index: %u)", data.max_node_index);
    }
    printf("\n");

    printf("  links: %zu\n", data.edges.size());

}

void
load_trips_data(
    File& file,
    const String& file_name,
    const NodesSet& nodes_set,
    Data& data)
{
    Trip trip;
    size_t line_num = 0;
    size_t o_source = 0;
    bool has_source = false;
    String line_str;
    NodesSet sources_set;
    std::vector<StringView> items;
    std::vector<StringView> items_2;

    auto throw_error = [&file_name, &line_num](const String& message)
    {
        throw ParserError(message, file_name, line_num);
    };

    data.total_flow = 0.0;
    while (std::getline(file, line_str))
    {
        ++line_num;

        auto line = StringView(line_str).trim();

        // skip empty, comment and metadata lines
        if (line.empty() || line.starts_with('~') || line.starts_with('<'))
        {
            continue;
        }

        if (line.starts_with('O') || line.starts_with('o'))
        {
            line.split(' ', items, 2);
            if (items.size() != 2)
            {
                throw_error(fmt::format("wrong data - '{}', cant parse it as trip source (origin)", line.to_string()));
            }

            o_source = to_long(items[1], "trip source (origin)", file_name, line_num);
            if (nodes_set.find(o_source) == nodes_set.end())
            {
                throw_error(fmt::format("trip source '{}' does not present in the net", o_source));
            }

            sources_set.insert(o_source);
            has_source = true;
            continue;
        }

        if (!has_source)
        {
            throw_error("no opening trip source");
        }

        line.split(';', items);
        for (size_t i = 0; i < items.size(); ++i)
        {
            if (items[i].trim().empty())
            {
                continue;
            }

            items[i].split(':', items_2);
            if (items_2.size() < 2)
            {
                throw_error(fmt::format("wrong data - '{}', can't parse it as 'trip_target : trip_flow'",
                                        items[i].to_string()));
            }

            trip.source = o_source;
            trip.target = to_long(items_2[0], "trip target", file_name, line_num);
            if (nodes_set.find(trip.target) == nodes_set.end())
            {
                throw_error(fmt::format("trip target '{}' does not present in the net", trip.target));
            }

            trip.flow = to_double(items_2[1], "trip flow", file_name, line_num);

            trip.flow *= 0.001; // FIXME remove - test with Yura
//             trip.flow *= 0.1; // FIXME remove - test with Yura

            // FIXME add checks:
            // 1) trips doubles - merge them
            // 2) drop if flow <= 0
            // 3) drop if source == target

            data.trips.push_back(trip);
            data.total_flow += trip.flow;
        }
    }

    data.sources = std::vector<uint32_t>(sources_set.begin(), sources_set.end());

    printf("  zones: %zu\n", data.sources.size());
    printf("  trips: %zu\n", data.trips.size());
}

void
load_flow_data(
    File& file,
    bool is_csv,
    const String& file_name,
    const NodesSet& nodes_set,
    const EdgesMap& edges_map,
    Data& data)
{
    bool has_header = false;
    size_t line_num = 0;
    String line_str;
    std::vector<StringView> items;

    auto throw_error = [&file_name, &line_num](const String& message)
    {
        throw ParserError(message, file_name, line_num);
    };

    data.flow.resize(data.edges.size());
    data.flow.setZero();

    while (std::getline(file, line_str))
    {
        ++line_num;

        auto line = StringView(line_str).trim();

        // skip empty and comment lines
        if (line.empty() || line.starts_with('~'))
        {
            continue;
        }

        if (!has_header)
        {
            if (std::isalpha(line.data()[0]))
            {
                has_header = true;
                continue;
            }

            throw_error("it looks like there is no header in flow file");
        }

        if (is_csv)
        {
            line.split(';', items);
            if (items.size() < 4) // FIXME change order in saver software (flow optimize)
            {
                throw_error("line contains less than 3 columns");
            }
        }
        else
        {
            line.split('\t', items);
            if (items.size() < 3)
            {
                throw_error("line contains less than 3 columns");
            }
        }


        const uint32_t source = to_long(items[0], "flow source", file_name, line_num);
        if (nodes_set.find(source) == nodes_set.end())
        {
            throw_error(fmt::format("flow source '{}' does not present in the net", source));
        }

        const uint32_t target = to_long(items[1], "flow target", file_name, line_num);
        if (nodes_set.find(target) == nodes_set.end())
        {
            throw_error(fmt::format("flow target '{}' does not present in the net", source));
        }

        const auto key = edge_key(source, target);
        const auto it = edges_map.find(key);
        if (it == edges_map.end())
        {
            throw_error(fmt::format("edge '{} -> {}' does not present in the net", source, target));
        }

        if (is_csv)
        {
            data.flow[it->second] = to_double(items[3], "flow value", file_name, line_num);
        }
        else
        {
            data.flow[it->second] = to_double(items[2], "flow value", file_name, line_num);
        }
    }

    if (is_csv)
    {
        printf("   flow: %ld (CSV)\n", data.flow.size());
    }
    else
    {
        printf("   flow: %ld\n", data.flow.size());
    }
}

void
load_tntp_data(const String& path, const String& name, Data& data)
{
    String prefix = fmt::format("{0}/{1}/{1}", path, name);
    String net_name = prefix + "_net.tntp";
    String trips_name = prefix + "_trips.tntp";
    String flow_name = prefix + "_flow.tntp";
    String flow_csv_name = prefix + "_flow.csv";

    printf("Loading TNTP data from %s/%s\n", path.c_str(), name.c_str());

    File net_file(net_name, std::ios::in);
    if (net_file.is_open() == false)
    {
        throw ParserError(fmt::format("unable to open file '{}'", net_name));
    }

    File trips_file(trips_name, std::ios::in);
    if (net_file.is_open() == false)
    {
        throw ParserError(fmt::format("unable to open file '{}'", trips_name));
    }

    NodesSet nodes_set;
    EdgesMap edges_map;

    load_net_data(net_file, net_name, nodes_set, edges_map, data);
    load_trips_data(trips_file, trips_name, nodes_set, data);

    File flow_file(flow_csv_name, std::ios::in);
    if (flow_file.is_open())
    {
        load_flow_data(flow_file, true, flow_name, nodes_set, edges_map, data);
    }
    else
    {
        File flow_file(flow_name, std::ios::in);
        if (flow_file.is_open())
        {
            load_flow_data(flow_file, false, flow_name, nodes_set, edges_map, data);
        }
    }
}

}
