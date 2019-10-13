#pragma once

#include <cstdio>
#include <cstdint>

namespace t_opt
{

class Method;
class Point;
struct State;

enum class LoggerDevice : uint8_t
{
    Stdout,
    CsvFile,
};

enum class LoggerMode : uint8_t
{
    Header,
    Value,
};

class Logger
{
public:
    static Logger
    stdout(const Method& method);

    static Logger
    csv(const Method& method, bool resume);

    inline LoggerDevice
    device() const
    {
        return m_device;
    }

    void
    flush();

    void
    put_new_line();

    void
    put_carriage_return();

    void
    put_delimiter();

    void
    put_value(const char* value);

    void
    put_method_name(LoggerMode mode);

    void
    put_time(LoggerMode mode, double time);

    void
    put_iter(LoggerMode mode, size_t iter);

    void
    put_point(LoggerMode mode, const Point& point);

    void
    put_state(LoggerMode mode, const State& state);

    void
    put_ls_step(LoggerMode mode, double ls_step);


private:
    Logger(const Method& method, LoggerDevice device, std::FILE* writer, const char* delimiter);

    void
    put_g_norms(LoggerMode mode, const Point& point);

    const Method& m_method;
    LoggerDevice m_device;
    std::FILE* m_writer;
    const char* m_delimiter;
};

}
