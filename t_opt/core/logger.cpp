#include "logger.hpp"

#include "types.hpp"
#include "method.hpp"

namespace t_opt
{

Logger::Logger(const Method& method, LoggerDevice device, std::FILE* writer, const char* delimiter)
    : m_method(method)
    , m_device(device)
    , m_writer(writer)
    , m_delimiter(delimiter)
{
}

Logger
Logger::stdout(const Method& method)
{
    return Logger(method, LoggerDevice::Stdout, ::stdout, " ");
}

Logger
Logger::csv(const Method& method, bool resume)
{
    String file_name(method.name());
    std::transform(file_name.begin(), file_name.end(), file_name.begin(), ::tolower);
    file_name += ".log";

    auto writer = resume ? fopen(file_name.c_str(), "a") : fopen(file_name.c_str(), "w"); // FIXME check result

    return Logger(method, LoggerDevice::CsvFile, writer, ";");
}

void
Logger::flush()
{
    std::fflush(m_writer);
}

void
Logger::put_new_line()
{
    std::fputs("\n", m_writer);
}

void
Logger::put_carriage_return()
{
    switch (m_device)
    {
        case LoggerDevice::Stdout:
            std::fputs("\r", m_writer);
            return;

        case LoggerDevice::CsvFile:
            // FIXME print warning ?
            return;
    };
}

void
Logger::put_delimiter()
{
    std::fputs(m_delimiter, m_writer);
}

void
Logger::put_value(const char* value)
{
    std::fputs(value, m_writer);
}

void
Logger::put_method_name(LoggerMode mode)
{
    static const char* COLUMN_NAME = "method";
    static const size_t COLUMN_WIDTH = strlen(COLUMN_NAME);

    const int width = std::max(m_method.name().length(), COLUMN_WIDTH);

    // Don't put delimiter since method name is first column

    switch (mode)
    {
        case LoggerMode::Header :
            switch (m_device)
            {
                case LoggerDevice::Stdout :
                    fprintf(m_writer, "%-*s", width, COLUMN_NAME);
                    break;

                case LoggerDevice::CsvFile :
                    fprintf(m_writer, "%s", COLUMN_NAME);
                    break;
            }
            break;

        case LoggerMode::Value :
            switch (m_device)
            {
                case LoggerDevice::Stdout :
                    fprintf(m_writer, "%-*s", width, m_method.name().c_str());
                    break;

                case LoggerDevice::CsvFile :
                    fprintf(m_writer, "%s", m_method.name().c_str());
                    break;
            }
            break;
    };
}

void
Logger::put_time(LoggerMode mode, double time)
{
    static const char* COLUMN_NAME = "time";

    static const int COLUMN_PRECISION = 3;
    static const int COLUMN_WIDTH = COLUMN_PRECISION + 5;

    put_delimiter();

    switch (mode)
    {
        case LoggerMode::Header :
            switch (m_device)
            {
                case LoggerDevice::Stdout :
                    fprintf(m_writer, "%*s", COLUMN_WIDTH, COLUMN_NAME);
                    break;

                case LoggerDevice::CsvFile :
                    fprintf(m_writer, "%s", COLUMN_NAME);
                    break;
            }
            break;

        case LoggerMode::Value :
            switch (m_device)
            {
                case LoggerDevice::Stdout :
                    fprintf(m_writer, "%*.*f", COLUMN_WIDTH, COLUMN_PRECISION, time);
                    break;

                case LoggerDevice::CsvFile :
                    fprintf(m_writer, "%.6e", time);
                    break;
            }
            break;
    };
}

void
Logger::put_iter(LoggerMode mode, size_t iter)
{
    static const char* COLUMN_NAME = "iter";
    static const int COLUMN_WIDTH = 9;

    put_delimiter();

    switch (mode)
    {
        case LoggerMode::Header :
            switch (m_device)
            {
                case LoggerDevice::Stdout :
                    fprintf(m_writer, "%*s", COLUMN_WIDTH, COLUMN_NAME);
                    break;

                case LoggerDevice::CsvFile :
                    fprintf(m_writer, "%s", COLUMN_NAME);
                    break;
            }
            break;

        case LoggerMode::Value :
            switch (m_device)
            {
                case LoggerDevice::Stdout :
                    fprintf(m_writer, "%*zu", COLUMN_WIDTH, iter);
                    break;

                case LoggerDevice::CsvFile :
                    fprintf(m_writer, "%zu", iter);
                    break;
            }
            break;
    };
}

void
Logger::put_point(LoggerMode mode, const Point& point)
{
    static const char* F_NAME = "f";

    static const int F_PRECISION = 16;
    static const int F_WIDTH = F_PRECISION + 7;

    put_delimiter();

    switch (mode)
    {
        case LoggerMode::Header :
            switch (m_device)
            {
                case LoggerDevice::Stdout :
                    fprintf(m_writer, " %-*s", F_WIDTH - 1, F_NAME);
                    if (m_method.use(ProblemProperty::Gradient))
                    {
                        put_g_norms(mode, point);
                    }
                    break;

                case LoggerDevice::CsvFile :
                    fprintf(m_writer, "%s", F_NAME);
                    if (m_method.use(ProblemProperty::Gradient))
                    {
                        put_g_norms(mode, point);
                    }
                    break;
            }
            break;

        case LoggerMode::Value :
            switch (m_device)
            {
                case LoggerDevice::Stdout :
                    fprintf(m_writer, "% -*.*e", F_WIDTH, F_PRECISION, point.f);
                    if (m_method.use(ProblemProperty::Gradient))
                    {
                        put_g_norms(mode, point);
                    }
                    break;

                case LoggerDevice::CsvFile :
//                     fprintf(m_writer, "%.*e", F_GN_PRECISION, point.f);
                    fprintf(m_writer, "%.*e", 6, point.f);
                    if (m_method.use(ProblemProperty::Gradient))
                    {
                        put_g_norms(mode, point);
                    }
                    break;
            }
            break;
    };
}

void
Logger::put_state(LoggerMode mode, const State& state)
{
    static const char* F_COLUMN_NAME = "f_count";
    static const char* G_COLUMN_NAME = "g_count";
    static const int COLUMN_WIDTH = 9;

    put_delimiter();

    switch (mode)
    {
        case LoggerMode::Header :
            switch (m_device)
            {
                case LoggerDevice::Stdout :
                    fprintf(m_writer, "%*s", COLUMN_WIDTH, F_COLUMN_NAME);
                    if (m_method.use(ProblemProperty::Gradient))
                    {
                        put_delimiter();
                        fprintf(m_writer, "%*s", COLUMN_WIDTH, G_COLUMN_NAME);
                    }
                    break;

                case LoggerDevice::CsvFile :
                    fprintf(m_writer, "%s", F_COLUMN_NAME);
                    if (m_method.use(ProblemProperty::Gradient))
                    {
                        put_delimiter();
                        fprintf(m_writer, "%s", G_COLUMN_NAME);
                    }
                    break;
            }
            break;

        case LoggerMode::Value :
            switch (m_device)
            {
                case LoggerDevice::Stdout :
                    fprintf(m_writer, "%*zu", COLUMN_WIDTH, state.f_count);
                    if (m_method.use(ProblemProperty::Gradient))
                    {
                        put_delimiter();
                        fprintf(m_writer, "%*zu", COLUMN_WIDTH, state.g_count);
                    }
                    break;

                case LoggerDevice::CsvFile :
                    fprintf(m_writer, "%zu", state.f_count);
                    if (m_method.use(ProblemProperty::Gradient))
                    {
                        put_delimiter();
                        fprintf(m_writer, "%zu", state.g_count);
                    }
                    break;
            }
            break;
    };
}

void
Logger::put_ls_step(LoggerMode mode, double ls_step)
{
    static const char* NAME = "ls_step";
    static const int PRECISION = 4;
    static const int WIDTH = PRECISION + 7;

    put_delimiter();

    switch (mode)
    {
        case LoggerMode::Header :
            switch (m_device)
            {
                case LoggerDevice::Stdout :
                    fprintf(m_writer, " %-*s", WIDTH - 1, NAME);
                    break;

                case LoggerDevice::CsvFile :
                    fprintf(m_writer, "%s", NAME);
                    break;
            }
            break;

        case LoggerMode::Value :
            switch (m_device)
            {
                case LoggerDevice::Stdout :
                    fprintf(m_writer, "% -*.*e", WIDTH, PRECISION, ls_step);
                    break;

                case LoggerDevice::CsvFile :
                    fprintf(m_writer, "%.*e", 6, ls_step);
                    break;
            }
            break;
    };
}

void
Logger::put_g_norms(LoggerMode mode, const Point& point)
{
    static const char* NRM_1_NAME = "g_nrm_1";
    static const char* NRM_2_NAME = "g_nrm_2";
    static const char* NRM_INF_NAME = "g_nrm_inf";

    static const int NRM_PRECISION = 6;
    static const int NRM_WIDTH = NRM_PRECISION + 7;

    put_delimiter();

    switch (mode)
    {
        case LoggerMode::Header :
            switch (m_device)
            {
                case LoggerDevice::Stdout :
                    fprintf(m_writer, " %-*s", NRM_WIDTH - 1, NRM_1_NAME);
                    put_delimiter();
                    fprintf(m_writer, " %-*s", NRM_WIDTH - 1, NRM_2_NAME);
                    put_delimiter();
                    fprintf(m_writer, " %-*s", NRM_WIDTH - 1, NRM_INF_NAME);
                    break;

                case LoggerDevice::CsvFile :
                    fprintf(m_writer, "%s", NRM_1_NAME);
                    put_delimiter();
                    fprintf(m_writer, "%s", NRM_2_NAME);
                    put_delimiter();
                    fprintf(m_writer, "%s", NRM_INF_NAME);
                    break;
            }
            break;

        case LoggerMode::Value :
            switch (m_device)
            {
                case LoggerDevice::Stdout :
                    fprintf(m_writer, "% -*.*e", NRM_WIDTH, NRM_PRECISION, point.g_nrm_1);
                    put_delimiter();
                    fprintf(m_writer, "% -*.*e", NRM_WIDTH, NRM_PRECISION, point.g_nrm2);
                    put_delimiter();
                    fprintf(m_writer, "% -*.*e", NRM_WIDTH, NRM_PRECISION, point.g_nrm_inf);
                    break;

                case LoggerDevice::CsvFile :
                    fprintf(m_writer, "%.*e", NRM_PRECISION, point.g_nrm_1);
                    put_delimiter();
                    fprintf(m_writer, "%.*e", NRM_PRECISION, point.g_nrm2);
                    put_delimiter();
                    fprintf(m_writer, "%.*e", NRM_PRECISION, point.g_nrm_inf);
                    break;
            }
            break;
    };
}

}
