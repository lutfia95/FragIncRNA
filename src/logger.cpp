#include "logger.hpp"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>

std::ofstream Logger::file_;
std::mutex Logger::mutex_;

namespace
{

std::string current_timestamp()
{
    using namespace std::chrono;
    auto now = system_clock::now();
    std::time_t t = system_clock::to_time_t(now);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    return oss.str();
}

void log_impl(std::string const & level, std::string const & msg,
              std::ofstream & file, std::mutex & mtx)
{
    std::lock_guard<std::mutex> lock{mtx};
    std::string line = "[" + level + "] " + current_timestamp() + " " + msg;
    // no stdout here: keep stdout clean for progress bar only
    if (file.is_open())
        file << line << '\n';
}


} // namespace

void Logger::init(std::string const & path)
{
    std::lock_guard<std::mutex> lock{mutex_};
    if (file_.is_open())
        file_.close();
    file_.open(path, std::ios::out | std::ios::trunc);
}

void Logger::info(std::string const & msg)
{
    log_impl("INFO", msg, file_, mutex_);
}

void Logger::warn(std::string const & msg)
{
    log_impl("WARN", msg, file_, mutex_);
}

void Logger::error(std::string const & msg)
{
    log_impl("ERROR", msg, file_, mutex_);
}

void Logger::print_stdout(std::string const & msg, bool newline)
{
    std::lock_guard<std::mutex> lock{mutex_};
    std::cout << msg;
    if (newline)
        std::cout << '\n';
    std::cout << std::flush;
}
