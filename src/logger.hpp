#pragma once

#include <fstream>
#include <mutex>
#include <string>

class Logger
{
public:
    static void init(std::string const & path);
    static void info(std::string const & msg);
    static void warn(std::string const & msg);
    static void error(std::string const & msg);

private:
    static std::ofstream file_;
    static std::mutex mutex_;
};
