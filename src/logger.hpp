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
    static void print_stdout(std::string const & msg, bool newline = false);

private:
    static std::ofstream file_;
    static std::mutex mutex_;
};
