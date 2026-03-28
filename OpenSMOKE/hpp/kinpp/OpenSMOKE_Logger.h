#pragma once
 
#include <fstream>
#include <string>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <stdexcept>
 
// Log severity levels
enum class LogLevel {
    INFO,
    WARNING,
    ERROR
};
 
/**
 * Logger ? writes timestamped, levelled messages to a log file.
 *
 * Usage:
 *   Logger log("simulation.log");
 *   log.info("Simulation started");
 *   log.warning("Step 42: convergence slow");
 *   log.error("NaN detected in field U");
 *
 *   // Or use the generic write() with an explicit level:
 *   log.write(LogLevel::INFO, "Checkpoint reached at t = 3.14");
 */
class Logger {
public:
 
    // Open (or create) the log file. Throws std::runtime_error on failure.
    // If append == true, new messages are added to an existing file;
    // otherwise the file is truncated at startup.
    explicit Logger(const std::string& filename, bool append = false)
        : m_filename(filename)
    {
        auto mode = std::ios::out | (append ? std::ios::app : std::ios::trunc);
        m_file.open(filename, mode);
        if (!m_file.is_open())
            throw std::runtime_error("Logger: cannot open file \"" + filename + "\"");
 
        write(LogLevel::INFO, "---- Log opened ----");
    }
 
    // Flush and close the log file on destruction.
    ~Logger() {
        if (m_file.is_open()) {
            write(LogLevel::INFO, "---- Log closed ----");
            m_file.close();
        }
    }
 
    // Non-copyable, non-movable (owns an open file stream).
    Logger(const Logger&)            = delete;
    Logger& operator=(const Logger&) = delete;
    Logger(Logger&&)                 = delete;
    Logger& operator=(Logger&&)      = delete;
 
    // ------------------------------------------------------------------ //
    //  Public interface                                                    //
    // ------------------------------------------------------------------ //
 
    void info   (const std::string& msg) { write(LogLevel::INFO,    msg); }
    void warning(const std::string& msg) { write(LogLevel::WARNING, msg); }
    void error  (const std::string& msg) { write(LogLevel::ERROR,   msg); }
 
    // Generic entry point ? use when the level is determined at runtime.
    void write(LogLevel level, const std::string& msg) {
        if (!m_file.is_open()) return;
        m_file << "[" << timestamp() << "] "
               << "[" << levelTag(level) << "] "
               << msg << "\n";
        m_file.flush();          // ensure the message survives a crash
    }
 
    // Convenience: write a section separator with an optional title.
    void section(const std::string& title = "") {
        std::string line(60, '-');
        if (!title.empty()) {
            // Centre the title inside the dashes
            std::string padded = "  " + title + "  ";
            size_t start = (line.size() > padded.size())
                           ? (line.size() - padded.size()) / 2 : 0;
            line.replace(start, padded.size(), padded);
        }
        m_file << line << "\n";
        m_file.flush();
    }
 
    const std::string& filename() const { return m_filename; }
 
private:
 
    // Returns the current wall-clock time as "YYYY-MM-DD HH:MM:SS"
    static std::string timestamp() {
        auto now  = std::chrono::system_clock::now();
        auto time = std::chrono::system_clock::to_time_t(now);
        std::tm tm{};
#if defined(_WIN32)
        localtime_s(&tm, &time);
#else
        localtime_r(&time, &tm);
#endif
        std::ostringstream ss;
        ss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
        return ss.str();
    }
 
    static const char* levelTag(LogLevel level) {
        switch (level) {
            case LogLevel::INFO:    return "INFO   ";
            case LogLevel::WARNING: return "WARNING";
            case LogLevel::ERROR:   return "ERROR  ";
        }
        return "UNKNOWN";
    }
 
    std::string   m_filename;
    std::ofstream m_file;
};
