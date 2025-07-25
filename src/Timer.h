#ifndef TIMER_H_
#define TIMER_H_

#include <chrono>
#include <ctime>
#include <sstream>

namespace timer
{
    class Timer {
        private:
            std::chrono::time_point<std::chrono::high_resolution_clock> start;
            std::chrono::time_point<std::chrono::high_resolution_clock> end;
        public:
            Timer() { start = std::chrono::high_resolution_clock::now(); }
            void stop() { end = std::chrono::high_resolution_clock::now(); }
            int64_t time() const { return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); }
            static std::string current_time() {
                std::time_t currentTime = std::time(nullptr);
                std::tm* localTime = std::localtime(&currentTime);
                std::ostringstream oss;
                oss << std::put_time(localTime, "%Y-%m-%d %H:%M:%S");
                return oss.str();
            }
    };
}

#endif
