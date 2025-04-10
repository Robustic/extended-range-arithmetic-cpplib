#ifndef TIMER_H_
#define TIMER_H_

#include <chrono>

namespace floatingExp2Integer
{
    class Timer {
        private:
            std::chrono::time_point<std::chrono::high_resolution_clock> start;
            std::chrono::time_point<std::chrono::high_resolution_clock> end;
        public:
            Timer() { start = std::chrono::high_resolution_clock::now(); }
            void stop() { end = std::chrono::high_resolution_clock::now(); }
            int64_t time() const { return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); }
    };
}

#endif
