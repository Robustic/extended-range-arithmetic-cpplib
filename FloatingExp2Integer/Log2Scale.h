#ifndef LOG2SCALE_H_
#define LOG2SCALE_H_

#include <vector>
#include <cmath>

namespace floatingExp2Integer
{
    class Log2Scale {
    private:
        double log2scale;
    public:
        Log2Scale() { log2scale = 0.0; }
        Log2Scale(double d) { log2scale = std::log2(d); }

        void log2_to(const double from) { log2scale = from; }
        static void log2s_to(const std::vector<double>& from, std::vector<floatingExp2Integer::Log2Scale>& to) {
            for (size_t i = 0; i < to.size(); i++) {
                to[i].log2_to(from[i]);
            }
        }

        static void as_vector(const std::vector<floatingExp2Integer::Log2Scale>& from, std::vector<double>& to) {
            for (size_t i = 0; i < to.size(); i++) {
                to[i] = from[i].log2scale;
            }
        }

        double as_double() { return std::exp2(log2scale); }
        double as_log2() { return log2scale; }

        Log2Scale& operator+=(Log2Scale d) {
            log2scale += d.log2scale;
            return *this;
        }

        Log2Scale& operator*=(Log2Scale d) {
            log2scale *= d.log2scale;
            return *this;
        }

        //void sum(const std::vector<floatingExp2Integer::Log2Scale>& dblValues) {
        //    
        //}

        //void multiply(const std::vector<floatingExp2Integer::Log2Scale>& dblValues) {
        //    
        //}
    };
}

#endif
