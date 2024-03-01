#ifndef CORE_UTILS_H
#define CORE_UTILS_H
#include <iterator>
#include <algorithm>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <random>

#include "../core_types.h"

/** Vector printing to streams */
template <typename T>
std::ostream& operator<< (std::ostream& os, const std::vector<T>& v) {
    if ( !v.empty() ) {
        //os << '[';
        std::copy (v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
        //os << "\b\b]\n";
    }
    return os;
}

template<typename T>
class Logspace {
    private:
        T curValue, base;

    public:
        Logspace(T first, T base) : curValue(first), base(base) {}

        T operator()() {
            T retval = curValue;
            curValue *= base;
            return retval;
        }
};

v_type logspace(double start, double stop, int num, double base);
template<typename T>
v_type linspace(T start_in, T end_in, int num_in)
{

  v_type ls(num_in);

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return ls; }
  if (num == 1)
    {
      ls[0] = start;
      return ls;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      //ls.push_back(start + delta * i);
      ls[i] = start + delta * i;
    }
  ls[num_in-1] = end;
  return ls;
}




template <typename T>
std::string to_string(const T& value,int n) {
    std::stringstream ss;
    ss <<  std::setprecision(n) << value;
    return ss.str();
}

template <typename T>
std::string to_string(const T& value) {
    std::stringstream ss;
    ss << value;
    return ss.str();
}

std::string to_upper(const std::string s);

template <typename F>
double central_difference(F f, double val, double h) {
    return (f(val+h)-f(val-h))/(2.0*h);
}
size_t number_of_digits (size_t i);
size_t max_number(size_t i);    // max_number(6) -> 999999

template<typename T>
T random(T range_from, T range_to) {
    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    std::uniform_real_distribution<T>    distr(range_from, range_to);
    return distr(generator);
}

/** Find element in std::vector and return its index.
 *
 * If not found returns -1.
 */
template <typename T>
int get_index(std::vector<T> v, T K) {
    auto it = find(v.begin(), v.end(), K);
    int index=-1;
    if (it != v.end()) {
        index = it - v.begin();
    }
    return index;
}

//std::vector<double> calc_rho(Structure &st);
#endif
