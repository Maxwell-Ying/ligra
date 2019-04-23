#ifndef MY_UTIL_H
#define MY_UTIL_H

#include <chrono>
#include <ctime>
#include <random>

#define DEBUG 1

#ifdef DEBUG
#define V_V_V_V_V_V_logstart         if(1){
#define V_V_V_V_V_V_logif(cond)      if(cond){
#else
#define V_V_V_V_V_V_logstart         {
#define V_V_V_V_V_V_logif(cond)      {
#endif
#define V_V_V_V_V_V_logend           }
#define V_V_V_V_V_V_logprepare       ;
#define V_V_V_V_V_V_logpreparefinish ;


/*     chrono defines       */
class ChronoTimer {
public:
  ChronoTimer() : beg_(clock_::now()) {}
  void start() { beg_ = clock_::now(); }
  double elapsed() const {
    return std::chrono::duration_cast<nano_>(clock_::now() - beg_).count();
  }
  double elapsed_milli() const {
    return std::chrono::duration_cast<milli_>(clock_::now() - beg_).count();
  }
  double moment() const {
    auto tp = clock_::now();
    return std::chrono::duration_cast<nano_>(tp.time_since_epoch()).count();
  }
  double moment_milli() const {
    auto tp = clock_::now();
    return std::chrono::duration_cast<milli_>(tp.time_since_epoch()).count();
  }

private:
  using clock_ = std::chrono::high_resolution_clock;
  using nano_ = std::chrono::nanoseconds;
  using milli_ = std::chrono::milliseconds;
  std::chrono::time_point<clock_> beg_;
};

uintT get_new_from(uintT current_max, uintT bound = 5) {
  if (bound == 0) {
    return current_max;
  }
  uintT current_min;
  if (current_max >= bound) {
    current_min = current_max - bound;
  } else {
    current_min = 0;
  }
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(current_min, current_max);

  return dis(gen);
}

#endif