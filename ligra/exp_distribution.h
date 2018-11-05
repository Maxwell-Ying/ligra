#include <random>
#include <ctime>
using namespace std;
struct expDist {
  default_random_engine e;
  
  expDist() {
    e.seed(time(0) * time(0));
  }
  
//   float getRaw() {
//     return expd(e);
//   }
  int getRand(int maxrand) {
    int ret;
    exponential_distribution<float> expd(0.5);
    do {
      ret = (int) lround(10 * expd(e));
    }
    while (ret >= maxrand);
    return maxrand - ret;
  }
};
