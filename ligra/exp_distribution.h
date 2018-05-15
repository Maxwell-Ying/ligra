#include <random>

using namespace std;
struct expDist {
  default_random_engine e;
  
  
  
//   float getRaw() {
//     return expd(e);
//   }
  int getRand(int maxrand) {
    int ret;
    exponential_distribution<float> expd(1.0);
    do {
      ret = (int) lround(10 * expd(e));
    }
    while (ret >= maxrand);
    return maxrand - ret;
  }
};