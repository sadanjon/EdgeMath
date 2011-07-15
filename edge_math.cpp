#include <cmath>
#include "edge_math.h"

namespace edge {
namespace math {

static unsigned int factorial(int t) {
    unsigned int result = t;
    for (int f = t-1; f > 0; --f) {
        result *= f;
    }

    return result;
}


} // namespace math
} // namespace edge


