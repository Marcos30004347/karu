#include "algebra/compute/lib/Commom.hpp"

namespace karu::algebra::compute {
	
unsigned int roundToPowerOfTwo(
  unsigned int v
)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v;
}

}
