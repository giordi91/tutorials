#include <cstdint>

namespace cpp_tools {
namespace algorithms {

inline uint32_t findHighestBit(uint32_t v) {
  uint32_t r = 0;

  while (v >>= 1) {
    r++;
  }
  return r;
}

uint32_t simpleMultSlow(uint32_t a, uint32_t b) {

  uint32_t result = 0;
  for (uint32_t bi = 0; bi <= 31; ++bi) {

    uint32_t bmask = 1 << bi;
    uint32_t bvalue = (b & bmask) >> bi;
    uint32_t currValue = 0;
    for (uint32_t ai = 0; ai <= 31; ++ai) {
      uint32_t amask = 1 << ai;
      uint32_t avalue = (a & amask) >> ai;

      currValue |= bvalue ? avalue << (ai + bi) : 0;
    }
    result += currValue;
  }

  return result;
}

uint32_t simpleMultFaster(uint32_t a, uint32_t b) {

  uint32_t abit = findHighestBit(a);
  uint32_t bbit = findHighestBit(b);
  if (abit + bbit > 31)
    return 0;

  uint32_t result = 0;
  for (uint32_t bi = 0; bi <= bbit; ++bi) {

    uint32_t bmask = 1 << bi;
    uint32_t bvalue = (b & bmask) >> bi;
    uint32_t aPosValue = a << bi;
    result += bvalue ? aPosValue : 0;
  }

  return result;
}

uint32_t karatsubaOneLevel(uint32_t x, uint32_t y, uint32_t size) {
  // generating mask
  int halfSize = size >> 1;

  uint32_t lowerHalfMask = (1 << (halfSize + 1)) - 1;
  uint32_t upperHalfMask = ~lowerHalfMask;

  // extracting upper lower part
  uint32_t a = (x & upperHalfMask) >> halfSize;
  uint32_t b = x & lowerHalfMask;
  uint32_t c = (y & upperHalfMask) >> halfSize;
  uint32_t d = y & lowerHalfMask;

  // computing steps
  uint32_t step1 = (a * c);
  uint32_t step2 = (b * d);
  uint32_t step3 = (a + b) * (c + d);
  ;
  uint32_t gauss = step3 - step2 - step1;
  return (2 << size) * step1 + (2 << halfSize) * gauss + step2;
}

uint32_t karatsuba(uint32_t x, uint32_t y, uint32_t size) {

  if (size <= 4) {
    return simpleMultFaster(x, y);
    // return x * y;
  }
  // generating mask
  int halfSize = size >> 1;
  uint32_t lowerHalfMask = (1 << (halfSize + 1)) - 1;
  uint32_t upperHalfMask = ~lowerHalfMask;

  // extracting upper lower part
  uint32_t a = (x & upperHalfMask) >> halfSize;
  uint32_t b = x & lowerHalfMask;
  uint32_t c = (y & upperHalfMask) >> halfSize;
  uint32_t d = y & lowerHalfMask;

  // computing steps
  uint32_t step1 = karatsuba(a, c, halfSize);
  uint32_t step2 = karatsuba(b, d, halfSize);
  uint32_t step3 = karatsuba((a + b), (c + d), halfSize);
  uint32_t gauss = step3 - step2 - step1;
  return (step1 << (size)) + (gauss << (halfSize)) + step2;
}

template <uint32_t SIZE> uint32_t karatsubaTemplate(uint32_t x, uint32_t y) {
  // TODO(giordi) might want to use constexpr here if we want to jump
  // on c++17
  if (SIZE <= 4) {
    return x * y;
  }
  // generating mask
  const uint32_t halfSize = SIZE >> 1;
  uint32_t lowerHalfMask = (1 << (halfSize + 1)) - 1;
  uint32_t upperHalfMask = ~lowerHalfMask;

  // extracting upper lower part
  uint32_t a = (x & upperHalfMask) >> halfSize;
  uint32_t b = x & lowerHalfMask;
  uint32_t c = (y & upperHalfMask) >> halfSize;
  uint32_t d = y & lowerHalfMask;

  // computing steps
  uint32_t step1 = karatsubaTemplate<halfSize>(a, c);
  uint32_t step2 = karatsubaTemplate<halfSize>(b, d);
  uint32_t step3 = karatsubaTemplate<halfSize>((a + b), (c + d));
  uint32_t gauss = step3 - step2 - step1;
  if (SIZE != 32) {
    // NOTE(ignore the warning here, the compiler is not able to figure out that
    // i am not shifting by anything bigger than 16
    return (step1 << (SIZE)) + (gauss << (halfSize)) + step2;
  } else {
    return (gauss << (halfSize)) + step2;
  }
}

constexpr uint32_t karatsubaConstExpr(uint32_t x, uint32_t y, uint32_t SIZE) {
  if (SIZE <= 4) {
    // return ((x == 1) & (y == 1)) ? 2 : simpleMultFaster( x, y);
    return ((x == 1) & (y == 1)) ? 2 : x * y;
  }
  // generating mask
  uint32_t halfSize = SIZE >> 1;
  uint32_t lowerHalfMask = (1 << (halfSize + 1)) - 1;
  uint32_t upperHalfMask = ~lowerHalfMask;

  // extracting upper lower part
  uint32_t a = (x & upperHalfMask) >> halfSize;
  uint32_t b = x & lowerHalfMask;
  uint32_t c = (y & upperHalfMask) >> halfSize;
  uint32_t d = y & lowerHalfMask;

  // computing steps
  uint32_t step1 = karatsubaConstExpr(a, c, halfSize);
  uint32_t step2 = karatsubaConstExpr(b, d, halfSize);
  uint32_t step3 = karatsubaConstExpr((a + b), (c + d), halfSize);
  uint32_t gauss = step3 - step2 - step1;
  return (step1 << (SIZE)) + (gauss << (halfSize)) + step2;
}

} // namespace algorithms
} // namespace cpp_tools
