#pragma once

#include <bitset>
#include <cmath>
#include <cstdint>
#include <iostream>

#ifdef MSVC
#include <intrin.h>
#endif
#ifdef CLANG
#include <x86intrin.h>
#endif

/**
 * Struct that allows us to easily access the different parts of the floating
 * point
 */
union SWFloat {
  struct {
    uint32_t mantissa : 23;
    uint32_t exponent : 8;
    uint32_t sign : 1;
  };
  float original;
};

inline uint32_t findHighestBit(uint32_t v) {
#ifdef MSVC
  return 31 - __lzcnt(v);
#endif
#ifdef CLANG
  return 31 - _lzcnt_u32(v);
#endif
}


inline uint32_t findHighestBitFromRight(uint32_t v) {
  // NOTE : the counter returned is counting from 0
  uint32_t r = 0;
  while (v >>= 1) {
    r++;
  }
  return r;
}

inline uint32_t findHighestBitLeft(uint32_t v) {

  // we need to handle the special case where the whole number is zero
  // if this is the case we set r to zero
  uint32_t r = (!v) ? 0 : 31;
  // we set v to -1, that will yield to 0xFFFFFFFF which will make
  // the loop exit at the first iteratin not touching r
  v = !v ? -1 : v;

  // here we keep checking if the bit from left to right is 0 or not
  while ((v & (1 << r)) == 0) {
    --r;
  }
  return r;
}

std::ostream &operator<<(std::ostream &os, const SWFloat &ff) {
  os << std::bitset<1>(ff.sign) << "-" << std::bitset<8>(ff.exponent) << "-"
     << std::bitset<23>(ff.mantissa);
  return os;
}

static const int DENORMAL = 111111111;
inline int normalize32BitMantissaInPlaceJumps(int &mantissa) {

  // NOTE : expected 32 bits mantissa with grs bits at the end
  // making a copy of the mantissa so we can freely manipulate it
  int tempMantissa = mantissa;

  // extracting the sticky bits
  int sticky = mantissa & 1;

  // here we find where the highest set bit is, keep in mind that find highest
  // bit returns the  numeration starting from 0, once we have that we need to
  // know how fare the highest bit is from the position it should be (26)
  int bit = 26 - findHighestBit(tempMantissa);

  // if the bit is greater than 23 we are out of range and we set everything as
  // denormal
  if (bit > 23) {
    mantissa = 0;
    return DENORMAL;
  }

  // at this point instead we have a good value to normalize
  // we shift left or right based on where the bit is
  if (bit > 0) {
    // tempMantissa = mantissa & ((1 << 26) - 1);
    // int bit = 26 - findHighestBit(tempMantissa);
    mantissa = mantissa << bit;
  } else {
    tempMantissa = mantissa & ((1 << 26) - 1);
    mantissa = mantissa >> abs(bit);
  }
  mantissa |= sticky;
  return bit;
}

inline uint32_t normalize32BitMantissaInPlace(uint32_t &mantissa) {

  // NOTE : expected 32 bits mantissa with grs bits at the end
  // making a copy of the mantissa so we can freely manipulate it
  uint32_t tempMantissa = mantissa;

  // extracting the sticky bits
  uint32_t sticky = mantissa & 1;

  // here we find where the highest set bit is, keep in mind that find highest
  // bit returns the  numeration starting from 0, once we have that we need to
  // know how fare the highest bit is from the position it should be (26)
  int bit = 26 - findHighestBit(tempMantissa);

  // need to be really careful here or MSVC won't generate instructions
  // conditional move, if here i use mantissa rather than tempMantissa it
  // wont work, clang is fine
  mantissa = bit > 23 ? 0 : tempMantissa;
  uint32_t returnValue = bit > 23 ? DENORMAL : bit;
  // at this point instead we have a good value to normalize
  // we shift left or right based on where the bit is

  uint32_t mantissaLeft = mantissa << bit;
  uint32_t mantissaRight = mantissa >> abs(bit);

  mantissa = (bit > 0) & (bit < 23) ? mantissaLeft : mantissaRight;
  mantissa |= sticky;
  return returnValue;
}

inline uint32_t insertHiddenOne(SWFloat value) {
  // casting the 23 bit mantissa to a 32 bit int
  // and inserting the 1 at the 24th slot
  uint32_t mantissa = (value.mantissa);
  return (mantissa |= (1 << 23));
}

inline uint32_t extendStickyGRSbits(uint32_t mantissa) { return mantissa << 3; }
inline uint64_t extendStickyGRSbits(uint64_t mantissa) { return mantissa << 3; }

inline uint32_t extractGRSbits(uint32_t mantissa) { return mantissa & 7u; }

inline uint32_t roundMantissa(uint32_t mantissa) {
  uint32_t grs = extractGRSbits(mantissa);
  // now that we extracted the values we can shift the mantissa by
  // the 3 extra bits
  uint32_t cleanedMantissa = mantissa >> 3;

  // now if the grs  is equal to 100 base two, aka 4 base 10,
  // we need to check if we need to round up, this happens only
  // if the LSB of the mantissa is one, we extract that with cleanedMantissa &1
  uint32_t toAdd = (grs == 4u) & (cleanedMantissa & 1) ? 1 : 0;

  // now we check just if the guard bit is one, and we did not already rounded
  // up thanks to the LSB we increase the rounding, this is basically computing
  // different branches of the if statement withouth banches

  // the first part of the check :((grs & 4u) == 4)
  // we make sure that the g bit is 1, we do that by making a mask 100 (which is
  // 4 in decimal, this mean only the 3rd bit will survive if set, the grs& 4u
  // operation will give us either a 4 or a 0. the second part we check if grs
  // is !=4 , meaning is not 100, because if so we alrady took care of it above
  toAdd += (((grs & 4u) == 4) & (grs != 4)) ? 1 : 0;
  // here we compute the mantissa if we have rounding, we are doing this
  // explicitelly because otherwise MSVC creates a jump instead clang doesnt
  // https://godbolt.org/g/StQajo
  uint32_t mantissaIf4u = cleanedMantissa + toAdd;

  // condtitional move based on guard bit

  cleanedMantissa = (grs & 4u) ? mantissaIf4u : cleanedMantissa;
  return cleanedMantissa;
}

inline int roundMantissaOneJump(int mantissa) {
  int grs = extractGRSbits(mantissa);
  // now that we extracted the values we can shift the mantissa by
  // the 3 extra bits
  int cleanedMantissa = mantissa >> 3;

  // if the guard bit was set we are going to try to round
  if ((grs & 4u)) {
    // in the case the grs configuration is 100, we only round up if the LSB of
    // the mantissa is 1
    int toAdd = (grs == 4) & (cleanedMantissa & 1) ? 1 : 0;
    toAdd += (grs == 4) ? 0 : 1;
    cleanedMantissa += toAdd;
  }
  return cleanedMantissa;
}

inline int roundMantissaTwoJump(int mantissa) {
  int grs = extractGRSbits(mantissa);
  // now that we extracted the values we can shift the mantissa by
  // the 3 extra bits
  int cleanedMantissa = mantissa >> 3;
  // if the guard bit was set we are going to try to round
  if ((grs & 4u)) {
    // in the case the grs configuration is 100, we only round up if the LSB of
    // the mantissa is 1
    if (grs == 4) {
      if (cleanedMantissa & 1) {
        cleanedMantissa += 1;
      }
    } else {
      cleanedMantissa += 1;
    }
  }
  return cleanedMantissa;
}

inline void shiftExponent(uint32_t &mantissa, uint32_t exponent) {
  uint32_t sticky = 0;

  for (uint32_t i = 0; i < exponent; ++i) {
    mantissa = mantissa >> 1;

    sticky |= mantissa & 1;
    mantissa |= sticky;
  }
}

inline void shiftExponent64(uint64_t &mantissa, int exponent) {
  uint64_t sticky = 0;
  for (int i = 0; i < exponent; ++i) {
    mantissa = mantissa >> 1;

    sticky |= mantissa & 1;
    mantissa |= sticky;
  }
}

inline uint32_t countMantissaBits(uint32_t mantissa) {
  while (!(mantissa & 1)) {
    mantissa = mantissa >> 1;
  }
  return findHighestBit(mantissa);
}

inline SWFloat swFloatAddition(SWFloat a, SWFloat b) {

  // the first step is to have both floating point on the
  // same exponents, once that is done we can perform the addition
  int deltaExponent = int(a.exponent) - int(b.exponent);

  // manipulating the mantissa to have the hidden one and the
  // grs bits for rounding, this will yield a 27 bits mantissa
  uint32_t amantissa = insertHiddenOne(a);
  uint32_t bmantissa = insertHiddenOne(b);
  amantissa = extendStickyGRSbits(amantissa);
  bmantissa = extendStickyGRSbits(bmantissa);

  // if the delta is negative, a has a lower exponent and needs
  // to be raised up, otherwise we do it to b
  if (deltaExponent < 0) {
    a.exponent = a.exponent + abs(deltaExponent);
    shiftExponent(amantissa, abs(deltaExponent));

  } else {
    b.exponent = b.exponent + abs(deltaExponent);
    shiftExponent(bmantissa, abs(deltaExponent));
  }

  SWFloat res;
  if (b.sign == a.sign) {
    // if the sign is the same we can just go ahead and perform the addition
    uint32_t addedMantissa = amantissa + bmantissa;
    int bit = normalize32BitMantissaInPlace(addedMantissa);
    addedMantissa = roundMantissa(addedMantissa);
    res.exponent = a.exponent - bit;
    res.mantissa = addedMantissa;
    res.sign = a.sign;

  } else {
    int mantissa = 0;
    int sign = 0;
    int exponent = a.exponent;

    // we compute both result based on if the aman or bman are neg
    // here we do a mult -1 rather than just - unary operator due to the
    // fact that visual studio was spittiong out a warning about still being
    // unsigned after the -
    int negativeA = (-1 * amantissa) + bmantissa;
    int negativeB = (-1 * bmantissa) + amantissa;
    // we use a conditional move to pick the corret computed mantissa
    mantissa = a.sign != 0 ? negativeA : negativeB;

    // checking if we need to flip the mantissa, if it is we hardcode the sign
    // to 1 and flip the mantissa
    bool shouldFlip = mantissa < 0;
    int flippedMantissa = -mantissa;
    mantissa = shouldFlip ? flippedMantissa : mantissa;
    sign = shouldFlip ? 1 : sign;

    uint32_t unsignedMantissa = static_cast<uint32_t>(mantissa);
    int bit = normalize32BitMantissaInPlace(unsignedMantissa);
    unsignedMantissa = roundMantissa(unsignedMantissa);
    // checking for denormal
    sign = (bit == DENORMAL) ? 0 : sign;
    exponent = (bit == DENORMAL) ? 0 : exponent - bit;

    // returning the result
    res.exponent = exponent;
    res.mantissa = unsignedMantissa;
    res.sign = sign;
  }

  return res;
}

inline uint64_t simpleMultFaster64(uint32_t a, uint32_t b) {

  // extending a and b to 64 bit
  uint64_t aextend = static_cast<uint64_t>(a);
  uint64_t bextend = static_cast<uint64_t>(b);

  // finding the highest bit that we will use to
  // loop only the necessary amount and not a fixed 32 bit
  uint32_t bbit = (findHighestBit(b));
  uint64_t result = 0;

  for (uint32_t bi = 0; bi <= bbit; ++bi) {
    uint64_t bmask = 1ll << bi;
    uint64_t bvalue = (bextend & bmask) >> bi;
    result += bvalue ? aextend << bi : 0;
  }

  return result;
}

SWFloat inline swFloatMultiplication(SWFloat a, SWFloat b) {

  uint32_t amantissa32 = insertHiddenOne(a);
  uint32_t bmantissa32 = insertHiddenOne(b);

  int aexp = int(a.exponent) - 127;
  int bexp = int(b.exponent) - 127;
  int exponent = (aexp + bexp) + 127;

  uint64_t mantissaMult = simpleMultFaster64(amantissa32, bmantissa32);
  uint64_t manSticky = extendStickyGRSbits(mantissaMult);
  shiftExponent64(manSticky, 46 - 23);

  uint32_t manSticky32 = static_cast<uint32_t>(manSticky);

  // normalizing and rounding
  int bit = normalize32BitMantissaInPlace(manSticky32);
  manSticky32 = roundMantissa(manSticky32);
  exponent -= bit;

  // here we shift by 3 again, since normalizing expects grs bits
  // we normalize again because there is the rare case where
  // the rounding might ripple all the way up and make the floating
  // point not normalized again (it actually came up during testing
  manSticky32 = manSticky32 << 3;
  bit = normalize32BitMantissaInPlace(manSticky32);
  manSticky32 = roundMantissa(manSticky32);
  exponent -= bit;

  SWFloat res;
  res.sign = (a.sign ^ b.sign) ? 1 : 0;
  res.exponent = exponent;
  res.mantissa = manSticky32;
  return res;
}

SWFloat inline swFloatDivision(SWFloat a, SWFloat b) {

  uint32_t amantissa = insertHiddenOne(a);
  uint32_t bmantissa = insertHiddenOne(b);

  // computing exponent
  int aexp = int(a.exponent) - 127;
  int bexp = int(b.exponent) - 127;
  int exponent = (aexp - bexp) + 127;

  uint32_t bmantbit = countMantissaBits(bmantissa) + 1;

  // generating mask
  uint32_t mask = ((1 << 24) - 1);
  uint32_t divmant = (bmantissa & mask) >> (24u - bmantbit);

  // now we need to extract the first nth bit from the amant
  uint64_t result = 0;
  uint32_t start = (amantissa >> (24u - bmantbit));

  int startingIndex = 24 - bmantbit - 1;
  for (int i = 0; i < (50); ++i, --startingIndex) {

    uint32_t currHigh = findHighestBit(start) + 1;
    if (currHigh >= bmantbit) {

      if (divmant <= start) {
        // it means we fit ay least once so we can subtract
        result = result << 1;
        result |= 1;
        start -= divmant;
      } else {
        result = result << 1;
      }
    } else {
      result = result << 1;
    }

    uint32_t extractedIfPositive = (amantissa >> abs(startingIndex)) & 1;
    uint32_t newExtracted = startingIndex >= 0 ? extractedIfPositive : 0;
    start = (start << 1) | newExtracted;
  }

  // extract extra res
  uint32_t extra = result & ((1 << 24) - 1);
  result = result >> 23;

  uint32_t result32 = static_cast<uint32_t>(result);
  uint32_t result32Sticky = result32 | 1;
  result32 = extra ? result32Sticky : result32;

  int bit = normalize32BitMantissaInPlace(result32);
  exponent -= bit;
  result32 = roundMantissa(result32);

  SWFloat res;
  res.sign = (a.sign ^ b.sign) ? 1 : 0;
  res.exponent = exponent;
  res.mantissa = result32;
  return res;
}
