//compile with g++ -std-c++11 -mavx2 -mbmi2 -O3 uv.cpp -o uvtest


#include <immintrin.h>
#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <chrono>

using namespace std;

constexpr float UV_OFFSET = 0.01f;
constexpr float UV_OFFSET_HALF = UV_OFFSET / 2.0f;
constexpr float UV_OFFSET_HALF_AVX = -UV_OFFSET / 2.0f;
const float off_buff[8] = { UV_OFFSET, 
                   UV_OFFSET_HALF_AVX, 
                   UV_OFFSET_HALF_AVX,
                   UV_OFFSET, 
                   UV_OFFSET_HALF_AVX,
                   UV_OFFSET_HALF_AVX, 
                   0.0f,0.0f};
const float storemask[8] = {0,0,0,0,0,0,1,1};



void offsetUVs(const float uv[2], float offset_uv[2])
{
    float u = uv[0];  
    float v = uv[1];  
    float w = 1.0f - u - v;
    if ( u < v && u <w)
    {
        offset_uv[0] = u + UV_OFFSET;
        offset_uv[1] = v - UV_OFFSET_HALF;
        return;
    }

    if ( v < u && v <w)
    {
        
        offset_uv[0] = u - UV_OFFSET_HALF;
        offset_uv[1] = v + UV_OFFSET;
        return;
    }

    offset_uv[0] = u - UV_OFFSET_HALF;
    offset_uv[1] = v - UV_OFFSET_HALF;

}

 
inline __m256 compress256(__m256 src, unsigned int mask /* from movmskps */)
{
    //mask is a interger on which each bit, represent wheter or not we should keep the result.
    //in my case I have a float[8], where 
    uint64_t expanded_mask = _pdep_u64(mask, 0x0101010101010101);  // unpack each bit to a byte
    expanded_mask *= 0xFF;    // mask |= mask<<1 | mask<<2 | ... | mask<<7;
    // ABC... -> AAAAAAAABBBBBBBBCCCCCCCC...: replicate each bit to fill its byte
    //
    // the identity shuffle for vpermps, packed to one index per byte
    const uint64_t identity_indices = 0x0706050403020100;    
    //extract on lower end the wanted bytes, basically removes and compact remaining
    //based on the mask, so the result should be the indices we need compacted
    //on the lower side, which will we use later to get the values we need
    uint64_t wanted_indices = _pext_u64(identity_indices, expanded_mask);

    //convertes 64 bits to a 128 register, zeroing out upper 64 register
    __m128i bytevec = _mm_cvtsi64_si128(wanted_indices);
    //expands a each byte to a 32 bit
    __m256i shufmask = _mm256_cvtepu8_epi32(bytevec);

    // 8-32 bit 
    return _mm256_permutevar8x32_ps(src, shufmask);
}

//https://godbolt.org/g/FYgupd
//https://godbolt.org/g/9I0H24
//http://stackoverflow.com/questions/36932240/avx2-what-is-the-most-efficient-way-to-pack-left-based-on-a-mask
void offsetUVsNoBranch2(const float uv[2], float offset_uv[2])
{
    //ref to make life easier should boil down to no op, compiler
    //will optimize it away
    const float& u = uv[0];
    const float& v = uv[1];

    __m128d uvd =  _mm_loadu_pd(reinterpret_cast<const double*>(uv));
    __m256d uvreg  = _mm256_broadcastsd_pd(uvd);
    
    __m256 offset = _mm256_loadu_ps(off_buff);
    __m256 to_be_masked =  _mm256_add_ps(_mm256_castpd_ps(uvreg),offset);

    //building the mask
    float w = 1.0f - (u + v);
    int isu = u<v && u<w;
    int isv = v<u && v <w;
    int isw = !(isu || isv);
    
    //extending bool and orring rather then mult and add?
    unsigned int m =3*isu + 12*isv + 48*isw;
    __m256 res = compress256(to_be_masked, m);
    
    __m256i storemaskreg =  _mm256_loadu_si256(reinterpret_cast<const __m256i*>(storemask));
    _mm256_maskstore_ps(offset_uv,storemaskreg,res);
}

int main()
{
    float uv[2];
    float uvoff[2];
    const uint32_t ITERATIONS = 1000000;
    
    auto br1 = chrono::high_resolution_clock::now();
    for (uint32_t i=0; i< ITERATIONS; ++i)
    {
        uv[0] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        uv[1] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        offsetUVs(uv,uvoff);
    }
    auto br2 = chrono::high_resolution_clock::now();
    auto brt = chrono::duration_cast<chrono::microseconds>(br2 - br1).count();
    std::cout <<"with branch: "<< brt<<" micro" <<std::endl;



    auto nbr1 = chrono::high_resolution_clock::now();
    for (uint32_t i=0; i< ITERATIONS; ++i)
    {
        uv[0] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        uv[1] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        offsetUVsNoBranch2(uv,uvoff);
    }
    auto nbr2 = chrono::high_resolution_clock::now();
    auto nbrt = chrono::duration_cast<chrono::microseconds>(nbr2 - nbr1).count();
    std::cout << "branchless:  "<<nbrt<<" micro"<<std::endl;
    std::cout << "difference: " <<(int) ((1.0f - float(nbrt)/float(brt))*100.0f)<<"%"<<std::endl;
    return 0;
}
