#define CPPECC_IMPLEMENTATION
#include "cppecc.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static const cppecc_u64 Increment = 0x14057B7EF767814FULL;
static const cppecc_u64 Multiplier = 0x5851F42D4C957F2DULL;
static cppecc_u64 state_;

cppecc_u32 rotr32(cppecc_u32 x, cppecc_u32 r)
{
    return (x >> r) | (x << ((~r + 1) & 31U));
}

float toF32_1(cppecc_u32 x)
{
    static const cppecc_u32 m0 = 0x3F800000U;
    static const cppecc_u32 m1 = 0x007FFFFFU;
    x = m0 | (x & m1);
    return (*(float*)&x) - 1.000000000f;
}

void srand_(cppecc_u64 seed)
{
    state_ = Increment + seed;
}

cppecc_u32 urand()
{
    cppecc_u64 x = state_;
    cppecc_u32 count = (cppecc_u32)(x >> 59);
    state_ = x * Multiplier + Increment;
    x ^= x >> 18;
    return rotr32((cppecc_u32)(x >> 27), count);
}

float frand()
{
    return toF32_1(urand());
}

struct Result
{
    long long encodeTime_;
    long long decodeTime_;
    cppecc_s32 numErrors_;
    cppecc_s64 avgErrors_;
};

cppecc_s32 random_range(cppecc_s32 minx, cppecc_s32 maxx)
{
    return (cppecc_s32)(frand()*(maxx-minx))+minx;
}

struct timespec get_micro()
{
    struct timespec clock;
    clock_gettime(CLOCK_MONOTONIC, &clock);
    return clock;
}

long long duration_micro(struct timespec start, struct timespec end)
{
    long long nsec_start = start.tv_nsec;
    long long nsec_end = start.tv_nsec;
    if(nsec_end<nsec_start){
        nsec_end += (1000LL * 1000LL * 1000LL);
        end.tv_sec -= 1;
    }
    long long sec = end.tv_sec - end.tv_sec;
    long long nsec = nsec_end - nsec_start;
    return sec*(1000LL*1000LL) + nsec/(1000LL);
}

struct Result reed_solomon(cppecc_s32 messageSize, cppecc_s32 eccSize, cppecc_s32 maxErrors, cppecc_s32 count)
{
    const cppecc_s32 MaxECC = eccSize>>1;
    cppecc_u8* message = (cppecc_u8*)malloc(messageSize);
    cppecc_u8* encoded = (cppecc_u8*)malloc(messageSize+eccSize);
    cppecc_u8* decoded = (cppecc_u8*)malloc(messageSize+eccSize);
    struct RSContext context;
    struct Result result = {0};

    gf_initialize(&context, eccSize);

    long long encodeTime = 0;
    long long decodeTime = 0;
    struct timespec start;
    struct timespec end;
    for(cppecc_s32 i = 0; i < count; ++i) {
        for(cppecc_s32 j=0; j<messageSize; ++j){
            message[j] = (cppecc_u8)(urand()&0xFFU);
        }
        memcpy(encoded, message, messageSize);
        start = get_micro();
        rs_encode(&context, messageSize, &encoded[0], eccSize);
        end = get_micro();
        encodeTime += duration_micro(start, end);

        memcpy(decoded, encoded, messageSize+eccSize);

        cppecc_s32 numErrors = random_range(0, maxErrors);
        result.avgErrors_ += numErrors;
        for(cppecc_s32 j=0; j<numErrors; ++j){
            decoded[random_range(0, messageSize+eccSize-1)] ^= urand();
        }
        numErrors = 0;
        for(cppecc_s32 j=0; j<(messageSize+eccSize); ++j){
            if(encoded[j] != decoded[j]){
                ++numErrors;
            }
            encoded[j] ^= decoded[j];
        }

        start = get_micro();
        cppecc_s32 corrected = rs_decode(&context, messageSize, &decoded[0], eccSize);
        end = get_micro();
        decodeTime += duration_micro(start, end);

        if(numErrors<=MaxECC){
            if(corrected<0){
                result.numErrors_ += 1;
            }else{
                for(cppecc_s32 j=0; j<messageSize; ++j){
                    if(message[j] != decoded[j]){
                        printf("message size: %d, ecc size: %d\n", messageSize, eccSize);
                        assert(false);
                    }
                }
            }
        }else{
            if(0<=corrected){
                result.numErrors_ += 1;
            }
        }
    }

    result.encodeTime_ = encodeTime / count;
    result.decodeTime_ = decodeTime / count;
    result.avgErrors_ /= count;
    free(decoded);
    free(encoded);
    free(message);
    return result;
}

int main(void)
{
    static const int Patterns = 1024;
    static const int Count = 4096;
    int MinMessageSize = 16;
    srand_(time(NULL));

    struct Result result;
    for(cppecc_s32 i = 0; i < Patterns; ++i) {
        cppecc_s32 messageSize = random_range(MinMessageSize, CPPECC_GF_NW1-CPPECC_MAX_ECC_SIZE);
        cppecc_s32 eccMaxSize = (CPPECC_MAX_ECC_SIZE<(messageSize-1))? CPPECC_MAX_ECC_SIZE : (messageSize-1);
        cppecc_s32 eccSize = random_range(1, eccMaxSize);
        result = reed_solomon(messageSize, eccSize, eccSize, Count);
        printf("result: message size: %d ecc size: %d encode time (micro): %lld decode time (micro): %lld errors: %d/%ld\n", messageSize, eccSize, result.encodeTime_, result.decodeTime_, result.numErrors_, result.avgErrors_);
    }
    return 0;
}
