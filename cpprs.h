#ifndef INC_CPPRS_H_
#define INC_CPPRS_H_
/**
@file cpprs.h
@author t-sakai

# License
This software is distributed under two licenses, choose whichever you like.

## MIT License
Copyright (c) 2020 Takuro Sakai

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Public Domain
This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org>
```

# Abount Reed-Solomon Codes
https://www.cs.cmu.edu/~guyb/realworld/cpprs/reed_solomon_codes.html
https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders

# Usage
Put '#define CPPRS_IMPLEMENTATION' before including this file to create the implementation.
*/
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#ifndef __cplusplus
#    include <stdbool.h>
#endif

#ifdef __cplusplus
#    define CPPRS_NAMESPACE_BEGIN(name) \
        namespace name \
        {
#    define CPPRS_NAMESPACE_END(name) }
#    define CPPRS_NAMESPACE_EMPTY_BEGIN \
        namespace \
        {
#    define CPPRS_NAMESPACE_EMPTY_END }
#    define CPPRS_STATIC
#    define CPPRS_REINTERPRET_CAST(type) reinterpret_cast<type>
#    define CPPRS_STATIC_CAST(type) static_cast<type>
#else
#    define CPPRS_NAMESPACE_BEGIN(name)
#    define CPPRS_NAMESPACE_END(name)
#    define CPPRS_NAMESPACE_EMPTY_BEGIN
#    define CPPRS_NAMESPACE_EMPTY_END
#    define CPPRS_STATIC static
#    define CPPRS_REINTERPRET_CAST(type) (type)
#    define CPPRS_STATIC_CAST(type) (type)
#endif

CPPRS_NAMESPACE_BEGIN(cpprs)

#ifndef CPPRS_NULL
#    ifdef __cplusplus
#        if 201103L <= __cplusplus || 1700 <= _MSC_VER
#            define CPPRS_NULL nullptr
#        else
#            define CPPRS_NULL 0
#        endif
#    else
#        define CPPRS_NULL (void*)0
#    endif
#endif

#ifndef CPPRS_TYPES
#    define CPPRS_TYPES
typedef int8_t rs_s8;
typedef int16_t rs_s16;
typedef int32_t rs_s32;
typedef int64_t rs_s64;

typedef uint8_t rs_u8;
typedef uint16_t rs_u16;
typedef uint32_t rs_u32;
typedef uint64_t rs_u64;

typedef float rs_f32;
typedef double rs_f64;

typedef size_t rs_size_t;
#endif // CPPRS_TYPES

#ifndef CPPRS_ASSERT
#    define CPPRS_ASSERT(exp) assert(exp)
#endif

#ifdef __cplusplus
static const rs_u32 RS_GF_W = 8;
static const rs_u32 RS_GF_NW = (1 << RS_GF_W);
static const rs_u32 RS_GF_NW1 = RS_GF_NW - 1;

//static const rs_u32 RS_GF_PRIMITIVE = 0x11DU;

static const rs_s32 RS_MAX_BODY_SIZE = RS_GF_NW1;
static const rs_s32 RS_MAX_ECC_SIZE = 52;
static const rs_s32 RS_MAX_ECC_RATE = 10;

static const rs_s32 RS_ERROR = -1;

#else
#    define RS_GF_W (8)
#    define RS_GF_NW (1 << RS_GF_W)
#    define RS_GF_NW1 (RS_GF_NW - 1)
//#define RS_GF_PRIMITIVE (0x11DU)

#    define RS_MAX_BODY_SIZE (RS_GF_NW1)
#    define RS_MAX_ECC_SIZE (52)
#    define RS_MAX_ECC_RATE (10)

#    define RS_ERROR (-1)
#endif

struct RSContext
{
    rs_u8 generator_[RS_MAX_ECC_SIZE + 1];
    rs_u8 syndromes_[RS_MAX_ECC_SIZE + 1];
    rs_u8 sigma_[RS_MAX_ECC_SIZE];
    rs_u8 errorPositions_[RS_MAX_ECC_SIZE];
    rs_u8 omega_[RS_MAX_ECC_SIZE + RS_MAX_ECC_SIZE];

    rs_u8 temp0_[RS_GF_NW];
};

/**
 @brief Initialize the generation polynomial from numSymbols.
 @param [in,out] context
 @param [in] numSymbols
 */
void gf_initialize(RSContext* context, rs_s32 numSymbols);

rs_u8 gf_add(rs_u8 a, rs_u8 b);

rs_u8 gf_mul(rs_u8 a, rs_u8 b);

rs_u8 gf_div(rs_u8 a, rs_u8 b);

rs_u8 gf_pow(rs_u8 x, rs_s32 p);

rs_u8 gf_inverse(rs_u8 x);

void gf_poly_scale(rs_s32 size, rs_u8 result[], const rs_u8 p[], rs_u8 x);

rs_s32 gf_poly_add(rs_u8 result[], rs_s32 psize, const rs_u8 p[], rs_s32 qsize, const rs_u8 q[]);

rs_s32 gf_poly_mul(rs_u8 result[], rs_s32 psize, const rs_u8 p[], rs_s32 qsize, const rs_u8 q[]);

rs_u8 gf_poly_eval(rs_s32 size, const rs_u8 poly[], rs_u8 x);

rs_s32 gf_poly_div(rs_u8 result[], rs_s32 sizeDividend, const rs_u8 dividend[], rs_s32 sizeDivisor, const rs_u8 divisor[]);

void rs_generator_poly(rs_s32 size, rs_u8 result[], rs_u8 tmp[]);

rs_s32 rs_modified_berlekamp_massey(RSContext* context, rs_u8 result[], rs_s32 numSyndromes, const rs_u8 syndromes[]);
rs_s32 rs_chien_search(rs_u8 result[2], rs_u8 start, rs_u8 end, rs_u8 a, rs_u8 b);
rs_s32 rs_chien_search(rs_u8 result[], rs_u8 size, rs_u8 numSigma, const rs_u8 sigma[]);
void rs_error_correct_forney(rs_u8 result[], rs_s32 length, rs_s32 numErrors, const rs_u8 pos[], rs_s32 numSigma, const rs_u8 sigma[], rs_s32 numOmega, const rs_u8 omega[]);

/**
 @brief Add redundant data to the original message for error correction.
 @param [in, out] context ... Used for the generation polynomial and buffers.
 @param [in] size ... message size
 @param [in, out] message[] ... The size should be 'size + numSymbols'. Output will be made from the original message and added redundant symbols (that size is numSymbols).
 @param numSymbols ... size of redundant symbols, that is equivalent to capability of error corrections.
 */
void rs_encode(RSContext* context, rs_s32 size, rs_u8 message[], rs_s32 numSymbols);

/**
 @brief Try to recover the original message from redundant symbols. But, when the number of errors exceeds the capability of Reed-Solomon codes, the message never be recoverted correctly.
 @param [in, out] context ... Used for buffers
 @param [in] size ... message size
 @param [in, out] message[] ... The size should be 'size + numSymbols'. Output's corrupted symbols will be corrected.
 @param numSymbols ... size of redundant symbols, that is equivalent to capability of error corrections.
 @return The number of corrected symbols. When the number of erros exceeds the capability (so, it's numSymbols/2), the original message will not be recovered.
 */
rs_s32 rs_decode(RSContext* context, rs_s32 size, rs_u8 message[], rs_s32 numSymbols);

CPPRS_NAMESPACE_END(cpprs)
#endif // INC_CPPRS_H_

#ifdef CPPRS_IMPLEMENTATION
CPPRS_NAMESPACE_BEGIN(cpprs)

CPPRS_NAMESPACE_EMPTY_BEGIN

// clang-format off
static const rs_u8 gflog[RS_GF_NW]={
    0x0U,0x0U,0x1U,0x19U,0x2U,0x32U,0x1AU,0xC6U,0x3U,0xDFU,0x33U,0xEEU,0x1BU,0x68U,0xC7U,0x4BU,
    0x4U,0x64U,0xE0U,0xEU,0x34U,0x8DU,0xEFU,0x81U,0x1CU,0xC1U,0x69U,0xF8U,0xC8U,0x8U,0x4CU,0x71U,
    0x5U,0x8AU,0x65U,0x2FU,0xE1U,0x24U,0xFU,0x21U,0x35U,0x93U,0x8EU,0xDAU,0xF0U,0x12U,0x82U,0x45U,
    0x1DU,0xB5U,0xC2U,0x7DU,0x6AU,0x27U,0xF9U,0xB9U,0xC9U,0x9AU,0x9U,0x78U,0x4DU,0xE4U,0x72U,0xA6U,
    0x6U,0xBFU,0x8BU,0x62U,0x66U,0xDDU,0x30U,0xFDU,0xE2U,0x98U,0x25U,0xB3U,0x10U,0x91U,0x22U,0x88U,
    0x36U,0xD0U,0x94U,0xCEU,0x8FU,0x96U,0xDBU,0xBDU,0xF1U,0xD2U,0x13U,0x5CU,0x83U,0x38U,0x46U,0x40U,
    0x1EU,0x42U,0xB6U,0xA3U,0xC3U,0x48U,0x7EU,0x6EU,0x6BU,0x3AU,0x28U,0x54U,0xFAU,0x85U,0xBAU,0x3DU,
    0xCAU,0x5EU,0x9BU,0x9FU,0xAU,0x15U,0x79U,0x2BU,0x4EU,0xD4U,0xE5U,0xACU,0x73U,0xF3U,0xA7U,0x57U,
    0x7U,0x70U,0xC0U,0xF7U,0x8CU,0x80U,0x63U,0xDU,0x67U,0x4AU,0xDEU,0xEDU,0x31U,0xC5U,0xFEU,0x18U,
    0xE3U,0xA5U,0x99U,0x77U,0x26U,0xB8U,0xB4U,0x7CU,0x11U,0x44U,0x92U,0xD9U,0x23U,0x20U,0x89U,0x2EU,
    0x37U,0x3FU,0xD1U,0x5BU,0x95U,0xBCU,0xCFU,0xCDU,0x90U,0x87U,0x97U,0xB2U,0xDCU,0xFCU,0xBEU,0x61U,
    0xF2U,0x56U,0xD3U,0xABU,0x14U,0x2AU,0x5DU,0x9EU,0x84U,0x3CU,0x39U,0x53U,0x47U,0x6DU,0x41U,0xA2U,
    0x1FU,0x2DU,0x43U,0xD8U,0xB7U,0x7BU,0xA4U,0x76U,0xC4U,0x17U,0x49U,0xECU,0x7FU,0xCU,0x6FU,0xF6U,
    0x6CU,0xA1U,0x3BU,0x52U,0x29U,0x9DU,0x55U,0xAAU,0xFBU,0x60U,0x86U,0xB1U,0xBBU,0xCCU,0x3EU,0x5AU,
    0xCBU,0x59U,0x5FU,0xB0U,0x9CU,0xA9U,0xA0U,0x51U,0xBU,0xF5U,0x16U,0xEBU,0x7AU,0x75U,0x2CU,0xD7U,
    0x4FU,0xAEU,0xD5U,0xE9U,0xE6U,0xE7U,0xADU,0xE8U,0x74U,0xD6U,0xF4U,0xEAU,0xA8U,0x50U,0x58U,0xAFU,
};
// clang-format on

// clang-format off
static const rs_u8 gfexp[RS_GF_NW]={
    0x1U,0x2U,0x4U,0x8U,0x10U,0x20U,0x40U,0x80U,0x1DU,0x3AU,0x74U,0xE8U,0xCDU,0x87U,0x13U,0x26U,
    0x4CU,0x98U,0x2DU,0x5AU,0xB4U,0x75U,0xEAU,0xC9U,0x8FU,0x3U,0x6U,0xCU,0x18U,0x30U,0x60U,0xC0U,
    0x9DU,0x27U,0x4EU,0x9CU,0x25U,0x4AU,0x94U,0x35U,0x6AU,0xD4U,0xB5U,0x77U,0xEEU,0xC1U,0x9FU,0x23U,
    0x46U,0x8CU,0x5U,0xAU,0x14U,0x28U,0x50U,0xA0U,0x5DU,0xBAU,0x69U,0xD2U,0xB9U,0x6FU,0xDEU,0xA1U,
    0x5FU,0xBEU,0x61U,0xC2U,0x99U,0x2FU,0x5EU,0xBCU,0x65U,0xCAU,0x89U,0xFU,0x1EU,0x3CU,0x78U,0xF0U,
    0xFDU,0xE7U,0xD3U,0xBBU,0x6BU,0xD6U,0xB1U,0x7FU,0xFEU,0xE1U,0xDFU,0xA3U,0x5BU,0xB6U,0x71U,0xE2U,
    0xD9U,0xAFU,0x43U,0x86U,0x11U,0x22U,0x44U,0x88U,0xDU,0x1AU,0x34U,0x68U,0xD0U,0xBDU,0x67U,0xCEU,
    0x81U,0x1FU,0x3EU,0x7CU,0xF8U,0xEDU,0xC7U,0x93U,0x3BU,0x76U,0xECU,0xC5U,0x97U,0x33U,0x66U,0xCCU,
    0x85U,0x17U,0x2EU,0x5CU,0xB8U,0x6DU,0xDAU,0xA9U,0x4FU,0x9EU,0x21U,0x42U,0x84U,0x15U,0x2AU,0x54U,
    0xA8U,0x4DU,0x9AU,0x29U,0x52U,0xA4U,0x55U,0xAAU,0x49U,0x92U,0x39U,0x72U,0xE4U,0xD5U,0xB7U,0x73U,
    0xE6U,0xD1U,0xBFU,0x63U,0xC6U,0x91U,0x3FU,0x7EU,0xFCU,0xE5U,0xD7U,0xB3U,0x7BU,0xF6U,0xF1U,0xFFU,
    0xE3U,0xDBU,0xABU,0x4BU,0x96U,0x31U,0x62U,0xC4U,0x95U,0x37U,0x6EU,0xDCU,0xA5U,0x57U,0xAEU,0x41U,
    0x82U,0x19U,0x32U,0x64U,0xC8U,0x8DU,0x7U,0xEU,0x1CU,0x38U,0x70U,0xE0U,0xDDU,0xA7U,0x53U,0xA6U,
    0x51U,0xA2U,0x59U,0xB2U,0x79U,0xF2U,0xF9U,0xEFU,0xC3U,0x9BU,0x2BU,0x56U,0xACU,0x45U,0x8AU,0x9U,
    0x12U,0x24U,0x48U,0x90U,0x3DU,0x7AU,0xF4U,0xF5U,0xF7U,0xF3U,0xFBU,0xEBU,0xCBU,0x8BU,0xBU,0x16U,
    0x2CU,0x58U,0xB0U,0x7DU,0xFAU,0xE9U,0xCFU,0x83U,0x1BU,0x36U,0x6CU,0xD8U,0xADU,0x47U,0x8EU,0xCCU,
};
// clang-format off

CPPRS_NAMESPACE_EMPTY_END

//
void gf_initialize(RSContext* context, rs_s32 numSymbols)
{
    //const rs_u32 prim = GF_PRIMITIVE;
    //const rs_u32 nw1 = RS_GF_NW1;
    //gflog[0] = 0;
    //rs_u32 b=1;
    //for(rs_u32 log=0; log<nw1; ++log){
    //    gflog[b] = log;
    //    gfexp[log] = b;
    //    b <<= 1;
    //    if(0 != (b & RS_GF_NW)){
    //        b = (b ^ prim) & RS_GF_NW1;
    //    }
    //}
    rs_generator_poly(numSymbols, context->generator_, context->temp0_);
}

rs_u8 gf_add(rs_u8 a, rs_u8 b)
{
    return a ^ b;
}

//
rs_u8 gf_mul(rs_u8 a, rs_u8 b)
{
    if(0 == a || 0 == b) {
        return 0;
    }
    rs_u32 sum = CPPRS_STATIC_CAST(rs_s32)(gflog[a]) + gflog[b];
    if(RS_GF_NW1 <= sum) {
        sum -= RS_GF_NW1;
    }
    return gfexp[sum];
}

rs_u8 gf_mulexp(rs_u8 a, rs_u8 b)
{
    if(0 == a) {
        return 0;
    }
    rs_u32 sum = CPPRS_STATIC_CAST(rs_s32)(gflog[a]) + b;
    if(RS_GF_NW1 <= sum) {
        sum -= RS_GF_NW1;
    }
    return gfexp[sum];
}

//
rs_u8 gf_div(rs_u8 a, rs_u8 b)
{
    if(0 == a) {
        return 0;
    }
    if(0 == b) {
        return CPPRS_STATIC_CAST(rs_u8)(-1);
    }
    rs_s32 diff = CPPRS_STATIC_CAST(rs_s32)(gflog[a]) - gflog[b];
    if(diff < 0) {
        diff += RS_GF_NW1;
    }
    return gfexp[diff];
}

rs_u8 gf_divexp(rs_u8 a, rs_u8 b)
{
    if(0 == a) {
        return 0;
    }
    rs_s32 diff = CPPRS_STATIC_CAST(rs_s32)(gflog[a]) - b;
    if(diff < 0) {
        diff += RS_GF_NW1;
    }
    return gfexp[diff];
}

rs_u8 gf_pow(rs_u8 x, rs_s32 p)
{
    return gfexp[(gflog[x] * p) % RS_GF_NW1];
}

rs_u8 gf_inverse(rs_u8 x)
{
    return gfexp[RS_GF_NW1 - gflog[x]];
}

void gf_poly_scale(rs_s32 size, rs_u8 result[], const rs_u8 p[], rs_u8 x)
{
    for(rs_s32 i = 0; i < size; ++i) {
        result[i] = gf_mul(p[i], x);
    }
}

rs_s32 gf_poly_add(rs_u8 result[], rs_s32 psize, const rs_u8 p[], rs_s32 qsize, const rs_u8 q[])
{
    rs_s32 size = (psize < qsize) ? qsize : psize;
    rs_s32 d = (psize < qsize) ? qsize - psize : 0;
    for(rs_s32 i = 0; i < d; ++i) {
        result[i] = 0;
    }

    for(rs_s32 i = 0; i < psize; ++i) {
        result[i + size - psize] = p[i];
    }

    for(rs_s32 i = 0; i < qsize; ++i) {
        result[i + size - qsize] ^= q[i];
    }
    return size;
}

rs_s32 gf_poly_mul(rs_u8 result[], rs_s32 psize, const rs_u8 p[], rs_s32 qsize, const rs_u8 q[])
{
    rs_s32 total = psize + qsize - 1;
    for(rs_s32 i = 0; i < total; ++i) {
        result[i] = 0;
    }
    for(rs_s32 i = 0; i < qsize; ++i) {
        for(rs_s32 j = 0; j < psize; ++j) {
            result[i + j] ^= gf_mul(p[j], q[i]);
        }
    }
    return total;
}

rs_s32 gf_poly_mul(rs_u8 result[], rs_s32 psize, const rs_u8 p[], rs_s32 qsize, const rs_u8 q[], rs_s32 l)
{
    rs_s32 total = l;
    for(rs_s32 i = 0; i < total; ++i) {
        result[i] = 0;
    }
    psize = (psize < l) ? psize : l;
    for(rs_s32 i = 0; i < psize; ++i) {
        if(0 == p[i]) {
            continue;
        }
        rs_u8 logp = gflog[p[i]];
        rs_s32 qs = l - i;
        qs = qsize < qs ? qsize : qs;
        for(rs_s32 j = 0; j < qs; ++j) {
            if(0 == q[j]) {
                continue;
            }
            result[i + j] ^= gf_mulexp(q[j], logp);
        }
    }
    return total;
}

rs_u8 gf_poly_eval(rs_s32 size, const rs_u8 poly[], rs_u8 x)
{
    rs_u8 y = poly[0];
    for(rs_s32 i = 1; i < size; ++i) {
        y = gf_mul(y, x) ^ poly[i];
    }
    return y;
}

#if 0
rs_s32 gf_poly_div(rs_u8 result[], rs_s32 sizeDividend, const rs_u8 dividend[], rs_s32 sizeDivisor, const rs_u8 divisor[])
{
    for(rs_s32 i=0; i<sizeDividend; ++i){
        result[i] = dividend[i];
    }
    rs_s32 size = sizeDividend - sizeDivisor + 1;
    for(rs_s32 i=0; i<size; ++i){
        rs_u8 coef = result[i];
        if( 0 != coef){
            for(rs_s32 j=1; j<sizeDivisor; ++j){
                if(0 != divisor[j]){
                    result[i+j] ^= gf_mul(divisor[j], coef);
                }
            }
        }
    }

    rs_s32 sizeResult = sizeDividend - size;
    for(rs_s32 i=0; i<sizeResult; ++i){
        result[i] = result[i+size];
    }
    return sizeResult;
}
#endif

rs_u8 gf_omega_value(rs_s32 numOmega, const rs_u8 omega[], rs_u8 l)
{
    rs_u8 w = l;
    rs_u8 o = omega[0];
    for(rs_s32 i = 1; i < numOmega; i++) {
        o ^= gf_mulexp(omega[i], w);
        w = (w + l) % RS_GF_NW1;
    }
    return o;
}

rs_u8 gf_sigma_dash_value(rs_s32 numSigma, const rs_u8 sigma[], rs_u8 l)
{
    rs_s32 size = numSigma - 1;
    rs_u8 l2 = (2 * l) % RS_GF_NW1;
    rs_u8 w = l2;
    rs_u8 d = sigma[1];
    for(rs_s32 i = 3; i <= size; i += 2) {
        d ^= gf_mulexp(sigma[i], w);
        w = (w + l2) % RS_GF_NW1;
    }
    return d;
}

void rs_generator_poly(rs_s32 size, rs_u8 result[], rs_u8 tmp[])
{
    result[0] = 1;
    rs_u8 poly[2] = {1, 0};
    rs_s32 polySize = 1;
    for(rs_s32 i = 0; i < size; ++i) {
        poly[1] = gfexp[i];
        polySize = gf_poly_mul(tmp, polySize, result, 2, poly);
        for(rs_s32 j = 0; j < polySize; ++j) {
            result[j] = tmp[j];
        }
    }
}

void rs_encode(RSContext* context, rs_s32 size, rs_u8 message[], rs_s32 numSymbols)
{
    CPPRS_ASSERT(CPPRS_NULL != context);
    CPPRS_ASSERT(CPPRS_STATIC_CAST(rs_u32)(size + numSymbols) < RS_GF_NW);
    CPPRS_ASSERT(numSymbols <= RS_MAX_ECC_SIZE);

    rs_u8* result = context->temp0_;
    for(rs_s32 i = 0; i < size; ++i) {
        result[i] = message[i];
    }
    rs_u8* generator = context->generator_;
    rs_s32 total = size + numSymbols;
    for(rs_s32 i = size; i < total; ++i) {
        result[i] = 0;
    }

    for(rs_s32 i = 0; i < size; ++i) {
        if(0 != result[i]) {
            for(rs_s32 j = 1; j <= numSymbols; ++j) {
                result[i + j] ^= gf_mul(generator[j], result[i]);
            }
        }
    }
    for(rs_s32 i = size; i < total; ++i) {
        message[i] = result[i];
    }
}

#if 0
bool rs_has_erros(rs_s32 size, const rs_u8 syndromes[])
{
    for(rs_s32 i = 1; i <= size; ++i) {
        if(0 != syndromes[i]) {
            return true;
        }
    }
    return false;
}

void rs_forney_syndromes(rs_u8 forney_syndromes[], rs_s32 message_size, rs_s32 numSyndromes, const rs_u8 syndromes[], rs_s32 numErasePos, const rs_u8 erasePos[])
{
    for(rs_s32 i = 0; i < numSyndromes; ++i) {
        forney_syndromes[i] = syndromes[i];
    }
    for(rs_s32 i = 0; i < numErasePos; ++i) {
        rs_u8 x = gfexp[message_size - 1 - erasePos[i]];
        for(rs_s32 j = 0; j < numSyndromes - 1; ++j) {
            forney_syndromes[j] = gf_mul(forney_syndromes[j], x) ^ forney_syndromes[j + 1];
        }
    }
}
#endif

rs_s32 rs_modified_berlekamp_massey(RSContext* context, rs_u8 result[], rs_s32 numSyndromes, const rs_u8 syndromes[])
{
    rs_u8 b0[RS_MAX_ECC_SIZE + 1] = {
        0,
        1,
    };
    rs_u8 b1[RS_MAX_ECC_SIZE + 1] = {
        1,
    };

    rs_u8* sg0 = b0;
    rs_u8* sg1 = b1;
    rs_u8* work = context->temp0_;
    for(rs_s32 i = 0; i <= numSyndromes; ++i) {
        work[i] = 0;
    }

    rs_s32 s0 = 1;
    rs_s32 s1 = 0;

    rs_s32 k = -1;

    for(rs_s32 i = 0; i < numSyndromes; ++i) {
        rs_s32 s = syndromes[i];
        for(rs_s32 j = 1; j <= s1; ++j) {
            s ^= gf_mul(sg1[j], syndromes[i - j]);
        }
        if(0 != s) {
            rs_u8 l = gflog[s];
            for(rs_s32 j = 0; j <= i; ++j) {
                work[j] = sg1[j] ^ gf_mulexp(sg0[j], l);
            }
            rs_s32 d = i - k;
            if(s1 < d) {
                for(rs_s32 j = 0; j <= s0; ++j) {
                    sg0[j] = gf_divexp(sg1[j], l);
                }
                k = i - s1;
                s0 = d;
                s1 = d;
            }
            rs_u8* tmp = sg1;
            sg1 = work;
            work = tmp;
        } //if(0 != s)
        for(rs_s32 j = s0 - 1; 0 <= j; --j) {
            sg0[j + 1] = sg0[j];
        }
        sg0[0] = 0;
        ++s0;
    } //for(rs_s32 i=0

    if(0 == sg1[s1]) {
        return -1;
    }
    rs_s32 size = s1 + 1;
    for(rs_s32 i = 0; i < size; ++i) {
        result[i] = sg1[i];
    }
    return size;
}

rs_s32 rs_chien_search(rs_u8 result[2], rs_u8 start, rs_u8 end, rs_u8 a, rs_u8 b)
{
    for(rs_u8 i = start; i < end; ++i) {
        rs_u8 z0 = gfexp[i];
        rs_u8 z1 = CPPRS_STATIC_CAST(rs_u8)(a) ^ z0;

        if(b == gf_mulexp(z1, i)) {
            rs_u8 index = gflog[z1];
            if(index <= i || end <= index) {
                return -1;
            }
            result[0] = z1;
            result[1] = z0;
            return 2;
        }
    }
    return -1;
}

rs_s32 rs_chien_search(rs_u8 result[], rs_u8 size, rs_u8 numSigma, const rs_u8 sigma[])
{
    rs_u8 s0 = numSigma - 1;
    rs_u8 sum = sigma[1];
    rs_u8 mul = sigma[s0];
    if(1 == s0) {
        if(size <= gflog[sum]) {
            return -1;
        }
        result[0] = sum;
        return 1;
    }
    if(2 == s0) {
        return rs_chien_search(result, 0, size, sum, mul);
    }

    rs_u8 temp0[4];
    rs_u8 index = s0 - 1;
    for(rs_u8 i = 0, z = RS_GF_NW1; i < size; ++i, --z) {
        rs_u8 temp = 1;
        for(rs_u8 j = 1, wz = z; j <= s0; ++j, wz = (wz + z) % RS_GF_NW1) {
            temp ^= gf_mulexp(sigma[j], wz);
        }
        if(0 != temp) {
            continue;
        }
        rs_u8 p = gfexp[i];
        sum ^= p;
        mul = gf_div(mul, p);
        result[index--] = p;
        if(1 == index) {
            rs_s32 t = rs_chien_search(temp0, i + 1, size, sum, mul);
            if(t < 0) {
                return -1;
            }
            result[0] = temp0[0];
            result[1] = temp0[1];
            return s0;
        }
    }
    return -1;
}

void rs_error_correct_forney(rs_u8 result[], rs_s32 length, rs_s32 numErrors, const rs_u8 pos[], rs_s32 numSigma, const rs_u8 sigma[], rs_s32 numOmega, const rs_u8 omega[])
{
    for(rs_s32 i = 0; i < numErrors; ++i) {
        rs_u8 l = RS_GF_NW1 - gflog[pos[i]];
        rs_u8 d = gf_sigma_dash_value(numSigma, sigma, l);
        rs_u8 o = gf_omega_value(numOmega, omega, l);
        rs_s32 p = length - 1 - gflog[pos[i]];
        result[p] ^= gf_divexp(gf_div(o, d), l);
    }
}

rs_s32 rs_decode(RSContext* context, rs_s32 size, rs_u8 message[], rs_s32 numSymbols)
{
    rs_s32 messageSize = size + numSymbols;
    CPPRS_ASSERT(CPPRS_STATIC_CAST(rs_u32)(messageSize) < RS_GF_NW);
    CPPRS_ASSERT(numSymbols <= RS_MAX_ECC_SIZE);

    rs_u8* syndromes = context->syndromes_;
    rs_s32 hasError = 0;
    for(rs_s32 i = 0; i < numSymbols; ++i) {
        syndromes[i] = gf_poly_eval(messageSize, message, gfexp[i]);
        hasError |= syndromes[i];
    }

    if(0 == hasError) {
        return 0;
    }

    rs_u8* sigma = context->sigma_;
    rs_s32 numSigma = rs_modified_berlekamp_massey(context, sigma, numSymbols, syndromes);
    if(numSigma < 0) {
        return RS_ERROR;
    }

    rs_u8* errorPositions = context->errorPositions_;
    rs_s32 numErrorPositions = rs_chien_search(errorPositions, CPPRS_STATIC_CAST(rs_u8)(messageSize), CPPRS_STATIC_CAST(rs_u8)(numSigma), sigma);
    if(numErrorPositions < 0) {
        return RS_ERROR;
    }
    rs_u8* omega = context->omega_;
    rs_s32 numOmega = gf_poly_mul(omega, numSymbols, syndromes, numSigma, sigma, numSigma - 1);

    rs_error_correct_forney(message, messageSize, numErrorPositions, errorPositions, numSigma, sigma, numOmega, omega);
    return numSigma - 1;
}

CPPRS_NAMESPACE_EMPTY_BEGIN
CPPRS_NAMESPACE_EMPTY_END

CPPRS_NAMESPACE_END(cpprs)
#endif // CPPRS_IMPLEMENTATION
