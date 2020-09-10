#define CPPRS_IMPLEMENTATION
#include "cpprs.h"

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <random>
#include <chrono>

using namespace cpprs;

struct Result
{
    long long encodeTime_;
    long long decodeTime_;
    rs_s32 numErrors_;
};

void print(rs_s32 size, std::vector<rs_u8>& v)
{
    for(rs_s32 i=0; i<size; ++i){
        std::cout << int(v[i]) << ", ";
    }
}

#if 0
void error_check()
{
    static const rs_s32 MessageSize = 16;
    static const rs_s32 ECCSize = 10;
    rs_u8 message[MessageSize] = {110, 211, 97, 221, 35, 153, 52, 124, 191, 109, 194, 65, 59, 242, 74, 22};
    rs_u8 diff[MessageSize+ECCSize] = {0, 0, 0, 92, 0, 237, 0, 0, 0, 8, 153, 0, 0, 0, 0, 0, 0, 0, 0, 0, 161, 0, 0, 0, 0, 0};
    rs_u8 encoded[MessageSize + ECCSize];
    rs_u8 decoded[MessageSize + ECCSize];

    RSContext context;
    gf_initialize(&context, ECCSize);

    std::copy(message, message+MessageSize, encoded);
    rs_encode(&context, MessageSize, &encoded[0], ECCSize);
    for(rs_s32 j = 0; j < (MessageSize + ECCSize); ++j) {
        decoded[j] = encoded[j] ^ diff[j];
    }

    rs_s32 corrected = rs_decode(&context, MessageSize, &decoded[0], ECCSize);
    for(rs_s32 j = 0; j < MessageSize; ++j) {
        if(message[j] != decoded[j]) {
            assert(false);
        }
    }
}
#endif

Result reed_solomon(rs_s32 messageSize, rs_s32 eccSize, rs_s32 maxErrors, rs_s32 count)
{
    const rs_s32 MaxECC = eccSize>>1;
    std::vector<rs_u8> message;
    message.resize(messageSize);
    std::vector<rs_u8> encoded;
    encoded.resize(messageSize+eccSize);
    std::vector<rs_u8> decoded;
    decoded.resize(messageSize+eccSize);

    RSContext context;
    Result result = {};

    cpprs::gf_initialize(&context, eccSize);

    std::random_device seed;
    std::mt19937 engine(seed());
    std::uniform_int_distribution<> randErrors(0, maxErrors);
    std::uniform_int_distribution<> randPositions(0, messageSize+eccSize-1);

    long long encodeTime = 0;
    long long decodeTime = 0;
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    for(rs_s32 i = 0; i < count; ++i) {
        for(rs_s32 j=0; j<messageSize; ++j){
            message[j] = static_cast<rs_u8>(engine()&0xFFU);
        }
        std::copy(message.begin(), message.end(), encoded.begin());

        start = std::chrono::high_resolution_clock::now();
        cpprs::rs_encode(&context, messageSize, &encoded[0], eccSize);
        end = std::chrono::high_resolution_clock::now();
        encodeTime += std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();

        std::copy(encoded.begin(), encoded.end(), decoded.begin());

        cpprs::rs_s32 numErrors = randErrors(engine);
        for(cpprs::rs_s32 j=0; j<numErrors; ++j){
            decoded[randPositions(engine)] ^= (engine());
        }
        numErrors = 0;
        for(rs_s32 j=0; j<(messageSize+eccSize); ++j){
            if(encoded[j] != decoded[j]){
                ++numErrors;
            }
            encoded[j] ^= decoded[j];
        }

        start = std::chrono::high_resolution_clock::now();
        cpprs::rs_s32 corrected = cpprs::rs_decode(&context, messageSize, &decoded[0], eccSize);
        end = std::chrono::high_resolution_clock::now();
        decodeTime += std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();

        if(numErrors<=MaxECC){
            if(corrected<0){
                result.numErrors_ += 1;
            }else{
                for(rs_s32 j=0; j<messageSize; ++j){
                    if(message[j] != decoded[j]){
                        std::cout << "message size: " << messageSize << ", ecc size: " << eccSize << std::endl;
                        std::cout << "message: ";
                        print(messageSize, message);
                        std::cout << std::endl;
                        std::cout << "diff: ";
                        print(messageSize+eccSize, encoded);
                        std::cout << std::endl;
                        assert(false);
                    }
                }
            }
        }else{
            if(0<=corrected){
                result.numErrors_ += 1;
                //for(rs_s32 j = 0; j < messageSize; ++j) {
                //    if(message[j] != decoded[j]) {
                //        assert(false);
                //    }
                //}
            }
        }
    }

    result.encodeTime_ = encodeTime / count;
    result.decodeTime_ = decodeTime / count;
    return result;
}

int main(void)
{
    static int Patterns = 10;
    static const int Count = 4096;
    int MinMessageSize = 16;

    std::random_device seed;
    std::mt19937 engine(seed());
    std::uniform_int_distribution<> messageSizeRange(MinMessageSize, RS_GF_NW1-RS_MAX_ECC_SIZE);

    Result result;
    for(rs_s32 i = 0; i < Patterns; ++i) {
        rs_s32 messageSize = messageSizeRange(engine);
        rs_s32 eccMaxSize = std::min(RS_MAX_ECC_SIZE, messageSize-1);
        std::uniform_int_distribution<> eccSizeRange(1, eccMaxSize);
        rs_s32 eccSize = eccSizeRange(engine);
        result = reed_solomon(messageSize, eccSize, eccSize, Count);
        std::cout << "result: " << "message size: " << messageSize << " ecc size: " << eccSize << " encode time: " << result.encodeTime_ << " decode time: " << result.decodeTime_ << " errors: " << result.numErrors_ << std::endl;
    }
    return 0;
}
