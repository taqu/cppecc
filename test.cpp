#define CPPECC_IMPLEMENTATION
#include "cppecc.h"

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <random>
#include <chrono>

using namespace cppecc;

struct Result
{
    long long encodeTime_;
    long long decodeTime_;
    cppecc_s32 numErrors_;
    cppecc_s64 avgErrors_;
};

void print(cppecc_s32 size, std::vector<cppecc_u8>& v)
{
    for(cppecc_s32 i=0; i<size; ++i){
        std::cout << int(v[i]) << ", ";
    }
}

#if 0
void error_check()
{
    static const cppecc_s32 MessageSize = 16;
    static const cppecc_s32 ECCSize = 10;
    cppecc_u8 message[MessageSize] = {110, 211, 97, 221, 35, 153, 52, 124, 191, 109, 194, 65, 59, 242, 74, 22};
    cppecc_u8 diff[MessageSize+ECCSize] = {0, 0, 0, 92, 0, 237, 0, 0, 0, 8, 153, 0, 0, 0, 0, 0, 0, 0, 0, 0, 161, 0, 0, 0, 0, 0};
    cppecc_u8 encoded[MessageSize + ECCSize];
    cppecc_u8 decoded[MessageSize + ECCSize];

    RSContext context;
    gf_initialize(&context, ECCSize);

    std::copy(message, message+MessageSize, encoded);
    rs_encode(&context, MessageSize, &encoded[0], ECCSize);
    for(cppecc_s32 j = 0; j < (MessageSize + ECCSize); ++j) {
        decoded[j] = encoded[j] ^ diff[j];
    }

    cppecc_s32 corrected = rs_decode(&context, MessageSize, &decoded[0], ECCSize);
    for(cppecc_s32 j = 0; j < MessageSize; ++j) {
        if(message[j] != decoded[j]) {
            assert(false);
        }
    }
}
#endif

Result reed_solomon(cppecc_s32 messageSize, cppecc_s32 eccSize, cppecc_s32 maxErrors, cppecc_s32 count)
{
    const cppecc_s32 MaxECC = eccSize>>1;
    std::vector<cppecc_u8> message;
    message.resize(messageSize);
    std::vector<cppecc_u8> encoded;
    encoded.resize(messageSize+eccSize);
    std::vector<cppecc_u8> decoded;
    decoded.resize(messageSize+eccSize);

    RSContext context;
    Result result = {};

    gf_initialize(&context, eccSize);

    std::random_device seed;
    std::mt19937 engine(seed());
    std::uniform_int_distribution<> randErrors(0, maxErrors);
    std::uniform_int_distribution<> randPositions(0, messageSize+eccSize-1);

    long long encodeTime = 0;
    long long decodeTime = 0;
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    for(cppecc_s32 i = 0; i < count; ++i) {
        for(cppecc_s32 j=0; j<messageSize; ++j){
            message[j] = static_cast<cppecc_u8>(engine()&0xFFU);
        }
        std::copy(message.begin(), message.end(), encoded.begin());

        start = std::chrono::high_resolution_clock::now();
        rs_encode(&context, messageSize, &encoded[0], eccSize);
        end = std::chrono::high_resolution_clock::now();
        encodeTime += std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();

        std::copy(encoded.begin(), encoded.end(), decoded.begin());

        cppecc_s32 numErrors = randErrors(engine);
        result.avgErrors_ += numErrors;
        for(cppecc_s32 j=0; j<numErrors; ++j){
            decoded[randPositions(engine)] ^= (engine());
        }
        numErrors = 0;
        for(cppecc_s32 j=0; j<(messageSize+eccSize); ++j){
            if(encoded[j] != decoded[j]){
                ++numErrors;
            }
            encoded[j] ^= decoded[j];
        }

        start = std::chrono::high_resolution_clock::now();
        cppecc_s32 corrected = rs_decode(&context, messageSize, &decoded[0], eccSize);
        end = std::chrono::high_resolution_clock::now();
        decodeTime += std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();

        if(numErrors<=MaxECC){
            if(corrected<0){
                result.numErrors_ += 1;
            }else{
                for(cppecc_s32 j=0; j<messageSize; ++j){
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
            }
        }
    }

    result.encodeTime_ = encodeTime / count;
    result.decodeTime_ = decodeTime / count;
    result.avgErrors_ /= count;
    return result;
}

int main(void)
{
    static const int Patterns = 1024;
    static const int Count = 4096;
    int MinMessageSize = 16;

    std::random_device seed;
    std::mt19937 engine(seed());
    std::uniform_int_distribution<> messageSizeRange(MinMessageSize, CPPECC_GF_NW1-CPPECC_MAX_ECC_SIZE);

    Result result;
    for(cppecc_s32 i = 0; i < Patterns; ++i) {
        cppecc_s32 messageSize = messageSizeRange(engine);
        cppecc_s32 eccMaxSize = std::min(CPPECC_MAX_ECC_SIZE, messageSize-1);
        std::uniform_int_distribution<> eccSizeRange(1, eccMaxSize);
        cppecc_s32 eccSize = eccSizeRange(engine);
        result = reed_solomon(messageSize, eccSize, eccSize, Count);
        std::cout << "result: " << "message size: " << messageSize << " ecc size: " << eccSize << " encode time (micro): " << result.encodeTime_ << " decode time (micro): " << result.decodeTime_ << " errors: " << result.numErrors_ << '/' << result.avgErrors_ << std::endl;
    }
    return 0;
}
