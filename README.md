# cpprs
Reed-Solomon Codec in C++

# Usage
To create an implementation, put `CPPECC_IMPLEMENTATION`.

``` cpp
#define CPPECC_IMPLEMENTATION
#include "cppecc.h"

#include <cassert>
#include <algorithm>

using namespace cpprs;

int main(void)
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
    return 0;
}
```

# Limitations
I limit the size of redundant symbols (that means the error correction capacity) for my use. It's equivalent to about 10% error correnction capability.
You can change this with the constant `CPPECC_MAX_ECC_SIZE`.

## Generation Polynomial
I choose a primitive `0x11D` for the generation polynomial (from the specification of the QR code).
You may be able to change it, but should regenerate tables 'gfexp' and 'gflog'.

# Warning
I'm not a mathematician, an engineer. Use carefully, when you use this.

# License
This software is distributed under two licenses 'The MIT License' or 'Public Domain', choose whichever you like.

