#pragma once

#include <cstdint>
#include <cstdio>

namespace ansicon {

struct Attributes {
    Attributes()
        : file(stdout) {}

    Attributes(FILE *file)
        : file(file) {}

    ~Attributes() {
        Reset();
    }

    Attributes &Reset() {
        fprintf(file, "\x1B[0m");
        return *this;
    }

    Attributes &SetFGColor24(uint8_t r, uint8_t g, uint8_t b) {
        fprintf(file, "\x1B[38;2;%u;%u;%um", r, g, b);
        return *this;
    }

    Attributes &SetBGColor24(uint8_t r, uint8_t g, uint8_t b) {
        fprintf(file, "\x1B[48;2;%u;%u;%um", r, g, b);
        return *this;
    }

private:
    FILE *file;
};

} // namespace ansicon
