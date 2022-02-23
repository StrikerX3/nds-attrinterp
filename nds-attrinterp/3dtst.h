#pragma once

#include <cerrno>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

struct TestData {
    static constexpr uint32_t kMagic = 0x43535344; // DSSC
    static constexpr uint16_t kMaxVersion = 2;

    template <typename T>
    static T Read(std::istream &is) {
        T value{};
        is.read(reinterpret_cast<char *>(&value), sizeof(value));
        return value;
    };

    bool Read(std::filesystem::path path) {
        if (!std::filesystem::exists(path)) {
            return false;
        }
        if (!std::filesystem::is_regular_file(path)) {
            return false;
        }

        if (std::filesystem::file_size(path) < 256 * 192 * 2 + 0x40) {
            return false;
        }

        std::basic_ifstream<char> file{path, std::ios::binary};
        auto magic = Read<uint32_t>(file);
        if (magic != kMagic) {
            return false;
        }

        auto version = Read<uint16_t>(file);
        if (version > kMaxVersion) {
            return false;
        }

        // Only accept files that use the 128x128 coordinate grid texture
        file.seekg(0x8, std::ios::beg);
        auto params = Read<uint32_t>(file);
        if (version != 2 || ((params >> 4) & 7) != 0) {
            return false;
        }

        // Ignore files with translucency unless wireframe is enabled
        if (((params >> 2) & 1) == 0) {
            uint8_t alpha = ((params >> 8) & 0x1F);
            if (alpha > 0 && alpha < 31) {
                return false;
            }
        }

        valid = true;
        frameLoaded = false;
        this->path = path;
        this->magic = magic;
        this->version = version;
        antiAlias = (params >> 0) & 1;
        edgeMark = (params >> 1) & 1;
        wireframe = (params >> 2) & 1;
        quad = (params >> 3) & 1;
        texMode = (params >> 4) & 0xF;
        alpha = (params >> 8) & 0x1F;

        file.seekg(0x10, std::ios::beg);
        for (size_t i = 0; i < 4; i++) {
            for (size_t j = 0; j < 3; j++) {
                verts[i][j] = Read<int32_t>(file);
            }
        }

        return true;
    }

    void LoadFrame() {
        if (!frameLoaded) {
            std::basic_ifstream<char> file{path, std::ios::binary};
            file.seekg(0x40, std::ios::beg);
            file.read(reinterpret_cast<char *>(frame), sizeof(frame));
            frameLoaded = true;
        }
    }

    bool valid = false;
    bool frameLoaded = false;
    std::filesystem::path path;

    uint32_t magic;
    uint16_t version;

    bool antiAlias;
    bool edgeMark;
    bool wireframe;
    bool quad;
    uint8_t texMode;
    uint8_t alpha;

    int32_t verts[4][3];
    uint16_t frame[192][256];
};
