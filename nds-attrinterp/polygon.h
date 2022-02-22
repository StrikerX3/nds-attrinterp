#pragma once

#include "3d_slope.h"
#include "3dtst.h"

#include <array>

struct Polygon {
    std::array<std::array<int32_t, 2>, 4> verts;

    bool windingCW;

    size_t topVtxIndex;
    size_t btmVtxIndex;
    size_t vertexCount;

    int32_t topScreenY;
    int32_t btmScreenY;

    uint8_t alpha;

    Polygon(const TestData &data) {
        vertexCount = data.quad ? 4 : 3;
        alpha = data.wireframe ? 0 : data.alpha;
        topScreenY = btmScreenY = data.verts[0][1];
        topVtxIndex = btmVtxIndex = 0;
        for (size_t i = 0; i < vertexCount; i++) {
            verts[i][0] = data.verts[i][0];
            verts[i][1] = data.verts[i][1];
            if (verts[i][1] < topScreenY) {
                topScreenY = verts[i][1];
                topVtxIndex = i;
            }
            if (verts[i][1] > btmScreenY) {
                btmScreenY = verts[i][1];
                btmVtxIndex = i;
            }
        }

        int32_t x10 = verts[0][0] - verts[1][0];
        int32_t y10 = verts[0][1] - verts[1][1];

        int32_t x12 = verts[2][0] - verts[1][0];
        int32_t y12 = verts[2][1] - verts[1][1];

        int64_t cz = ((int64_t)x12 * y10 - (int64_t)y12 * x10);
        windingCW = (cz > 0);
    }
};
