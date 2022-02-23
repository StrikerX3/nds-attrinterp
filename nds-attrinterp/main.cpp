#include <fmt/format.h>

#include "3dtst.h"
#include "ansi_console.h"
#include "polygon.h"
#include "rasterizer.h"

#include <array>
#include <cstdint>
#include <vector>

constexpr bool printColors = true;
constexpr bool showMatches = false;
constexpr bool computeSlopes = true;

constexpr int16_t kW0 = 4096;
constexpr int16_t kW1 = 4096;

union Color15 {
    uint16_t u16;
    struct {
        uint16_t r : 5;
        uint16_t g : 5;
        uint16_t b : 5;
        uint16_t a : 1;
    };
};

struct TexCoord {
    int16_t s;
    int16_t t;
};

int BitRev5(int n) {
    return (((n >> 0) & 1) << 4) | (((n >> 1) & 1) << 3) | (((n >> 2) & 1) << 2) | (((n >> 3) & 1) << 1) |
           (((n >> 4) & 1) << 0);
}

TexCoord ToTexCoord(const Color15 &clr) {
    TexCoord texCoord;
    texCoord.s = BitRev5(clr.r) | (((BitRev5(clr.b) >> 1) & 0x3) << 5);
    texCoord.t = BitRev5(clr.g) | (((BitRev5(clr.b) >> 3) & 0x3) << 5);
    return texCoord;
}

Color15 ToColor(const TexCoord &texCoord) {
    Color15 clr;
    clr.r = BitRev5(texCoord.s & 31);
    clr.g = BitRev5(texCoord.t & 31);
    clr.b = BitRev5(((texCoord.s >> 5) | ((texCoord.t >> 5) * 4)) * 2 + 1);
    clr.a = 1;
    return clr;
}

struct EdgeParams {
    size_t cvIndex = ~0, nvIndex = ~0;
    int32_t x0 = -1, y0 = -1;
    int32_t x1 = -1, y1 = -1;
    int32_t sMin = -1, sMax = -1;
    int32_t tMin = -1, tMax = -1;
    int32_t exactS = -1, exactT = -1;
};

struct AttrParams {
    int32_t sMin, sMax;
    int32_t tMin, tMax;
};

struct SlopeParams {
    bool leftEdge;
    int32_t x0, y0;
    int32_t x1, y1;
    int32_t s0, t0;
    int32_t s1, t1;
    int32_t attrFirstY = -1;
    std::vector<AttrParams> attrParams;
};

void test(const std::filesystem::path &path) {
    TestData data;
    if (!data.Read(path)) {
        return;
    }
    data.LoadFrame();
    fmt::print("Testing {}...\n", path.string());

    Polygon polygon{data};
    polygon.w0 = kW0;
    polygon.w1 = kW1;
    for (size_t i = 0; i < polygon.vertexCount; i++) {
        fmt::print("  [{}] {}x{}\n", i, polygon.verts[i][0], polygon.verts[i][1]);
    }
    fmt::print("  {}\n", polygon.windingCW ? "clockwise" : "counter-clockwise");
    fmt::print("  Alpha {}\n", polygon.alpha);

    size_t matchCount = 0;
    size_t totalCount = 0;

    // Collect edge parameters per scanline
    std::array<EdgeParams, 192> leftEdgeParams{};
    std::array<EdgeParams, 192> rightEdgeParams{};

    // Texture coordinates for the 128x128 texture
    static constexpr int32_t kTexCoords[4][2] = {{0, 0}, {0, 2048}, {2048, 2048}, {2048, 0}};

    PolygonState state;
    state.SetupPolygon(polygon);
    size_t ys = std::max(0, polygon.topScreenY);
    size_t ye = std::min(191, polygon.btmScreenY);
    for (size_t y = ys; y <= ye; y++) {
        if (y == polygon.btmScreenY && state.tall) {
            continue;
        }

        // Update slopes if the polygon is taller than a single scanline
        if (state.tall) {
            state.CalcSlopes(y);
        }

        auto &leftEdge = state.leftEdge;
        auto &rightEdge = state.rightEdge;

        // Get vertex indices and edge slopes
        size_t cvlIndex = leftEdge.currVertex;
        size_t nvlIndex = leftEdge.nextVertex;
        size_t cvrIndex = rightEdge.currVertex;
        size_t nvrIndex = rightEdge.nextVertex;
        Slope *leftSlope = &leftEdge.slope;
        Slope *rightSlope = &rightEdge.slope;

        // Determine edge spans for this scanline
        int32_t fxls = leftEdge.slope.FracXStart(y);
        int32_t fxle = (leftEdge.slope.DX() != 0) ? leftEdge.slope.FracXEnd(y) : fxls;
        int32_t fxrs = rightEdge.slope.FracXStart(y);
        int32_t fxre = (rightEdge.slope.DX() != 0) ? rightEdge.slope.FracXEnd(y) : fxrs;
        if (leftEdge.slope.IsNegative()) {
            std::swap(fxls, fxle);
        }
        if (rightEdge.slope.IsNegative()) {
            std::swap(fxrs, fxre);
        }

        // Remove fractional parts
        int32_t xls = fxls >> Slope::kFracBits;
        int32_t xle = fxle >> Slope::kFracBits;
        int32_t xrs = fxrs >> Slope::kFracBits;
        int32_t xre = fxre >> Slope::kFracBits;

        // Swap edges if they are reversed
        if (xls > xrs || xle > xre) {
            std::swap(fxls, fxrs);
            std::swap(fxle, fxre);
            std::swap(xls, xrs);
            std::swap(xle, xre);
            std::swap(cvlIndex, cvrIndex);
            std::swap(nvlIndex, nvrIndex);
            std::swap(leftSlope, rightSlope);

            // Slope spans break when edges are reversed
            xls = xle;
            xre = xrs;
            fxls = fxle;
            fxre = fxrs;
            if (xls > xrs) {
                std::swap(xls, xrs);
                std::swap(xle, xre);
                std::swap(fxls, fxrs);
                std::swap(fxle, fxre);
            }
        }

        // Update interpolators
        leftEdge.slope.ComputeFactors(y);
        rightEdge.slope.ComputeFactors(y);

        // Interpolate texture coordinates for the left edge
        const int32_t sl = leftSlope->InterpolateAttribute(kTexCoords[cvlIndex][0], kTexCoords[nvlIndex][0]);
        const int32_t tl = leftSlope->InterpolateAttribute(kTexCoords[cvlIndex][1], kTexCoords[nvlIndex][1]);

        // Interpolate texture coordinates for the right edge
        const int32_t sr = rightSlope->InterpolateAttribute(kTexCoords[cvrIndex][0], kTexCoords[nvrIndex][0]);
        const int32_t tr = rightSlope->InterpolateAttribute(kTexCoords[cvrIndex][1], kTexCoords[nvrIndex][1]);

        // Setup X coordinate interpolator
        Interpolator xInterp;
        xInterp.Setup(xls, xre + (rightSlope->DX() != 0), kW0, kW1);

        // Determine which portions of the polygon to render
        bool drawLeftEdge;
        bool drawInterior = (polygon.alpha > 0) || (y == polygon.topScreenY) || (y == polygon.btmScreenY - 1);
        bool drawRightEdge;
        if (polygon.alpha < 31 || data.antiAlias || data.edgeMark || y == 191) {
            drawLeftEdge = true;
            drawRightEdge = true;
        } else {
            drawLeftEdge = leftSlope->IsNegative() || !leftSlope->IsXMajor();
            drawRightEdge = (!rightSlope->IsNegative() && rightSlope->IsXMajor()) || rightSlope->DX() == 0;
        }

        // Nudge right edge to the left by one pixel when it's perfectly vertical
        if (rightSlope->DX() == 0) {
            xrs--;
            xre--;
        }

        // Trim right edge in some rare cases where it starts before the left edge
        xrs = std::max(xrs, xle + 1);

        bool hadMismatch = false;

        auto cvt5to8 = [](uint8_t x) { return (x << 3) | (x >> 2); };

        auto testSpan = [&](int32_t xs, int32_t xe) {
            for (int32_t x = xs; x <= xe; x++) {
                // Update X interpolation factors
                xInterp.ComputeFactors(x);

                // Interpolate texture coordinates along X axis
                const int16_t fs = std::clamp(xInterp.InterpolateAttribute(sl, sr), 0, 2047);
                const int16_t ft = std::clamp(xInterp.InterpolateAttribute(tl, tr), 0, 2047);

                // Remove fractional part
                const int16_t s = fs >> 4;
                const int16_t t = ft >> 4;

                // Convert to color to aid comparisons
                const Color15 stClr = ToColor({s, t});

                // Apply alpha
                const Color15 clr = [&]() -> Color15 {
                    if (polygon.alpha == 0 || polygon.alpha == 31) {
                        return stClr;
                    }
                    uint8_t alpha = data.alpha + 1;
                    auto apply = [&](uint16_t c) -> uint16_t { return ((c * alpha) + (3 * (32 - alpha))) >> 5; };
                    return {.r = apply(stClr.r), .g = apply(stClr.g), .b = apply(stClr.b), .a = 1};
                }();

                // Compare against captured data
                totalCount++;
                const Color15 frameClr = {.u16 = data.frame[y][x]};
                const auto texCoord = ToTexCoord(frameClr);
                if (frameClr.r == 3 && frameClr.g == 3 && frameClr.b == 3) {
                    // Pixel present where there should be background.
                    // Very few cases occur, and they are never an attribute interpolation error.
                    fmt::print(
                        "  /!\\ {:>3d}x{:<3d} | {:>4d}x{:<4d} ({:>4s}x{:>4s}) >> {:>3d}x{:<3d} != background         ",
                        x, y, fs, ft, fmt::format("{:X}.{:X}", s, fs & 0xF), fmt::format("{:X}.{:X}", t, ft & 0xF), s,
                        t);
                    fmt::print("   [{:2d},{:2d},{:2d}] ", clr.r, clr.g, clr.b);
                    if constexpr (printColors) {
                        ansicon::Attributes attrs;
                        attrs.SetBGColor24(cvt5to8(clr.r), cvt5to8(clr.g), cvt5to8(clr.b));
                        fmt::print("  ");
                    }
                    fmt::print(" != ");
                    if constexpr (printColors) {
                        ansicon::Attributes attrs;
                        attrs.SetBGColor24(cvt5to8(3), cvt5to8(3), cvt5to8(3));
                        fmt::print("  ");
                    }
                    fmt::print(" [ 3, 3, 3]");
                    fmt::print("   {} / {}\n", xInterp.Num(), xInterp.Den());
                    hadMismatch = true;
                } else if (texCoord.s != s || texCoord.t != t) {
                    // Mismatched pixel.
                    fmt::print("  /!\\ {:>3d}x{:<3d} | {:>4d}x{:<4d} ({:>4s}x{:>4s}) >> {:>3d}x{:<3d} != {:>3d}x{:<3d} "
                               "{:<11s}",
                               x, y, fs, ft, fmt::format("{:X}.{:X}", s, fs & 0xF),
                               fmt::format("{:X}.{:X}", t, ft & 0xF), s, t, texCoord.s, texCoord.t,
                               fmt::format("({:d}x{:d})", s - texCoord.s, t - texCoord.t));
                    fmt::print("   [{:2d},{:2d},{:2d}] ", clr.r, clr.g, clr.b);
                    if constexpr (printColors) {
                        ansicon::Attributes attrs;
                        attrs.SetBGColor24(cvt5to8(clr.r), cvt5to8(clr.g), cvt5to8(clr.b));
                        fmt::print("  ");
                    }
                    fmt::print(" != ");
                    if constexpr (printColors) {
                        ansicon::Attributes attrs;
                        attrs.SetBGColor24(cvt5to8(frameClr.r), cvt5to8(frameClr.g), cvt5to8(frameClr.b));
                        fmt::print("  ");
                    }
                    fmt::print(" [{:2d},{:2d},{:2d}]", frameClr.r, frameClr.g, frameClr.b);
                    fmt::print("   {} / {}\n", xInterp.Num(), xInterp.Den());
                    hadMismatch = true;
                } else {
                    matchCount++;
                    if constexpr (showMatches) {
                        fmt::print(
                            "      {:>3d}x{:<3d} | {:>4d}x{:<4d} ({:>4s}x{:>4s}) >> {:>3d}x{:<3d} == {:>3d}x{:<3d}  "
                            "          ",
                            x, y, fs, ft, fmt::format("{:X}.{:X}", s, fs & 0xF), fmt::format("{:X}.{:X}", t, ft & 0xF),
                            s, t, texCoord.s, texCoord.t, s - texCoord.s, t - texCoord.t);
                        fmt::print("   [{:2d},{:2d},{:2d}] ", clr.r, clr.g, clr.b);
                        if constexpr (printColors) {
                            ansicon::Attributes attrs;
                            attrs.SetBGColor24(cvt5to8(clr.r), cvt5to8(clr.g), cvt5to8(clr.b));
                            fmt::print("  ");
                        }
                        fmt::print(" == ");
                        if constexpr (printColors) {
                            ansicon::Attributes attrs;
                            attrs.SetBGColor24(cvt5to8(frameClr.r), cvt5to8(frameClr.g), cvt5to8(frameClr.b));
                            fmt::print("  ");
                        }
                        fmt::print(" [{:2d},{:2d},{:2d}]", frameClr.r, frameClr.g, frameClr.b);
                        fmt::print("   {} / {}\n", xInterp.Num(), xInterp.Den());
                    }
                }
            };
        };

        bool anyPixelDrawn = false;

        // Test left edge
        if (drawLeftEdge) {
            const int32_t xStart = std::max(xls, 0);
            const int32_t xEnd = std::min(xle, (int32_t)255);
            anyPixelDrawn = (xStart <= xEnd);
            testSpan(xStart, xEnd);
        }

        // Test polygon interior
        if (drawInterior) {
            const int32_t xStart = std::max(xle + 1, 0);
            const int32_t xEnd = std::min(xrs, (int32_t)256);
            anyPixelDrawn |= (xStart < xEnd);
            testSpan(xStart, xEnd - 1);
        }

        // Test right edge
        if (drawRightEdge) {
            const int32_t xStart = std::max(xrs, 0);
            const int32_t xEnd = std::min(xre, (int32_t)255);
            anyPixelDrawn |= (xStart <= xEnd);
            testSpan(xStart, xEnd);
        }

        if (showMatches || hadMismatch) {
            auto &cvl = data.verts[cvlIndex];
            auto &nvl = data.verts[nvlIndex];
            auto &ls = *leftSlope;
            auto &li = ls.Interp();

            auto &cvr = data.verts[cvrIndex];
            auto &nvr = data.verts[nvrIndex];
            auto &rs = *rightSlope;
            auto &ri = rs.Interp();

            fmt::print(
                "    Left:  edge={:>3d}x{:<3d}..{:>3d}x{:<3d}  texcoords={:>4d}x{:<4d} = {:>3d}x{:<3d}   slope={}{} "
                "{:>3d}..{:<3d} ({:>6d}..{:<6d}) DX={:<10d}   interp={:>3d}..{:<3d} ({})  {}\n",
                cvl[0], cvl[1], nvl[0], nvl[1], sl, tl, sl >> 4, tl >> 4, (ls.IsNegative() ? '-' : '+'),
                (ls.IsXMajor() ? 'X' : 'Y'), ls.XStart(y), ls.XEnd(y), ls.FracXStart(y) & Slope::kFracMask,
                ls.FracXEnd(y) & Slope::kFracMask, ls.DX(), li.X0(), li.X1(), li.XMax(),
                (drawLeftEdge ? "drawn" : "hidden"));

            fmt::print(
                "    Right: edge={:>3d}x{:<3d}..{:>3d}x{:<3d}  texcoords={:>4d}x{:<4d} = {:>3d}x{:<3d}   slope={}{} "
                "{:>3d}..{:<3d} ({:>6d}..{:<6d}) DX={:<10d}   interp={:>3d}..{:<3d} ({})  {}\n",
                cvr[0], cvr[1], nvr[0], nvr[1], sr, tr, sr >> 4, tr >> 4, (rs.IsNegative() ? '-' : '+'),
                (rs.IsXMajor() ? 'X' : 'Y'), rs.XStart(y), rs.XEnd(y), rs.FracXStart(y) & Slope::kFracMask,
                rs.FracXEnd(y) & Slope::kFracMask, rs.DX(), ri.X0(), ri.X1(), ri.XMax(),
                (drawRightEdge ? "drawn" : "hidden"));

            fmt::print("    X interpolator: {:>3d}..{:<3d} ({})   interior {}\n", xInterp.X0(), xInterp.X1(),
                       xInterp.XMax(), (drawInterior ? "drawn" : "hidden"));
        }

        if constexpr (computeSlopes) {
            // Find properly matching parameters regardless of whether we found a match or not

            // TODO: optimize; this is a dumb brute-force algorithm
            // - could use the existing values as a starting point since they are very close to the correct values
            // - nudge up or down depending on whether the interpolation undershoots or overshoots the captured data
            // - only problem is when both left and right edge attributes are not constant
            //   - which of the sides to nudge?
            //     - maybe check the delta gradient
            //       - split dataset in exactly half (add middle entry to both if odd), compare halves
            //     - if overshooting and gradient increases (more +) towards right edge, decrement right edge
            //     - if overshooting and gradient decreases (more +) towards right edge, decrement left edge
            //     - if undershooting and gradient increases (more -) towards right edge, increment right edge
            //     - if undershooting and gradient decreases (more -) towards right edge, increment left edge

            // Get attribute range limits
            int32_t slLo = kTexCoords[cvlIndex][0];
            int32_t slHi = kTexCoords[nvlIndex][0];
            if (slLo > slHi) {
                std::swap(slLo, slHi);
            }

            int32_t srLo = kTexCoords[cvrIndex][0];
            int32_t srHi = kTexCoords[nvrIndex][0];
            if (srLo > srHi) {
                std::swap(srLo, srHi);
            }

            int32_t tlLo = kTexCoords[cvlIndex][1];
            int32_t tlHi = kTexCoords[nvlIndex][1];
            if (tlLo > tlHi) {
                std::swap(tlLo, tlHi);
            }

            int32_t trLo = kTexCoords[cvrIndex][1];
            int32_t trHi = kTexCoords[nvrIndex][1];
            if (trLo > trHi) {
                std::swap(trLo, trHi);
            }

            const int32_t start = xls;
            const int32_t end = xre + (rightSlope->DX() != 0);

            auto testAttr = [&](int32_t l, int32_t r, bool s) {
                Interpolator xInterp;
                xInterp.Setup(start, end, kW0, kW1);
                auto test = [&](int32_t xs, int32_t xe) {
                    for (int32_t x = xs; x <= xe; x++) {
                        // Update X interpolation factors
                        xInterp.ComputeFactors(x);

                        // Interpolate texture coordinate along X axis
                        const int16_t frac = std::clamp(xInterp.InterpolateAttribute(l, r), 0, 2047);

                        // Remove fractional part
                        const int16_t integ = frac >> 4;

                        // Compare against captured data
                        const Color15 frameClr = {.u16 = data.frame[y][x]};
                        const auto texCoord = ToTexCoord(frameClr);
                        const int32_t targetCoord = (s ? texCoord.s : texCoord.t);
                        if (frameClr.r == 3 && frameClr.g == 3 && frameClr.b == 3) {
                            // Ignore background pixels; not a problem with attribute interpolation
                        } else if (targetCoord != integ) {
                            // Found a mismatch
                            return false;
                        }
                    }
                    return true;
                };

                // Test left edge
                if (drawLeftEdge) {
                    const int32_t xStart = std::max(xls, 0);
                    const int32_t xEnd = std::min(xle, (int32_t)255);
                    if (!test(xStart, xEnd)) {
                        return false;
                    }
                }

                // Test polygon interior
                if (drawInterior) {
                    const int32_t xStart = std::max(xle + 1, 0);
                    const int32_t xEnd = std::min(xrs, (int32_t)256);
                    if (!test(xStart, xEnd - 1)) {
                        return false;
                    }
                }

                // Test right edge
                if (drawRightEdge) {
                    const int32_t xStart = std::max(xrs, 0);
                    const int32_t xEnd = std::min(xre, (int32_t)255);
                    if (!test(xStart, xEnd)) {
                        return false;
                    }
                }
                return true;
            };

            int32_t slMin = -1;
            int32_t slMax = -1;
            int32_t srMin = -1;
            int32_t srMax = -1;

            int32_t tlMin = -1;
            int32_t tlMax = -1;
            int32_t trMin = -1;
            int32_t trMax = -1;

            bool foundMin = false;
            bool foundMax = false;
            for (int32_t slTest = slLo; slTest <= slHi; slTest++) {
                for (int32_t srTest = srLo; srTest <= srHi; srTest++) {
                    if (testAttr(slTest, srTest, true)) {
                        // Found a match for S
                        if (foundMin) {
                            slMax = slTest;
                            srMax = srTest;
                        } else {
                            slMin = slMax = slTest;
                            srMin = srMax = srTest;
                            foundMin = true;
                        }
                    } else if (foundMin) {
                        foundMax = true;
                        break;
                    }
                }
                if (foundMax) {
                    break;
                }
            }

            foundMin = false;
            foundMax = false;
            for (int32_t tlTest = tlLo; tlTest <= tlHi; tlTest++) {
                for (int32_t trTest = trLo; trTest <= trHi; trTest++) {
                    if (testAttr(tlTest, trTest, false)) {
                        // Found a match for T
                        if (foundMin) {
                            tlMax = tlTest;
                            trMax = trTest;
                        } else {
                            tlMin = tlMax = tlTest;
                            trMin = trMax = trTest;
                            foundMin = true;
                        }
                    } else if (foundMin) {
                        foundMax = true;
                        break;
                    }
                }
                if (foundMax) {
                    break;
                }
            }

            auto fmtRes = [](int32_t valMin, int32_t valMax) -> std::string {
                if (valMin == -1) {
                    return "??";
                } else if (valMin == valMax) {
                    return fmt::format("{:d}", valMin);
                } else {
                    return fmt::format("{:d}..{:d}", valMin, valMax);
                }
            };

            if (hadMismatch) {
                fmt::print("       Corrected left texcoords:  {:>10s}x{:<10s}\n", fmtRes(slMin, slMax),
                           fmtRes(tlMin, tlMax));
                fmt::print("       Corrected right texcoords: {:>10s}x{:<10s}\n", fmtRes(srMin, srMax),
                           fmtRes(trMin, trMax));
            }

            // Store edge parameters, but only if at least one pixel was drawn on this scanline
            if (anyPixelDrawn) {
                leftEdgeParams[y].cvIndex = cvlIndex;
                leftEdgeParams[y].nvIndex = nvlIndex;
                leftEdgeParams[y].x0 = polygon.verts[cvlIndex][0];
                leftEdgeParams[y].y0 = polygon.verts[cvlIndex][1];
                leftEdgeParams[y].x1 = polygon.verts[nvlIndex][0];
                leftEdgeParams[y].y1 = polygon.verts[nvlIndex][1];

                rightEdgeParams[y].cvIndex = cvrIndex;
                rightEdgeParams[y].nvIndex = nvrIndex;
                rightEdgeParams[y].x0 = polygon.verts[cvrIndex][0];
                rightEdgeParams[y].y0 = polygon.verts[cvrIndex][1];
                rightEdgeParams[y].x1 = polygon.verts[nvrIndex][0];
                rightEdgeParams[y].y1 = polygon.verts[nvrIndex][1];

                leftEdgeParams[y].sMin = slMin;
                leftEdgeParams[y].sMax = slMax;
                leftEdgeParams[y].tMin = tlMin;
                leftEdgeParams[y].tMax = tlMax;

                rightEdgeParams[y].sMin = srMin;
                rightEdgeParams[y].sMax = srMax;
                rightEdgeParams[y].tMin = trMin;
                rightEdgeParams[y].tMax = trMax;

                if (!hadMismatch) {
                    leftEdgeParams[y].exactS = sl;
                    leftEdgeParams[y].exactT = tl;

                    rightEdgeParams[y].exactS = sr;
                    rightEdgeParams[y].exactT = tr;
                }
            }
        }
    }

    fmt::print("Accuracy: {} / {} ({:3.2f}%)\n", matchCount, totalCount, (double)matchCount / totalCount * 100.0);

    if (computeSlopes && matchCount != totalCount) {
        // Aggregate into slopes
        std::vector<SlopeParams> slopes;
        auto processEdges = [&](bool left, std::array<EdgeParams, 192> &edgeParams) {
            for (size_t y = ys; y <= ye; y++) {
                const auto &row = edgeParams[y];
                if (row.x0 != -1) {
                    auto newEdge = [&]() -> SlopeParams & {
                        auto &edge = slopes.emplace_back();
                        edge.x0 = row.x0;
                        edge.y0 = row.y0;
                        edge.x1 = row.x1;
                        edge.y1 = row.y1;
                        edge.s0 = kTexCoords[row.cvIndex][0];
                        edge.t0 = kTexCoords[row.cvIndex][1];
                        edge.s1 = kTexCoords[row.nvIndex][0];
                        edge.t1 = kTexCoords[row.nvIndex][1];
                        edge.leftEdge = left;
                        return edge;
                    };
                    auto matchingEdge = [&]() -> SlopeParams & {
                        if (slopes.empty()) {
                            return newEdge();
                        } else {
                            auto &edge = slopes.back();
                            if (edge.x0 != row.x0 || edge.y0 != row.y0 || edge.x1 != row.x1 || edge.y1 != row.y1) {
                                return newEdge();
                            } else {
                                return edge;
                            }
                        }
                    };
                    auto &edge = matchingEdge();
                    if (edge.attrFirstY == -1) {
                        edge.attrFirstY = y;
                    }
                    edge.attrParams.emplace_back(
                        AttrParams{.sMin = row.sMin, .sMax = row.sMax, .tMin = row.tMin, .tMax = row.tMax});
                }
            }
        };
        processEdges(true, leftEdgeParams);
        processEdges(false, rightEdgeParams);

        fmt::print("Slopes:\n");
        for (auto &slope : slopes) {
            const char side = (slope.leftEdge ? 'L' : 'R');
            const bool negative = (slope.x0 > slope.x1);
            const bool xmajor = (std::abs(slope.x1 - slope.x0) > std::abs(slope.y1 - slope.y0));
            const bool diagonal = (std::abs(slope.x1 - slope.x0) == std::abs(slope.y1 - slope.y0));
            Interpolator interp;
            interp.Setup(slope.y0, slope.y1, kW0, kW1);

            auto fmtRange = [](int32_t min, int32_t max) {
                if (min == max) {
                    return fmt::format("{}", min);
                } else {
                    return fmt::format("{}..{}", min, max);
                }
            };
            fmt::print("  {}{}{} {:>3d}x{:<3d}..{:>3d}x{:<3d}  dx={:<4d}  dy={:4d}\n", side, (negative ? '-' : '+'),
                       (xmajor ? 'X' : 'Y'), slope.x0, slope.y0, slope.x1, slope.y1, slope.x1 - slope.x0,
                       slope.y1 - slope.y0);
            fmt::print("    ranges   s={:<10s}  t={:<10s}\n", fmtRange(slope.s0, slope.s1),
                       fmtRange(slope.t0, slope.t1));
            for (int32_t i = 0; i < slope.attrParams.size(); i++) {
                int32_t y = i + slope.attrFirstY;

                const auto &attrs = slope.attrParams[i];
                fmt::print("    y={:<3d} -> s={:<10s}  t={:<10s}", y, fmtRange(attrs.sMin, attrs.sMax),
                           fmtRange(attrs.tMin, attrs.tMax));

                // L-X and R+X need +1
                // L-Y and R+Y need +1, but only perfect diagonals
                interp.ComputeFactors(y + ((xmajor || diagonal) && (slope.leftEdge == negative)));
                int32_t s = interp.InterpolateAttribute(slope.s0, slope.s1);
                int32_t t = interp.InterpolateAttribute(slope.t0, slope.t1);
                fmt::print("  interp s={:<4d}  t={:<4d}\n", s, t);
            }
        }
    }
}

bool read(const std::filesystem::path &path) {
    if (std::filesystem::is_directory(path)) {
        fmt::print("Reading directory {}...\n", path.string());
        for (const std::filesystem::directory_entry &dir_entry : std::filesystem::recursive_directory_iterator{path}) {
            if (dir_entry.is_regular_file()) {
                test(dir_entry.path());
            }
        }
    } else if (std::filesystem::is_regular_file(path)) {
        test(path);
    } else {
        fmt::print("{} does not exist or cannot be accessed\n", path.string());
        return false;
    }
    return true;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fmt::print("missing argument: path\n");
        return EXIT_FAILURE;
    }

    std::filesystem::path path = argv[1];
    if (read(path)) {
        return EXIT_SUCCESS;
    } else {
        return EXIT_FAILURE;
    }
}
