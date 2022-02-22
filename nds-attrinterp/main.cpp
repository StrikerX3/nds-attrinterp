#include <fmt/format.h>

#include "3dtst.h"
#include "ansi_console.h"
#include "polygon.h"
#include "rasterizer.h"

constexpr bool printColors = true;
constexpr bool showMatches = true;

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

void test(const std::filesystem::path &path) {
    TestData data;
    if (!data.Read(path)) {
        return;
    }
    fmt::print("Testing {}...\n", path.string());

    Polygon polygon{data};
    for (size_t i = 0; i < polygon.vertexCount; i++) {
        fmt::print("  [{}] {}x{}\n", i, polygon.verts[i][0], polygon.verts[i][1]);
    }
    fmt::print("  {}\n", polygon.windingCW ? "clockwise" : "counter-clockwise");
    fmt::print("  Alpha {}\n", polygon.alpha);

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
        leftEdge.slope.ComputeFactors(y, false);
        rightEdge.slope.ComputeFactors(y, true);

        // Texture coordinates for the 128x128 texture
        static constexpr int32_t kTexCoords[4][2] = {{0, 0}, {0, 2048}, {2048, 2048}, {2048, 0}};

        // Interpolate texture coordinates for the left edge
        const int32_t sl = leftSlope->InterpolateAttribute(kTexCoords[cvlIndex][0], kTexCoords[nvlIndex][0]);
        const int32_t tl = leftSlope->InterpolateAttribute(kTexCoords[cvlIndex][1], kTexCoords[nvlIndex][1]);

        // Interpolate texture coordinates for the right edge
        const int32_t sr = rightSlope->InterpolateAttribute(kTexCoords[cvrIndex][0], kTexCoords[nvrIndex][0]);
        const int32_t tr = rightSlope->InterpolateAttribute(kTexCoords[cvrIndex][1], kTexCoords[nvrIndex][1]);

        // Setup X coordinate interpolator
        Interpolator xInterp;
        xInterp.Setup(xls, xre + (!rightSlope->IsXMajor() && rightSlope->DX() != 0), 4096, 4096);

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

        // bool hadMismatch = false;

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

                auto cvt5to8 = [](uint8_t x) { return (x << 3) | (x >> 2); };

                // Compare against captured data
                const Color15 frameClr = {.u16 = data.frame[y][x]};
                const auto texCoord = ToTexCoord(frameClr);
                if (frameClr.r == 3 && frameClr.g == 3 && frameClr.b == 3) {
                    fmt::print(
                        "  /!\\ {:>3d}x{:<3d} | {:>4d}x{:<4d} ({:>4s}x{:>4s}) >> {:>3d}x{:<3d} != backgnd            ",
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
                    // hadMismatch = true;
                } else if (texCoord.s != s || texCoord.t != t) {
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
                    // hadMismatch = true;
                } else if constexpr (showMatches) {
                    fmt::print("      {:>3d}x{:<3d} | {:>4d}x{:<4d} ({:>4s}x{:>4s}) >> {:>3d}x{:<3d} == {:>3d}x{:<3d}  "
                               "          ",
                               x, y, fs, ft, fmt::format("{:X}.{:X}", s, fs & 0xF),
                               fmt::format("{:X}.{:X}", t, ft & 0xF), s, t, texCoord.s, texCoord.t, s - texCoord.s,
                               t - texCoord.t);
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
            };
        };

        // Test left edge
        if (drawLeftEdge) {
            const int32_t xStart = std::max(xls, 0);
            const int32_t xEnd = std::min(xle, (int32_t)255);
            testSpan(xStart, xEnd);
        }

        // Test polygon interior
        if (drawInterior) {
            const int32_t xStart = std::max(xle + 1, 0);
            const int32_t xEnd = std::min(xrs, (int32_t)256);
            testSpan(xStart, xEnd - 1);
        }

        // Test right edge
        if (drawRightEdge) {
            const int32_t xStart = std::max(xrs, 0);
            const int32_t xEnd = std::min(xre, (int32_t)255);
            testSpan(xStart, xEnd);
        }

        // if (hadMismatch) {
        auto &cvl = data.verts[cvlIndex];
        auto &nvl = data.verts[nvlIndex];
        auto &ls = *leftSlope;
        auto &li = ls.Interp();

        auto &cvr = data.verts[cvrIndex];
        auto &nvr = data.verts[nvrIndex];
        auto &rs = *rightSlope;
        auto &ri = rs.Interp();

        fmt::print("    Left:  edge={:>3d}x{:<3d}..{:>3d}x{:<3d}  texcoords={:>4d}x{:<4d}   slope={}{} "
                   "{:>3d}..{:<3d} ({:>10d}..{:<10d}) DX={:<10d}   interp={:>3d}..{:<3d} ({})\n",
                   cvl[0], cvl[1], nvl[0], nvl[1], sl, tl, (ls.IsNegative() ? '-' : '+'), (ls.IsXMajor() ? 'X' : 'Y'),
                   ls.XStart(y), ls.XEnd(y), ls.FracXStart(y), ls.FracXEnd(y), ls.DX(), li.X0(), li.X1(), li.XMax());

        fmt::print("    Right: edge={:>3d}x{:<3d}..{:>3d}x{:<3d}  texcoords={:>4d}x{:<4d}   slope={}{} "
                   "{:>3d}..{:<3d} ({:>10d}..{:<10d}) DX={:<10d}   interp={:>3d}..{:<3d} ({})\n",
                   cvr[0], cvr[1], nvr[0], nvr[1], sr, tr, (rs.IsNegative() ? '-' : '+'), (rs.IsXMajor() ? 'X' : 'Y'),
                   rs.XStart(y), rs.XEnd(y), rs.FracXStart(y), rs.FracXEnd(y), rs.DX(), ri.X0(), ri.X1(), ri.XMax());

        fmt::print("    X interpolator: {:>3d}..{:<3d} ({})\n", xInterp.X0(), xInterp.X1(), xInterp.XMax());
        //}
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
