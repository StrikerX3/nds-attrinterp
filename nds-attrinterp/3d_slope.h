#pragma once

#include "3d_interp.h"

#include <cstdint>
#include <utility>

/// <summary>
/// Computes 3D rasterization slopes based on Nintendo DS's hardware interpolation.
/// </summary>
/// <remarks>
/// The algorithm implemented by this class produces pixel-perfect slopes matching the Nintendo DS's 3D interpolator.
///
/// The hardware uses 32-bit integers with 18-bit fractional parts throughout the interpolation process, with one
/// notable exception in X-major slopes.
///
/// To calculate the X increment per scanline (DX), the hardware first computes the reciprocal of Y1-Y0 then multiplies
/// the result by X1-X0. This order of operations avoids a multiplication overflow at the cost of precision on the
/// division.
///
/// For X-major lines, the interpolator produces line spans for each scanline. The start of the span is calculated by
/// first offseting the Y coordinate to the Y0-Y1 range (subtracting Y0 from Y), then multiplying the offset Y by DX,
/// adding the X0 offset and finally a +0.5 bias. The end on the span is computed based on its starting coordinate,
/// discarding (masking out) the 9 least significant bits (which could be seen as rounding down, or the floor function),
/// then adding DX and subtracting 1.0. The exact algorithm is unknown, but it is possible that the starting coordinate
/// is shifted right by 9 for an intermediate calculation then shifted back left by 9 to restore the fractional part.
///
/// The formulae for determining the starting and ending X coordinates of a span of an X-major slope are:
///
///    DX = 1 / (Y1 - Y0) * (X1 - X0)
///    Xstart = (Y - Y0) * DX + X0 + 0.5
///    Xend = Xstart[discarding 9 LSBs] + DX - 1.0
///
/// Due to the 9 LSBs being discarded, certain X-major slopes (such as 69x49, 70x66, 71x49 and more) display a one-pixel
/// gap on hardware. This is calculated accurately with the formulae above.
///
/// Y-major slopes contain only one pixel per scanline. The formula for interpolating the X coordinate based on the Y
/// coordinate is very similar to that of the X-major interpolation, with the only difference being that the +0.5 bias
/// is not applied.
///
///    X = (Y - Y0) * DX + X0
///
/// Note that there is no need to compute a span as there's always going to be only one pixel per scanline. Also, there
/// are no one-pixel gaps on Y-major lines since the Nintendo DS's rasterizer is scanline-based and the interpolation is
/// computed on every scanline.
///
/// Negative slopes work in a similar fashion. In fact, negative slopes perfectly match their positive counterparts down
/// to the one-pixel gaps which happen in exactly the same spots. The gaps in negative slopes are to the left of a span,
/// while in positive slopes the gaps are to the right of a span, as shown below (rows 35 to 39 of 69x49 slopes):
///
///    Positive slope        Negative slope
///      ##  +---- mind the gap ----+  ##
///        # |                      | #
///         #V                      V#
///           #                    #
///            #                  #
///
/// The behavior for negative slopes is implemented in this class in the following manner, compared to positive slopes:
/// - The raw value of X0 coordinate used to compute the starting X coordinate of a span is decremented by 1, that is,
///   the value is subtracted an amount equal to 1.0 / 2^fractionalBits
/// - X0 and X1 are swapped; as a consequence, DX remains positive
/// - The starting X coordinate is the span's rightmost pixel; conversely, the ending X coordinate is its leftmost pixel
/// - The starting X coordinate is decremented by the computed Y*DX displacement instead of incremented
/// - The nine least significant bits of the ending X coordinate are rounded up the to the largest number less than 1.0
///   (511 as a raw integer with 9 fractional bits)
///
/// All other operations are otherwise identical.
/// </remarks>
class Slope {
public:
    /// <summary>
    /// The number of fractional bits (aka resolution) of the interpolator.
    /// </summary>
    /// <remarks>
    /// The Nintendo DS uses 18 fractional bits for interpolation.
    /// </remarks>
    static constexpr uint32_t kFracBits = 18;

    /// <summary>
    /// The value 1.0 with fractional bits.
    /// </summary>
    static constexpr uint32_t kOne = (1 << kFracBits);

    /// <summary>
    /// The bitmask for the fractional bits.
    /// </summary>
    static constexpr uint32_t kFracMask = (kOne - 1);

    /// <summary>
    /// The bias applied to the interpolation of X-major spans.
    /// </summary>
    static constexpr uint32_t kBias = (kOne >> 1);

    /// <summary>
    /// The mask applied during interpolation of X-major spans, removing half of the least significant fractional bits
    /// (rounded down).
    /// </summary>
    static constexpr uint32_t kMask = (~0u << (kFracBits / 2));

    /// <summary>
    /// Configures the slope to interpolate the line (X0,X1)-(Y0,Y1) using screen coordinates.
    /// </summary>
    /// <param name="x0">First X coordinate</param>
    /// <param name="y0">First Y coordinate</param>
    /// <param name="x1">Second X coordinate</param>
    /// <param name="y1">Second Y coordinate</param>
    /// <param name="w0">W value at the first coordinate</param>
    /// <param name="w1">W value at the second coordinate</param>
    void Setup(int32_t x0, int32_t y0, int32_t x1, int32_t y1, int32_t w0, int32_t w1) {
        // Always interpolate top to bottom
        if (y1 < y0) {
            std::swap(x0, x1);
            std::swap(y0, y1);
        }

        // Store reference coordinates
        m_x0 = x0 << kFracBits;
        m_y0 = y0;

        // Determine if this is a negative slope and adjust accordingly
        m_negative = (x1 < x0);
        if (m_negative) {
            m_x0--;
            std::swap(x0, x1);
        }

        // Compute coordinate deltas and determine if the slope is X-major
        int32_t dx = (x1 - x0);
        int32_t dy = (y1 - y0);
        m_xMajor = (dx > dy);

        // Precompute bias for X-major or diagonal slopes
        if (m_xMajor || dx == dy) {
            if (m_negative) {
                m_x0 -= kBias;
            } else {
                m_x0 += kBias;
            }
        }

        // Compute X displacement per scanline
        m_dx = dx;
        if (dy != 0) {
            m_dx *= kOne / dy; // This ensures the division is performed before the multiplication
        } else {
            m_dx *= kOne;
        }

        // Configure attribute/depth interpolator
        if (m_xMajor) {
            if (m_negative) {
                m_interp.Setup(x1, x0, w0, w1);
                // m_interp.Setup(XEnd(y0), XStart(y1 - 1), w0, w1);
                // m_interp.Setup(FracXEnd(y0), FracXEnd(y1), w0, w1);
            } else {
                m_interp.Setup(x0, x1, w0, w1);
                // m_interp.Setup(XStart(y0), XEnd(y1 - 1), w0, w1);
                // m_interp.Setup(FracXStart(y0), FracXStart(y1), w0, w1);
            }
        } else {
            m_interp.Setup(y0, y1, w0, w1);
        }
    }

    /// <summary>
    /// Update interpolation factors for the specified Y coordinate.
    /// </summary>
    /// <param name="y">The Y coordinate</param>
    /// <param name="end">true to compute factors at the end of the span, false for the start of the span. Only used
    /// with X-major slopes.</param>
    void ComputeFactors(int32_t y, bool end) {
        if (m_xMajor) {
            // X coordinate to use depends on the end flag and m_negative:
            //    neg  end  !=   coordinate
            //    no   no   no   start
            //    no   yes  yes  end
            //    yes  no   yes  end
            //    yes  yes  no   start
            // int32_t x = (!m_negative && end) ? XStart(y + 1) : (m_negative && !end) ? XEnd(y) : XStart(y);
            int32_t x = (m_negative != end) ? XEnd(y) : XStart(y);
            // int32_t x = (m_negative != end) ? FracXEnd(y) : FracXStart(y);
            m_interp.ComputeFactors(x);
        } else {
            m_interp.ComputeFactors(y);
        }
    }

    /// <summary>
    /// Compuptes perspective-correct interpolation of the attribute using the current interpolation factors.
    /// </summary>
    /// <param name="a0">The attribute value at interpolation factor 0</param>
    /// <param name="a1">The attribute value at interpolation factor 1</param>
    /// <returns>The interpolated attribute value</returns>
    int32_t InterpolateAttribute(int32_t a0, int32_t a1) const {
        return m_interp.InterpolateAttribute(a0, a1);
    }

    /// <summary>
    /// Interpolates the depth value using the appropriate method -- perspective-correct interpolation for W or linear
    /// interpolation for Z.
    /// </summary>
    /// <param name="z0">The depth value at interpolation factor 0</param>
    /// <param name="z1">The depth value at interpolation factor 0</param>
    /// <param name="w">true if the depth value is W, false if it is Z</param>
    /// <returns>The interpolated depth value</returns>
    int32_t InterpolateDepth(int32_t z0, int32_t z1, bool w) const {
        return m_interp.InterpolateDepth(z0, z1, w);
    }

    /// <summary>
    /// Computes the starting position of the span at the specified Y coordinate, including the fractional part.
    /// </summary>
    /// <param name="y">The Y coordinate, which must be between Y0 and Y1 specified in Setup.</param>
    /// <returns>The starting X coordinate of the specified scanline's span</returns>
    constexpr int32_t FracXStart(int32_t y) const {
        int32_t displacement = (y - m_y0) * m_dx;
        if (m_negative) {
            return m_x0 - displacement;
        } else {
            return m_x0 + displacement;
        }
    }

    /// <summary>
    /// Computes the ending position of the span at the specified Y coordinate, including the fractional part.
    /// </summary>
    /// <param name="y">The Y coordinate, which must be between Y0 and Y1 specified in Setup.</param>
    /// <returns>The ending X coordinate of the specified scanline's span</returns>
    constexpr int32_t FracXEnd(int32_t y) const {
        int32_t result = FracXStart(y);
        if (m_xMajor) {
            if (m_negative) {
                // The bit manipulation sequence (~mask - (x & ~mask)) acts like a ceiling function.
                // Since we're working in the opposite direction here, the "floor" is actually the ceiling.
                result = result + (~kMask - (result & ~kMask)) - m_dx + kOne;
            } else {
                result = (result & kMask) + m_dx - kOne;
            }
        }
        return result;
    }

    /// <summary>
    /// Computes the starting position of the span at the specified Y coordinate as a screen coordinate (dropping the
    /// fractional part).
    /// </summary>
    /// <param name="y">The Y coordinate, which must be between Y0 and Y1 specified in Setup.</param>
    /// <returns>The starting X screen coordinate of the scanline's span</returns>
    constexpr int32_t XStart(int32_t y) const {
        return FracXStart(y) >> kFracBits;
    }

    /// <summary>
    /// Computes the ending position of the span at the specified Y coordinate as a screen coordinate (dropping the
    /// fractional part).
    /// </summary>
    /// <param name="y">The Y coordinate, which must be between Y0 and Y1 specified in Setup.</param>
    /// <returns>The ending X screen coordinate of the scanline's span</returns>
    constexpr int32_t XEnd(int32_t y) const {
        return FracXEnd(y) >> kFracBits;
    }

    /// <summary>
    /// Retrieves the X coordinate increment per scanline.
    /// </summary>
    /// <returns>The X displacement per scanline (DX)</returns>
    constexpr int32_t DX() const {
        return m_dx;
    }

    /// <summary>
    /// Determines if the slope is X-major.
    /// </summary>
    /// <returns>true if the slope is X-major.</returns>
    constexpr bool IsXMajor() const {
        return m_xMajor;
    }

    /// <summary>
    /// Determines if the slope is negative (i.e. X decreases as Y increases).
    /// </summary>
    /// <returns>true if the slope is negative.</returns>
    constexpr bool IsNegative() const {
        return m_negative;
    }

    const Interpolator &Interp() const {
        return m_interp;
    }

private:
    int32_t m_x0;    // X0 coordinate (minus 1 if this is a negative slope)
    int32_t m_y0;    // Y0 coordinate
    int32_t m_dx;    // X displacement per scanline
    bool m_negative; // True if the slope is negative (X1 < X0)
    bool m_xMajor;   // True if the slope is X-major (X1-X0 > Y1-Y0)

    Interpolator m_interp; // Attribute interpolator
};
