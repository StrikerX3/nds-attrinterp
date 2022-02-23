#pragma once

#include <cstdint>

// According the article https://melonds.kuribo64.net/comments.php?id=85, the interpolation method used by the NDS is
// based on the typical perspective-correct formula, simplified to reduce the number of divisions and using 0..1 as a
// generic range for the attributes. The simplified, "less division-y" formula for perspective-correct interpolation is:
//
//           ((xmax - x) * A0 * W1) + (x * A1 * W0)
//     Ax = ----------------------------------------
//                ((xmax - x) * W1) + (x * W0)
//
// where:
// - x is the position within the span, units don't matter (as long as they're consistent with xmax)
// - xmax is the length of the span
// - A0 and A1 are the attributes at each end of the span
// - W0 and W1 are the associated W coordinates
// - Ax is the interpolated attribute at position x
//
// The NDS goes a step further and uses A0 = 0 and A1 = 1, so that the interpolation function spans the unit range,
// effectively becoming an "interpolation factor" function. The resulting formula is:
//
//                         (x * W0)
//     factor = ------------------------------
//               (x * W0) + ((xmax - x) * W1)
//
// The interpolator also uses a simple linear interpolation of the x parameter when W0 and W1 are equal and have their
// least significant bits clear, mostly to help with 2D rendering scenarios where the added work would not only be
// unnecessary, but also produce visible imperfections.
class Interpolator {
public:
    Interpolator() {
        Setup(0, 0, 0, 0);
    }

    void Setup(int32_t x0, int32_t x1, int32_t w0, int32_t w1) {
        m_x0 = x0;
        m_x1 = x1;
        m_xmax = x1 - x0;

        m_w0 = w0;
        m_w1 = w1;
    }

    // Computes the interpolation factors for the given X position within the range X0..X1
    void ComputeFactors(int32_t x) {
        if (m_xmax != 0) {
            x -= m_x0;
            m_x = x;
            // Compute numerator and denominator for perspective-correct interpolation
            m_perspNum = (int64_t)x * m_w0;
            m_perspDen = m_perspNum + (int64_t)(m_xmax - x) * m_w1;
        }
    }

    // Interpolates the Z or W coordinate in the range Z0..Z1
    int32_t InterpolateDepth(int32_t z0, int32_t z1, bool w) const {
        if (m_xmax == 0) {
            return z0; // No need to interpolate an empty range
        }
        if (z0 == z1) {
            return z0; // No need to interpolate an unchanging value
        }

        if (w) {
            // W-buffering uses perspective-correct interpolation
            return LerpPerspective(z0, z1);
        } else {
            // Z-buffering uses linear interpolation
            return LerpLinear(z0, z1);
        }
    }

    // Interpolates the attribute in the range A0..A1
    int32_t InterpolateAttribute(int32_t a0, int32_t a1) const {
        if (m_xmax == 0) {
            return a0; // No need to interpolate an empty range
        }
        if (a0 == a1) {
            return a0; // No need to interpolate an unchanging value
        }
        return LerpPerspective(a0, a1);
    }

    int32_t X0() const {
        return m_x0;
    }

    int32_t X1() const {
        return m_x1;
    }

    int32_t XMax() const {
        return m_xmax;
    }

    int64_t Num() const {
        return m_perspNum;
    }

    int64_t Den() const {
        return m_perspDen;
    }

private:
    int32_t m_x0;   // X0: initial point of the interpolation range
    int32_t m_x1;   // X1: final point of the interpolation range
    int32_t m_xmax; // X1 - X0 (xmax in the formula)
    int32_t m_x;    // Current X value (x in the formula)

    int32_t m_w0; // W0 value
    int32_t m_w1; // W1 value

    int64_t m_perspNum = 0; // Perspective-correct interpolation numerator
    int64_t m_perspDen = 0; // Perspective-correct interpolation denominator

    int32_t LerpLinear(int32_t a0, int32_t a1) const {
        return a0 + (int64_t)(a1 - a0) * m_x / m_xmax;
    }

    int32_t LerpPerspective(int32_t a0, int32_t a1) const {
        if (a0 > a1) {
            return a1 + (((int64_t)(a0 - a1) * (m_perspDen - m_perspNum)) / m_perspDen);
        } else {
            return a0 + (int64_t)(a1 - a0) * m_perspNum / m_perspDen;
        }
    }
};
