#pragma once

#include "polygon.h"

struct PolygonState {
    Polygon *polygon;
    bool tall; // Polygon occupies more than one scanline

    struct Edge {
        size_t currVertex;
        size_t nextVertex;
        size_t vertexCount;
        bool nextPositive;

        Slope slope;

        void Setup(const Polygon &polygon, bool right) {
            // Initialize vertex pointers
            nextVertex = polygon.topVtxIndex;
            vertexCount = polygon.vertexCount;
            nextPositive = (right == polygon.windingCW);
            AdvanceVertex();

            // Initialize slope
            CalcSlope(polygon, polygon.topScreenY);
        }

        void CalcSlope(const Polygon &polygon, const int32_t y) {
            // Find first vertex at or below the current Y screen coordinate; stop if it's the last vertex
            while (currVertex != polygon.btmVtxIndex && y >= polygon.verts[nextVertex][1]) {
                AdvanceVertex();
            }

            const auto &cv = polygon.verts[currVertex];
            const auto &nv = polygon.verts[nextVertex];
            slope.Setup(cv[0], cv[1], nv[0], nv[1], 1, 1);
        }

        static size_t IncrementIndex(size_t index, size_t count) {
            return (index == count - 1) ? 0 : (index + 1);
        }

        static size_t DecrementIndex(size_t index, size_t count) {
            return (index == 0) ? (count - 1) : (index - 1);
        }

        static size_t NextIndex(size_t index, size_t count, bool nextPositive) {
            return nextPositive ? IncrementIndex(index, count) : DecrementIndex(index, count);
        }

        void AdvanceVertex() {
            currVertex = nextVertex;
            nextVertex = NextIndex(currVertex, vertexCount, nextPositive);
        }
    };

    void CalcSlopes(const int32_t y) {
        if (leftEdge.currVertex != polygon->btmVtxIndex && y >= polygon->verts[leftEdge.nextVertex][1]) {
            leftEdge.CalcSlope(*polygon, y);
        }
        if (rightEdge.currVertex != polygon->btmVtxIndex && y >= polygon->verts[rightEdge.nextVertex][1]) {
            rightEdge.CalcSlope(*polygon, y);
        }
    }

    void SetupPolygon(Polygon &polygon) {
        this->polygon = &polygon;
        tall = (polygon.topScreenY != polygon.btmScreenY);

        leftEdge.Setup(polygon, false);
        rightEdge.Setup(polygon, true);
    }

    Edge leftEdge;
    Edge rightEdge;
};
