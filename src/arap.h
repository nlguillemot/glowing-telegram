#pragma once

#define DEFAULT_NUM_ARAP_ITERATIONS 3

void arap(
    const float* vertexBindPosePositionXYZs, float* vertexInitialGuessPositionXYZs, int numVertices,
    const int* vertexIncidentHalfedgeIDs,
    const int* halfedgeVFNPs, // (vertexID, faceID, nextHalfedgeID, prevHalfedgeID, opposite = halfedgeID ^ 1)
    const float* edgeWeights, // one fore every 2 halfedges
    int* constrainedVertexIDs, int numConstrainedVertexIDs,
    int numIterations);