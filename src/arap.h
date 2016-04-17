#pragma once

#define DEFAULT_NUM_ARAP_ITERATIONS 3

// convenience for knowing how many floats are needed for a system matrix for N vertices
#define ARAP_PACKED_SYSTEM_SIZE(n) (n * (n + 1) / 2)

// The system matrix only needs to get built and factorized once up front.
void arap_factorize_system(
    int numVertices, int numEdges,
    const int* halfedgeVFNPs, // (vertexID, faceID, nextHalfedgeID, prevHalfedgeID), and opposite = halfedgeID ^ 1
    const float* edgeWeights, // one for every 2 halfedges
    float* factorizedPackedSystemMatrix);

void arap(
    const float* vertexBindPosePositionXYZs, float* vertexInitialGuessPositionXYZs, int numVertices,
    const int* vertexIncidentHalfedgeIDs,
    const int* halfedgeVFNPs, // (vertexID, faceID, nextHalfedgeID, prevHalfedgeID), and opposite = halfedgeID ^ 1
    const float* edgeWeights, // one for every 2 halfedges
    const float* factorizedPackedSystemMatrix, // created from arap_factorize_system
    int* constrainedVertexIDs, int numConstrainedVertexIDs,
    int numIterations);