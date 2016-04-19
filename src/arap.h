#pragma once

// Recommended number of iterations by Sorkine, O.
#define DEFAULT_NUM_ARAP_ITERATIONS 3

// convenience for knowing how many floats are needed for a system matrix for N vertices
// You need at most ARAP_PACKED_SYSTEM_SIZE_IN_FLOATS(N) floats for N vertices
// You need at least ARAP_PACKED_SYSTEM_SIZE_IN_FLOATS(NFree) floats for a system with N unconstrained vertices
// Pick whatever lower bound is most convenient...
#define ARAP_PACKED_SYSTEM_SIZE_IN_FLOATS(n) (((n) * ((n) + 1)) / 2)

// The system matrix only needs to get built and factorized once up front.
void arap_factorize_system(
    int numVertices, // Total number of vertices in the mesh
    const int* vertexConstraintStatuses, // 0 == not constrained, 1 === constrained
    const int* vertexIncidentHalfedgeIDs, // the halfedgeID incident to each vertex
    const int* halfedgeVFNPs, // (vertexID, faceID, nextHalfedgeID, prevHalfedgeID), and opposite = halfedgeID ^ 1
    const float* edgeWeights, // one for every 2 halfedges
    float* factorizedPackedSystemMatrix); // storage for the output. See ARAP_PACKED_SYSTEM_SIZE_IN_FLOATS for size requirements.

void arap(
    int numVertices, // Total number of vertices in the mesh
    const float* vertexBindPosePositionXYZs, // original rest pose of the mesh
    float* vertexInitialGuessPositionXYZs, // current guess. for constrained vertices, set their constrained positions in the guess.
    const int* vertexConstraintStatuses, // 0 == not constrained, 1 === constrained
    const int* vertexIncidentHalfedgeIDs, // the halfedgeID incident to each vertex
    const int* halfedgeVFNPs, // (vertexID, faceID, nextHalfedgeID, prevHalfedgeID), and opposite = halfedgeID ^ 1
    const float* edgeWeights, // one for every 2 halfedges
    const float* factorizedPackedSystemMatrix, // created from arap_factorize_system
    int numIterations); // number of iterations of the algorithm (more iterations = more rigid). Recommended DEFAULT_NUM_ARAP_ITERATIONS.
