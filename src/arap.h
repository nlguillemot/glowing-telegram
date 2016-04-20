#pragma once

// Recommended number of iterations by Sorkine, O.
#define DEFAULT_NUM_ARAP_ITERATIONS 3

struct arap_system;

// The system matrix only needs to get built and factorized once up front.
arap_system* create_arap_system_matrix(
    int numVertices, // Total number of vertices in the mesh
    const int* vertexConstraintStatuses, // 0 == not constrained, 1 === constrained
    const int* vertexIncidentHalfedgeIDs, // the halfedgeID incident to each vertex
    const int* halfedgeVFNPs, // (vertexID, faceID, nextHalfedgeID, prevHalfedgeID), and opposite = halfedgeID ^ 1
    const float* edgeWeights); // one for every 2 halfedges

void destroy_arap_system_matrix(arap_system* sys);

void arap(
    arap_system* systemMatrix, // from create_arap_system_matrix
    const float* vertexBindPosePositionXYZs, // original rest pose of the mesh
    float* vertexInitialGuessPositionXYZs, // current guess. for constrained vertices, set their constrained positions in the guess.
    const int* vertexIncidentHalfedgeIDs, // the halfedgeID incident to each vertex
    const int* halfedgeVFNPs, // (vertexID, faceID, nextHalfedgeID, prevHalfedgeID), and opposite = halfedgeID ^ 1
    const float* edgeWeights, // one for every 2 halfedges
    int numIterations); // number of iterations of the algorithm (more iterations = more rigid). Recommended DEFAULT_NUM_ARAP_ITERATIONS.
