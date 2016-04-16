#pragma once

#include <vector>
#include <cstdint>

struct HalfedgeMesh
{
    struct Vertex
    {
        int HalfedgeID; // A halfedge coming out of this vertex
        // Vertex 3D positions are to be stored separately
    };

    struct Face
    {
        int HalfedgeID; // One halfedge this face is adjacent to
    };

    struct Halfedge
    {
        int VertexID; // Vertex this halfedge is pointing to
        int FaceID; // Face this halfedge is adjacent to
        int NextHalfedgeID;
        int PrevHalfedgeID;

        // Opposite halfedges are adjacent in the array:
        // oppositeID = myID ^ 1
    };

    std::vector<Vertex> Vertices;
    std::vector<Face> Faces;
    std::vector<Halfedge> Halfedges;
};

HalfedgeMesh HalfedgeFromIndexedTriangles(
    const float* srcVertexXYZs, int numSrcVertices,
    const uint32_t* srcTriangles, int numSrcTriangles);