#include "halfedge.h"

#include <algorithm>
#include <array>
#include <cassert>

HalfedgeMesh HalfedgeFromIndexedTriangles(
    const float* srcVertexXYZs, int numSrcVertices,
    const uint32_t* srcTriangles, int numSrcTriangles)
{
#ifdef _DEBUG
    {
        // positions should have perfect reuse, otherwise you'll have invisible seams
        std::vector<std::array<float,3>> positions((std::array<float,3>*)srcVertexXYZs, (std::array<float, 3>*)srcVertexXYZs + numSrcVertices);
        std::sort(begin(positions), end(positions));
        for (int i = 0; i < (int)positions.size() - 1; i++)
        {
            assert(positions[i] != positions[i + 1]);
        }
    }
#endif

    HalfedgeMesh dstHalfedge;

    // Initialize vertices
    dstHalfedge.Vertices.resize(numSrcVertices, HalfedgeMesh::Vertex{ -1 });

    // Initialize faces
    dstHalfedge.Faces.resize(numSrcTriangles, HalfedgeMesh::Face{ -1 });

    // grab all full edges
    struct Edge
    {
        int VertexID0;
        int VertexID1;
        int AdjTriID;
    };
    std::vector<Edge> edges;
    for (int srcTriID = 0; srcTriID < numSrcTriangles; srcTriID++)
    {
        const uint32_t* tri = &srcTriangles[srcTriID * 3];

        // add all edges of this triangle
        for (int e = 0; e < 3; e++)
        {
            edges.push_back(Edge{ (int)tri[e], (int)tri[(e + 1) % 3], (int)srcTriID });
        }
    }

    auto EdgeVertexCmp = [](Edge e0, Edge e1)
    {
        // ensure unique order
        if (e0.VertexID0 > e0.VertexID1) std::swap(e0.VertexID0, e0.VertexID1);
        if (e1.VertexID0 > e1.VertexID1) std::swap(e1.VertexID0, e1.VertexID1);

        // sort by edge vertices then by adjacent triangle
        if (e0.VertexID0 < e1.VertexID0) return true;
        if (e0.VertexID0 > e1.VertexID0) return false;

        if (e0.VertexID1 < e1.VertexID1) return true;
        if (e0.VertexID1 > e1.VertexID1) return false;

        if (e0.AdjTriID < e1.AdjTriID) return true;
        if (e0.AdjTriID > e1.AdjTriID) return false;

        return false;
    };

    // sort the edges to find all the adjacent faces to each edge
    std::sort(begin(edges), end(edges), EdgeVertexCmp);

    std::vector<int> boundaryHalfedgeIDs;
    {
        // the edges list is currently missing boundary edges,
        // which means not every edge is duplicated in the list.
        // to make this so, find and merge in missing boundary edges.
        std::vector<Edge> boundaryEdges;
        for (size_t edgeIdx = 0; edgeIdx < edges.size(); )
        {
            size_t nextEdgeIdx;
            for (nextEdgeIdx = edgeIdx + 1; nextEdgeIdx < edges.size(); nextEdgeIdx++)
            {
                bool sameEdge =
                    (edges[edgeIdx].VertexID0 == edges[nextEdgeIdx].VertexID0 && edges[edgeIdx].VertexID1 == edges[nextEdgeIdx].VertexID1) ||
                    (edges[edgeIdx].VertexID1 == edges[nextEdgeIdx].VertexID0 && edges[edgeIdx].VertexID0 == edges[nextEdgeIdx].VertexID1);
                if (!sameEdge)
                {
                    break;
                }
            }

            // no more than 2 faces can share the same edge
            assert(nextEdgeIdx - edgeIdx <= 2);

            if (nextEdgeIdx - edgeIdx == 1)
            {
                // found boundary edge
                boundaryEdges.push_back(Edge{ edges[edgeIdx].VertexID1, edges[edgeIdx].VertexID0, -1 });
            }

            // move on to processing the next edge
            edgeIdx = nextEdgeIdx;
        }

        if (!boundaryEdges.empty())
        {
            std::vector<Edge> mergedEdges(edges.size() + boundaryEdges.size());
            boundaryHalfedgeIDs.resize(boundaryEdges.size());

            // merge the two lists of edges (non-boundary edges and boundary edges)
            // while keeping track of the halfedgeIDs of the boundary edges
            int edgeIdx = 0;
            int boundaryEdgeIdx = 0;
            int mergedEdgeIdx = 0;
            int boundaryHalfedgeIdx = 0;
            while (mergedEdgeIdx < (int)mergedEdges.size())
            {
                bool mergeFromBoundaryEdges;
                if (edgeIdx < (int)edges.size() && boundaryEdgeIdx < (int)boundaryEdges.size())
                {
                    if (EdgeVertexCmp(edges[edgeIdx], boundaryEdges[boundaryEdgeIdx]))
                    {
                        mergeFromBoundaryEdges = false;
                    }
                    else
                    {
                        mergeFromBoundaryEdges = true;
                    }
                }
                else if (edgeIdx < (int)edges.size())
                {
                    mergeFromBoundaryEdges = false;
                }
                else
                {
                    mergeFromBoundaryEdges = true;
                }

                Edge edgeToMerge;
                if (mergeFromBoundaryEdges)
                {
                    edgeToMerge = boundaryEdges[boundaryEdgeIdx];
                    boundaryEdgeIdx++;

                    boundaryHalfedgeIDs[boundaryHalfedgeIdx] = mergedEdgeIdx;
                    boundaryHalfedgeIdx++;
                }
                else
                {
                    edgeToMerge = edges[edgeIdx];
                    edgeIdx++;
                }

                mergedEdges[mergedEdgeIdx] = edgeToMerge;
                mergedEdgeIdx++;
            }

            edges = std::move(mergedEdges);
        }
    }

    dstHalfedge.Halfedges.resize(edges.size());

    // build up halfedge data structure (except halfedge next/prev links)
    for (size_t edgeIdx = 0; edgeIdx < edges.size(); edgeIdx += 2)
    {
        HalfedgeMesh::Halfedge* halfs = &dstHalfedge.Halfedges[edgeIdx];

        // Initial setup of the two halfedges
        halfs[0].VertexID = edges[edgeIdx].VertexID1;
        halfs[0].FaceID = edges[edgeIdx].AdjTriID;
        halfs[1].VertexID = edges[edgeIdx + 1].VertexID1;
        halfs[1].FaceID = edges[edgeIdx + 1].AdjTriID;

        // update associated face/vertex if this is the first adjacent halfedge
        for (int halfIdx = 0; halfIdx < 2; halfIdx++)
        {
            int otherHalfIdx = (halfIdx + 1) % 2;
            HalfedgeMesh::Halfedge* half = &halfs[halfIdx];
            HalfedgeMesh::Halfedge* otherHalf = &halfs[otherHalfIdx];

            if (dstHalfedge.Vertices[half->VertexID].HalfedgeID == -1)
            {
                // Vertices store the *outgoing* halfedge, so use the otherHalf
                dstHalfedge.Vertices[half->VertexID].HalfedgeID = (int)(otherHalf - &dstHalfedge.Halfedges[0]);
            }
            if (half->FaceID != -1 && dstHalfedge.Faces[half->FaceID].HalfedgeID == -1)
            {
                dstHalfedge.Faces[half->FaceID].HalfedgeID = (int)(half - &dstHalfedge.Halfedges[0]);
            }
        }
    }

    // the only thing missing now is the next/prev links for the halfedges

    // first, insert the next/prev links for the non-boundary edges
    for (int triID = 0; triID < numSrcTriangles; triID++)
    {
        // find edge indices for this triangle
        int edgeIDs[3];
        for (int e = 0; e < 3; e++)
        {
            Edge edgeToFind{ (int)srcTriangles[triID * 3 + e], (int)srcTriangles[triID * 3 + (e + 1) % 3], (int)triID };
            std::vector<Edge>::iterator it = std::lower_bound(begin(edges), end(edges), edgeToFind, EdgeVertexCmp);
            edgeIDs[e] = (int)std::distance(begin(edges), it);
        }

        // find associated half-edges for this triangle
        int halfEdgeIDs[3];
        for (int e = 0; e < 3; e++)
        {
            for (int h = 0; h < 2; h++)
            {
                if (dstHalfedge.Halfedges[edgeIDs[e] + h].FaceID == triID)
                {
                    halfEdgeIDs[e] = edgeIDs[e] + h;
                    break;
                }
            }
        }

        // associate the next/prev indices for the half-edges
        for (int e = 0; e < 3; e++)
        {
            dstHalfedge.Halfedges[halfEdgeIDs[e]].NextHalfedgeID = halfEdgeIDs[(e + 1) % 3];
            dstHalfedge.Halfedges[halfEdgeIDs[e]].PrevHalfedgeID = halfEdgeIDs[(e + 2) % 3];
        }
    }

    // finally, need to hook up next/prev links for boundary half-edges
    for (int boundaryHalfedgeID : boundaryHalfedgeIDs)
    {
        // rotate around the vertex to find the previous half-edge of the boundary
        int prevHalfedgeID = boundaryHalfedgeID ^ 1;
        while (dstHalfedge.Halfedges[prevHalfedgeID].FaceID != -1)
        {
            // next then opposite (rotate by one triangle clockwise)
            prevHalfedgeID = dstHalfedge.Halfedges[prevHalfedgeID].NextHalfedgeID;
            prevHalfedgeID ^= 1;
        }
        dstHalfedge.Halfedges[boundaryHalfedgeID].PrevHalfedgeID = prevHalfedgeID;
        dstHalfedge.Halfedges[prevHalfedgeID].NextHalfedgeID = boundaryHalfedgeID;
    }

    // ok, all done!
    return dstHalfedge;
}