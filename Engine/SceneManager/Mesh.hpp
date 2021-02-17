#pragma once
#include "portable.hpp"

namespace Panda
{
    typedef ENUM(_IndexSize)
    {
        kIndexSize16 = 0x00000010,
        kIndexSize32 = 0x00000020,
        kIndexSizePlatformDependent = 0x10000000,    // start point of platform specific value
    } IndexSize;

    typedef ENUM(_PrimitiveType)
    {
        kPrimitiveTypeNode                  = 0x00000000, ///< No particular primitive type.
        kPrimitiveTypePointList             = 0x00000001, ///< For N>=0, vertex N renders a point.
        kPrimitiveTypeLineList              = 0x00000002, ///< For N>=0, vertices [N*2+0, N*2+1] render a line.
        kPrimitiveTypeLineStrip             = 0x00000003, ///< For N>=0, vertices [N, N+1] render a line.
        kPrimitiveTypeTriList               = 0x00000004, ///< For N>=0, vertices [N*3+0, N*3+1, N*3+2] render a triangle.
        kPrimitiveTypeTriFan                = 0x00000005, ///< For N>=0, vertices [0, (N+1)%M, (N+2)%M] render a triangle, where M is the vertex count.
        kPrimitiveTypeTriStrip              = 0x00000006, ///< For N>=0, vertices [N*2+0, N*2+1, N*2+2] and [N*2+2, N*2+1, N*2+3] render triangles.
        kPrimitiveTypePatch                 = 0x00000009, ///< Used for tessellation.
        kPrimitiveTypeLineListAdjacency     = 0x0000000a, ///< For N>=0, vertices [N*4..N*4+3] render a line from [1, 2]. Lines [0, 1] and [2, 3] are adjacent to the rendered line.
        kPrimitiveTypeLineStripAdjacency    = 0x0000000b, ///< For N>=0, vertices [N+1, N+2] render a line. Lines [N, N+1] and [N+2, N+3] are adjacent to the rendered line.
        kPrimitiveTypeTriListAdjacency      = 0x0000000c, ///< For N>=0, vertices [N*6..N*6+5] render a triangle from [0, 2, 4]. Triangles [0, 1, 2] [4, 2, 3] and [5, 0, 4] are adjacent to the rendered triangle.
        kPrimitiveTypeTriStripAdjacency     = 0x0000000d, ///< For N>=0, vertices [N*4..N*4+6] render a triangle from [0, 2, 4] and [4, 2, 6]. Odd vertices Nodd form adjacent triangles with indices min(Nodd+1,Nlast) and max(Nodd-3,Nfirst).
        kPrimitiveTypeRectList              = 0x00000011, ///< For N>=0, vertices [N*3+0, N*3+1, N*3+2] render a screen-aligned rectangle. 0 is upper-left, 1 is upper-right, and 2 is the lower-left corner.
        kPrimitiveTypeLineLoop              = 0x00000012, ///< Like <c>kPrimitiveTypeLineStrip</c>, but the first and last vertices also render a line.
        kPrimitiveTypeQuadList              = 0x00000013, ///< For N>=0, vertices [N*4+0, N*4+1, N*4+2] and [N*4+0, N*4+2, N*4+3] render triangles.
        kPrimitiveTypeQuadStrip             = 0x00000014, ///< For N>=0, vertices [N*2+0, N*2+1, N*2+3] and [N*2+0, N*2+3, N*2+2] render triangles.
        kPrimitiveTypePolygon               = 0x00000015, ///< For N>=0, vertices [0, N+1, N+2] render a triangle.
    } PrimitiveType;

    class Mesh
    {
        public:
            Mesh() :
                m_VertexBuffer(0),
                m_VertexBufferSize(0),
                m_VertexCount(0),
                m_VertexAttributeCount(0),
                m_IndexBuffer(0),
                m_IndexBufferSize(0),
                m_IndexCount(0),
                m_IndexType(static_cast<uint32_t>(IndexSize::kIndexSize16))
                {}

        public:
            void*       m_VertexBuffer;
            uint32_t    m_VertexBufferSize;
            uint32_t    m_VertexCount;
            uint32_t    m_VertexAttributeCount;
            uint32_t    m_Reserved0[4];
            void*       m_IndexBuffer;
            uint32_t    m_IndexBufferSize;
            uint32_t    m_IndexCount;
            uint32_t    m_IndexType;
    };

    class SimpleMesh : public Mesh
    {
        public:
            SimpleMesh() : Mesh(),
                m_PrimitiveType(PrimitiveType::kPrimitiveTypeTriList) {}

            uint32_t        m_VertexStride;
            PrimitiveType   m_PrimitiveType;
            uint32_t        m_Reversed[2];
    };
}
