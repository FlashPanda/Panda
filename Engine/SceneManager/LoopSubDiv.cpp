#include "LoopSubDiv.hpp"
#include "Core/MemoryManager.hpp"
#include "SceneObjectTriangle.hpp"

namespace Panda
{
    struct SDFace;
    struct SDVertex;

// LoopSubdiv Macros
#define NEXT(i) (((i) + 1) % 3)
#define PREV(i) (((i) + 2) % 3)

    struct SDVertex{
        SDVertex(const Vector3Df& p = Vector3Df({0.f, 0.f, 0.f})) : p(p) {}

        int32_t Valence();
        void OneRing(Vector3Df* p);

        Vector3Df p;
        SDFace* startFace = nullptr;    // A face which contains the vertex. We can use it to iterate all faces.
                                        // The pointer will be managered by unique_ptr when it is created.
        SDVertex* child = nullptr;      // The vertex which is origin vertex modified by subdivide operations.
        bool regular = false;
        bool boundary = false;
    };

    struct SDFace 
    {
        SDFace()
        {
            for (int32_t i = 0; i < 3; ++i)
            {
                v[i] = nullptr;
                f[i] = nullptr;
            }
            for (int32_t i = 0; i < 4; ++i)
                children[i] = nullptr;
        }

        int32_t Vnum(SDVertex* vert) const
        {
            for (int32_t i = 0; i < 3; ++i)
                if (v[i] == vert) 
                    return i;
            return -1;
        }
        SDFace* NextFace(SDVertex* vert)
        {
            return f[Vnum(vert)];
        }
        SDFace* PrevFace(SDVertex* vert)
        {
            return f[PREV(Vnum(vert))];
        }
        SDVertex* NextVert(SDVertex* vert)
        {
            return v[NEXT(Vnum(vert))];
        }
        SDVertex* PrevVert(SDVertex* vert)
        {
            return v[PREV(Vnum(vert))];
        }
        SDVertex* OtherVert(SDVertex* v0, SDVertex* v1)
        {
            for (int32_t i = 0; i < 3; ++i)
                if (v[i] != v0 && v[i] != v1)
                    return v[i];
            return nullptr;
        }

        SDVertex* v[3];
        SDFace* f[3];
        SDFace* children[4];
    };

    struct SDEdge
    {
        SDEdge(SDVertex* v0 = nullptr, SDVertex* v1 = nullptr)
        {
            v[0] = (std::min)(v0, v1);
            v[1] = (std::max)(v0, v1);
            f[0] = f[1] = nullptr;
            f0EdgeNum = -1;
        }

        // Comparison Function
        // Used in std::map data.
        bool operator<(const SDEdge& e2) const
        {
            if (v[0] == e2.v[0]) return v[1] < e2.v[1];
            return v[0] < e2.v[0];
        }

        SDVertex* v[2];
        SDFace* f[2];
        int32_t f0EdgeNum;
    };

    inline int32_t SDVertex::Valence()
    {
        SDFace* f = startFace;
        if (!boundary)
        {
            // Compute valence of interior vertex
            int32_t nf = 1;
            while ((f = f->NextFace(this)) != startFace) ++nf;
            return nf;
        }
        else
        {
            // COmpute valence of boundary vertex
            int32_t nf = 1;
            while ((f = f->NextFace(this)) != nullptr) ++nf;
            f = startFace;
            while ((f = f->PrevFace(this)) != nullptr) ++nf;
            return nf + 1;
        }
    }

    void SDVertex::OneRing(Vector3Df* p)
    {
        if (!boundary)
        {
            // Get one-ring vertices for interior  vertex
            SDFace* face = startFace;
            do{
                *p++ = face->NextVert(this)->p;
                face = face->NextFace(this);
            }while (face != startFace);
        }
        else 
        {
            // Get one-ring vertices for boundary vertex
            SDFace* face = startFace, *f2;
            while ((f2 = face->NextFace(this)) != nullptr) face = f2;
            *p++ = face->NextVert(this)->p;
            do{
                *p++ = face->PrevVert(this)->p;
                face = face ->PrevFace(this);
            }while (face != nullptr);
        }
    }

    static Vector3Df WeightOneRing(SDVertex* vert, float beta)
    {
        // Put _vert_ one-ring in _pRing_
        int32_t valence = vert->Valence();
        Vector3Df* pRing = reinterpret_cast<Vector3Df*>(alloca(sizeof(Vector3Df) * valence));
        vert->OneRing(pRing);
        Vector3Df p = (1 - valence * beta) * vert->p;
        for (int32_t i = 0; i < valence; ++i) p += beta * pRing[i];
        return p;
    }

    static Vector3Df WeightBoundary(SDVertex* vert, float beta)
    {
        // Put _vert_ one-ring in _pRing_
        int32_t valence = vert->Valence();
        Vector3Df* pRing = reinterpret_cast<Vector3Df*>(sizeof(Vector3Df) * valence);
        vert->OneRing(pRing);
        Vector3Df p = (1 - 2 * beta) * vert->p;
        p += beta * pRing[0];
        p += beta * pRing[valence - 1];
        return p;
    }

    inline float Beta(int32_t valence)
    {
        if (valence == 3)
            return 3.f / 16.f;
        else 
            return 3.f / (8.f * valence);
    }

    inline float LoopGamma(int32_t valence)
    {
        return 1.f / (valence + 3.f / (8.f * Beta(valence)));
    }

    // LoopSubdiv Funciton Definitions
    static std::vector<std::shared_ptr<SceneObjectShape>> LoopSubDivide(
		const std::shared_ptr<Matrix4f>& objectToWorld,
		const std::shared_ptr<Matrix4f>& worldToObject,
        int32_t nLevels, int32_t nIndices,
        const int32_t*vertexIndices, int32_t nVertices, const Vector3Df* p)
    {
        /// The two are just pointers, without memory allocated.
        std::vector<SDVertex*> vertices;
        std::vector<SDFace*> faces;

        // Allocate _LoopSubDiv_ vertices and faces
        std::unique_ptr<SDVertex[]> verts(new SDVertex[nVertices]); // The new vertices.
        for (int32_t i = 0; i < nVertices; ++i)
        {
            verts[i] = SDVertex(p[i]);
            vertices.push_back(&verts[i]); // Save the pointers which we can easily use.
        }
        int32_t nFaces = nIndices / 3;
        std::unique_ptr<SDFace[]> fs(new SDFace[nFaces]); // The new faces.
        for (int32_t i = 0; i < nFaces; ++i)
            faces.push_back(&fs[i]); // Save the pointers which we can easily use.

        // Set face to vertex pointers.
        const int32_t* vp = vertexIndices;
        for (int32_t i = 0; i < nFaces; ++i, vp += 3)
        {
            SDFace* f = faces[i];
            for (int32_t j = 0; j < 3; ++j)
            {
                SDVertex* v = vertices[vp[j]];
                f->v[j] = v;
                v->startFace = f;
            }
        }

        // Set neighbor pointers in _faces_
        std::set<SDEdge> edges;
        for (int32_t i = 0; i < nFaces; ++i)
        {
            SDFace* f = faces[i];
            for (int32_t edgeNum = 0; edgeNum < 3; ++edgeNum)
            {
                // Update neighbor pointer for _edgeNum_
                int32_t v0 = edgeNum, v1 = NEXT(edgeNum);
                SDEdge e(f->v[v0], f->v[v1]);
                if (edges.find(e) == edges.end())
                {
                    // Handle new edge
                    e.f[0] = f;
                    e.f0EdgeNum = edgeNum;
                    edges.insert(e);
                }
                else 
                {
                    // Handle previously seen edge
                    e = *edges.find(e);
                    e.f[0]->f[e.f0EdgeNum] = f;
                    f->f[edgeNum] = e.f[0];
                    edges.erase(e);
                }
            }
        }

        // Finish vertex initialization
        for (int32_t i = 0; i < nVertices; ++i)
        {
            SDVertex* v = vertices[i];
            SDFace* f = v->startFace;
            do {
                f = f->NextFace(v);
            }while(f && f != v->startFace);
            v->boundary = (f == nullptr);
            if (!v->boundary && v->Valence() == 6)
                v->regular = true;
            else if (v->boundary && v->Valence() == 4)
                v->regular = true;
            else
                v->regular = false;
        }


        // Refine _LoopSubDiv_ into triangles
        std::vector<SDFace*> f = faces;
        std::vector<SDVertex*> v = vertices;
        const int32_t ALIGNMENT = 4;
        for (int32_t i = 0; i < nLevels; ++i)
        {
            // Update _f_ and _v_ for next level of subdivision
            std::vector<SDFace*> newFaces;
            std::vector<SDVertex*> newVertices;

            for (SDVertex* vertex : v)
            {
                vertex->child = reinterpret_cast<SDVertex*>(g_pMemoryManager->Allocate(sizeof(SDVertex), ALIGNMENT));
                vertex->child->regular = vertex->regular;
                vertex->child->boundary = vertex->boundary;
                newVertices.push_back(vertex->child);
            }            
            for (SDFace* face : f)
            {
                for (int32_t k = 0; i < 4; ++k)
                {
                    face->children[k] = reinterpret_cast<SDFace*>(g_pMemoryManager->Allocate(sizeof(SDFace), ALIGNMENT));
                    newFaces.push_back(face->children[k]);
                }
            }

            // Update vertex positions and create new edge vertices
            
            // Update vertex poisitions for even vertices
            for (SDVertex* vertex : v)
            {
                if (!vertex ->boundary)
                {
                    // Apply one-ring rule for even vertex
                    if (vertex->regular)
                        vertex->child->p = WeightOneRing(vertex, 1.f / 16.f);
                    else
                        vertex->child->p = WeightOneRing(vertex, Beta(vertex->Valence()));
                }
                else
                {
                    // Apply boundary rule for even vertex 
                    vertex->child->p = WeightBoundary(vertex, 1.f / 8.f);
                }
            }

            // Compute new odd edge vertices
            std::map<SDEdge, SDVertex*> edgeVerts;
            for (SDFace* face : f)
            {
                for (int32_t k = 0; k < 3; ++k)
                {
                    // Compute odd vertex on _k_th edge
                    SDEdge edge(face->v[k], face->v[NEXT(k)]);
                    SDVertex* vert = edgeVerts[edge];
                    if (!vert)
                    {
                        // Create and initialize new odd vertex
                        vert = reinterpret_cast<SDVertex*>(g_pMemoryManager->Allocate(sizeof(SDVertex), ALIGNMENT));
                        newVertices.push_back(vert);
                        vert->regular = true;
                        vert->boundary = (face->f[k] == nullptr);
                        vert->startFace = face->children[3];

                        // Apply edge rules to compute new vertex position
                        if (vert->boundary)
                        {
                            vert->p = 0.5f * edge.v[0]->p;
                            vert->p += 0.5f * edge.v[1]->p;
                        }
                        else 
                        {
                            vert->p = 3.f / 8.f * edge.v[0]->p;
                            vert->p += 3.f / 8.f * edge.v[1]->p;
                            vert->p += 1.f / 8.f * face->OtherVert(edge.v[0], edge.v[1])->p;
                            vert->p += 1.f / 8.f * face->f[k]->OtherVert(edge.v[0], edge.v[1])->p;
                        }
                        edgeVerts[edge] = vert;
                    }
                }
            }

            // Update new mesh topology

            // Update even vertex face pointers
            for (SDVertex* vertex : v)
            {
                int32_t vertNum = vertex->startFace->Vnum(vertex);
                vertex->child->startFace = vertex->startFace->children[vertNum];
            }

            // Update face neighbor pointers
            for (SDFace* face : f)
            {
                for (int32_t j = 0; j < 3; ++j)
                {
                    // Update children _f_ pointers for siblings
                    face->children[3]->f[j] = face->children[NEXT(j)];
                    face->children[j]->f[NEXT(j)] = face->children[3];

                    // Update children _f_ pointers for neighbor children
                    SDFace* f2 = face->f[j];
                    face->children[j]->f[j] = f2? f2->children[f2->Vnum(face->v[j])] : nullptr;
                    f2 = face->f[PREV(j)];
                    face->children[j]->f[PREV(j)] = f2? f2->children[f2->Vnum(face->v[j])] : nullptr;
                }
            }

            // Update face vertex pointers
            for (SDFace* face : f)
            {
                for (int32_t j = 0; j < 3; ++j)
                {
                    // Update child vertex pointer to new even vertex
                    face->children[j]->v[j] = face->v[j]->child;

                    // Update child vertex pointer to new odd vertex
                    SDVertex* vert = edgeVerts[SDEdge(face->v[j], face->v[NEXT(j)])];
                    face->children[j]->v[NEXT(j)] = vert;
                    face->children[NEXT(j)]->v[j] = vert;
                    face->children[3]->v[j] = vert;
                }
            }

            // Prepare for next level of subdivision
            f = newFaces;
            v = newVertices;
        }

        // Push vertices to limit surface
        std::unique_ptr<Vector3Df[]> pLimit(new Vector3Df[v.size()]);
        for (size_t i = 0; i < v.size(); ++i)
        {
            if (v[i]->boundary)
                pLimit[i] = WeightBoundary(v[i], 1.f / 5.f);
            else 
                pLimit[i] = WeightOneRing(v[i], LoopGamma(v[i]->Valence()));
        }
        for (size_t i = 0; i < v.size(); ++i) v[i]->p = pLimit[i];

        // Compute vertex tangents on limit surface
        std::vector<Vector3Df> Ns;
        Ns.reserve(v.size());
        std::vector<Vector3Df> pRing(16, Vector3Df());
        for (SDVertex* vertex : v)
        {
            Vector3Df S(0.f), T(0.f);
            int32_t valence = vertex->Valence();
            if (valence > (int32_t)pRing.size()) pRing.resize(valence);
            vertex->OneRing(&pRing[0]);
            if (!vertex->boundary)
            {
                // Compute tangents of interior face
                for (int32_t j = 0; j < valence; ++j)
                {
                    S += std::cosf(2.f * PI * j / valence) * Vector3Df(pRing[j]);
                    T += std::sinf(2.f * PI * j / valence) * Vector3Df(pRing[j]);
                }
            }
            else 
            {
                // Compute tangents of boundary face
                S = pRing[valence - 1] - pRing[0];
                if (valence == 2)
                    T = Vector3Df(pRing[0] + pRing[1] - 2 * vertex->p);
                else if (valence == 3)
                    T = pRing[1] - vertex->p;
                else if (valence == 4) // regular
                    T = Vector3Df(-1 * pRing[0] + 2 * pRing[1] + 2 * pRing[2] + -1 * pRing[3] + -2 * vertex->p);
                else 
                {
                    float theta = PI / float (valence - 1);
                    T = Vector3Df (std::sinf(theta) * (pRing[0] + pRing[valence - 1]));
                    for (int32_t k = 1; k < valence - 1; ++k)
                    {
                        float wt = (2 * std::cosf(theta) - 2) * std::sinf((k) * theta);
                        T += Vector3Df(wt * pRing[k]);
                    }
                    T = -T;
                }
            }
            Ns.push_back(CrossProduct(S, T));
        }

        // Create Triangle mesh from subdivision mesh
        {
            size_t ntris = f.size();
            std::unique_ptr<int32_t[]> verts(new int32_t[3 * ntris]);
            int32_t* vp = verts.get();
            size_t totVerts = v.size();
            std::map<SDVertex*, int> usedVerts;
            for (size_t i = 0; i < totVerts; ++i) usedVerts[v[i]] = i;
            for (size_t i = 0; i < ntris; ++i)
            {
                for (int32_t j = 0; j < 3; ++j)
                {
                    *vp = usedVerts[f[i]->v[j]];
                    ++vp;
                }
            }

            return CreateTriangleMesh(objectToWorld, worldToObject, ntris,
                verts.get(), totVerts, pLimit.get(), nullptr, &Ns[0], nullptr, 
				std::shared_ptr<SceneObjectTexture>(nullptr), std::shared_ptr<SceneObjectTexture>(nullptr),
				nullptr);
        }
    }

    std::vector<std::shared_ptr<SceneObjectShape>> CreateLoopDubDiv(
		const std::shared_ptr<Matrix4f>& objectToWorld,
		const std::shared_ptr<Matrix4f>& worldToObject,
        const ParamSet& params)
    {
        int32_t nLevels = params.FindOneInt("levels", params.FindOneInt("nLevels", 3));
        int32_t nps, nIndices;
        const int32_t* vertexIndices = params.FindInt("indices", nIndices);
        const Vector3Df* P = params.FindPoint3Df("P", nps);
        if (!vertexIndices)
        {
            return std::vector<std::shared_ptr<SceneObjectShape>>();
        }
        if (!P)
        {
            return std::vector<std::shared_ptr<SceneObjectShape>>();
        }

        // Don't actually use this for now...
        std::string scheme = params.FindOneString("scheme", "loop");
        return LoopSubDivide(objectToWorld, worldToObject, nLevels, nIndices,
                             vertexIndices, nps, P);
    }
}