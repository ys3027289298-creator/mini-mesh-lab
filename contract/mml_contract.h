/*
 * mini-mesh-lab shared mesh contract.
 *
 * This header is the boundary contract between welded/, normal/ and bbox/.
 * It is intentionally storage-free: descriptors borrow caller-owned arrays.
 * See contract/CONTRACT.md for the ownership, indexing, invalid-value and
 * evolution rules.
 */
#ifndef MML_CONTRACT_H
#define MML_CONTRACT_H

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#define MML_ABI_VERSION 1u

/* Valid indices are [0, MML_INVALID_INDEX). */
#define MML_INVALID_INDEX UINT32_MAX

typedef uint32_t mml_index;

typedef struct mml_vec3 {
    float x, y, z;
} mml_vec3;

typedef struct mml_triangle {
    mml_index a, b, c;
} mml_triangle;

/* Every cross-module struct starts with this header. */
typedef struct mml_header {
    uint32_t size;        /* caller-visible sizeof(struct), in bytes */
    uint32_t abi_version; /* major/generation version */
} mml_header;

/* Initialize a current-version stack/heap struct and zero its fields. */
#define MML_STRUCT_INIT(ptr, Type)                                           \
    do {                                                                     \
        memset((ptr), 0, sizeof(Type));                                      \
        (ptr)->header.size = (uint32_t)sizeof(Type);                         \
        (ptr)->header.abi_version = MML_ABI_VERSION;                          \
    } while (0)

/* Byte offset just past Field; used to test whether a caller knew a field. */
#define MML_FIELD_END(Type, Field)                                           \
    ((uint32_t)(offsetof(Type, Field) + sizeof(((Type *)0)->Field)))

#define MML_HAS_FIELD(header_ptr, Type, Field)                               \
    ((header_ptr) != NULL &&                                                 \
     (header_ptr)->size >= MML_FIELD_END(Type, Field))

typedef enum mml_status {
    MML_OK = 0,
    MML_E_ARGUMENT = -1,
    MML_E_NULL = -2,
    MML_E_HEADER_SIZE = -3,
    MML_E_ABI_VERSION = -4,
    MML_E_CAPACITY = -5,
    MML_E_INDEX = -6,
    MML_E_NONFINITE = -7,
    MML_E_EMPTY = -8,
    MML_E_DEGENERATE = -9,
    MML_E_NOMEM = -10,
    MML_E_OVERFLOW = -11
} mml_status;

static inline const char *mml_status_str(mml_status s)
{
    switch (s) {
    case MML_OK: return "MML_OK";
    case MML_E_ARGUMENT: return "MML_E_ARGUMENT";
    case MML_E_NULL: return "MML_E_NULL";
    case MML_E_HEADER_SIZE: return "MML_E_HEADER_SIZE";
    case MML_E_ABI_VERSION: return "MML_E_ABI_VERSION";
    case MML_E_CAPACITY: return "MML_E_CAPACITY";
    case MML_E_INDEX: return "MML_E_INDEX";
    case MML_E_NONFINITE: return "MML_E_NONFINITE";
    case MML_E_EMPTY: return "MML_E_EMPTY";
    case MML_E_DEGENERATE: return "MML_E_DEGENERATE";
    case MML_E_NOMEM: return "MML_E_NOMEM";
    case MML_E_OVERFLOW: return "MML_E_OVERFLOW";
    }
    return "MML_E_UNKNOWN";
}

static inline mml_status mml_check_header(const mml_header *h,
                                          uint32_t minimum_size)
{
    if (h == NULL)
        return MML_E_NULL;
    if (h->size < minimum_size)
        return MML_E_HEADER_SIZE;
    if (h->abi_version != MML_ABI_VERSION)
        return MML_E_ABI_VERSION;
    return MML_OK;
}

static inline int mml_vec3_is_finite(mml_vec3 v)
{
    return isfinite(v.x) && isfinite(v.y) && isfinite(v.z);
}

enum mml_mesh_flags {
    /* When set, mml_mesh_view::vertex_normals points at vertex_count vec3s. */
    MML_MESH_HAS_NORMALS = 1u << 0
};

/*
 * Read-only mesh view. All pointers are borrowed only for the duration of
 * the callee call. Indices are local to this descriptor.
 *
 * V1's mandatory prefix ends at triangles. Future optional fields must be
 * appended after this field and gated with MML_HAS_FIELD().
 */
typedef struct mml_mesh_view {
    mml_header header;
    uint32_t flags;
    uint32_t vertex_count;
    uint32_t triangle_count;
    const mml_vec3 *positions;
    const mml_vec3 *vertex_normals; /* optional; invalid normal = (0,0,0) */
    const mml_triangle *triangles;
} mml_mesh_view;

#define MML_MESH_VIEW_V1_SIZE MML_FIELD_END(mml_mesh_view, triangles)

/*
 * Caller-owned output mesh storage. The descriptor does not own the arrays.
 * vertex/triangle_capacity are element capacities; counts are updated by the
 * callee.
 */
typedef struct mml_mesh_buffer {
    mml_header header;
    uint32_t flags;
    uint32_t vertex_count;
    uint32_t vertex_capacity;
    uint32_t triangle_count;
    uint32_t triangle_capacity;
    mml_vec3 *positions;
    mml_vec3 *vertex_normals; /* optional derived storage; weld leaves unset */
    mml_triangle *triangles;
} mml_mesh_buffer;

#define MML_MESH_BUFFER_V1_SIZE MML_FIELD_END(mml_mesh_buffer, triangles)

typedef struct mml_box {
    mml_header header;
    mml_vec3 min; /* meaningful only when valid != 0 */
    mml_vec3 max;
    uint32_t valid;
} mml_box;

#define MML_BOX_V1_SIZE MML_FIELD_END(mml_box, valid)

/* Validate a view and all indices/values needed by the three modules. */
static inline mml_status mml_mesh_view_validate(const mml_mesh_view *m)
{
    mml_status st;
    uint32_t i;

    if (m == NULL)
        return MML_E_NULL;

    st = mml_check_header(&m->header, MML_MESH_VIEW_V1_SIZE);
    if (st != MML_OK)
        return st;

    if (m->vertex_count > 0 && m->positions == NULL)
        return MML_E_NULL;
    if (m->triangle_count > 0 && m->triangles == NULL)
        return MML_E_NULL;
    if ((m->flags & MML_MESH_HAS_NORMALS) != 0) {
        if (m->vertex_count > 0 && m->vertex_normals == NULL)
            return MML_E_NULL;
    }

    for (i = 0; i < m->vertex_count; ++i) {
        if (!mml_vec3_is_finite(m->positions[i]))
            return MML_E_NONFINITE;
        if ((m->flags & MML_MESH_HAS_NORMALS) != 0 &&
            !mml_vec3_is_finite(m->vertex_normals[i]))
            return MML_E_NONFINITE;
    }

    for (i = 0; i < m->triangle_count; ++i) {
        const mml_triangle *t = &m->triangles[i];
        if (t->a >= m->vertex_count || t->b >= m->vertex_count ||
            t->c >= m->vertex_count)
            return MML_E_INDEX;
    }

    return MML_OK;
}

static inline mml_status mml_box_init(mml_box *box)
{
    if (box == NULL)
        return MML_E_NULL;
    MML_STRUCT_INIT(box, mml_box);
    box->valid = 0;
    return MML_OK;
}

#endif /* MML_CONTRACT_H */
