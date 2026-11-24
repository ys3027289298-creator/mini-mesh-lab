/* weld_det: deterministic vertex welding by quantized total order.
 *
 * Total order on vertices: lexicographic (qx, qy, qz, original_index),
 * where q* = llround(coord / step) are int64 quantization keys.
 *
 * Why this is a total order (and why the old eps-ball weld was not):
 *  - Quantization maps each coordinate to an integer by one IEEE-754
 *    correctly-rounded division plus llround (round-half-away-from-zero,
 *    exactly specified). On any IEEE-compliant platform the key is
 *    identical, independent of instruction selection and FP flags.
 *  - Integer lexicographic comparison is a strict total order: no NaN,
 *    no rounding ambiguity at comparison time.
 *  - The original index is unique per vertex, so the comparator never
 *    returns 0 for two distinct entries: the sorted sequence does not
 *    depend on qsort being stable or on its implementation.
 *  - Merge equivalence is "quantized keys equal", which IS an
 *    equivalence relation (reflexive/symmetric/transitive). The old
 *    "distance <= eps" relation is not transitive, so its grouping
 *    depended on scan order and on FP contraction (FMA) in dx*dx+...;
 *    different machines could legitimately produce different groups.
 *
 * Preconditions for the guarantee: coordinates are finite (no NaN/Inf),
 * step > 0 and finite, and |coord/step| < 2^62 so keys fit int64.
 *
 * Usage: weld_det <mesh.obj> <step>
 *   stdout: canonical welded mesh ("v x y z" / "f a b c", 0-based)
 *   stderr: stats nv, nu, collision_pairs, faces_in, faces_out
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXV 65536
#define MAXF 131072

typedef struct { double x, y, z; } V3;
typedef struct { int a, b, c; } Tri;

typedef struct {
    long long qx, qy, qz;
    int idx;                    /* original vertex index, unique */
} Key;

static int key_cmp(const void* pa, const void* pb) {
    const Key* a = (const Key*)pa;
    const Key* b = (const Key*)pb;
    if (a->qx != b->qx) return a->qx < b->qx ? -1 : 1;
    if (a->qy != b->qy) return a->qy < b->qy ? -1 : 1;
    if (a->qz != b->qz) return a->qz < b->qz ? -1 : 1;
    if (a->idx != b->idx) return a->idx < b->idx ? -1 : 1;
    return 0;
}

static V3 verts[MAXV];
static int nv = 0;
static Tri faces[MAXF];
static int nf = 0;

static int read_mesh(const char* path) {
    FILE* fp = fopen(path, "r");
    if (!fp) return -1;
    char line[512];
    while (fgets(line, sizeof line, fp)) {
        if (line[0] == 'v' && (line[1] == ' ' || line[1] == '\t')) {
            if (nv >= MAXV) { fclose(fp); return -2; }
            if (sscanf(line + 1, "%lf %lf %lf",
                       &verts[nv].x, &verts[nv].y, &verts[nv].z) != 3) {
                fclose(fp); return -3;
            }
            nv++;
        } else if (line[0] == 'f' && (line[1] == ' ' || line[1] == '\t')) {
            if (nf >= MAXF) { fclose(fp); return -2; }
            int a, b, c;
            if (sscanf(line + 1, "%d %d %d", &a, &b, &c) != 3) {
                fclose(fp); return -3;
            }
            faces[nf].a = a - 1; faces[nf].b = b - 1; faces[nf].c = c - 1;
            nf++;
        }
    }
    fclose(fp);
    return 0;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        fprintf(stderr, "usage: %s <mesh.obj> <step>\n", argv[0]);
        return 2;
    }
    double step = strtod(argv[2], NULL);
    if (!(step > 0.0) || !isfinite(step)) {
        fprintf(stderr, "bad step\n");
        return 2;
    }
    int rc = read_mesh(argv[1]);
    if (rc != 0) {
        fprintf(stderr, "read error %d\n", rc);
        return 2;
    }
    for (int i = 0; i < nv; i++) {
        if (!isfinite(verts[i].x) || !isfinite(verts[i].y) || !isfinite(verts[i].z)) {
            fprintf(stderr, "non-finite vertex %d\n", i);
            return 2;
        }
    }
    for (int i = 0; i < nf; i++) {
        if (faces[i].a < 0 || faces[i].a >= nv ||
            faces[i].b < 0 || faces[i].b >= nv ||
            faces[i].c < 0 || faces[i].c >= nv) {
            fprintf(stderr, "face %d out of range\n", i);
            return 2;
        }
    }

    static Key keys[MAXV];
    static int map[MAXV];       /* original index -> welded index */
    for (int i = 0; i < nv; i++) {
        keys[i].qx = llround(verts[i].x / step);
        keys[i].qy = llround(verts[i].y / step);
        keys[i].qz = llround(verts[i].z / step);
        keys[i].idx = i;
    }
    qsort(keys, nv, sizeof(Key), key_cmp);

    /* group consecutive equal quantized keys; representative is the
     * group's first entry in sorted order = smallest original index */
    int nu = 0;
    long long collision_pairs = 0;
    for (int i = 0; i < nv; ) {
        int j = i + 1;
        while (j < nv &&
               keys[j].qx == keys[i].qx &&
               keys[j].qy == keys[i].qy &&
               keys[j].qz == keys[i].qz) j++;
        long long k = j - i;
        collision_pairs += k * (k - 1) / 2;
        int rep = keys[i].idx;
        for (int t = i; t < j; t++) map[keys[t].idx] = nu;
        /* welded vertex sequence: sorted order, representative coords */
        printf("v %.17g %.17g %.17g\n", verts[rep].x, verts[rep].y, verts[rep].z);
        nu++;
        i = j;
    }

    int faces_out = 0;
    for (int i = 0; i < nf; i++) {
        int a = map[faces[i].a], b = map[faces[i].b], c = map[faces[i].c];
        if (a == b || b == c || a == c) continue;   /* degenerate after weld */
        printf("f %d %d %d\n", a, b, c);
        faces_out++;
    }

    fprintf(stderr, "stats nv=%d nu=%d collision_pairs=%lld faces_in=%d faces_out=%d\n",
            nv, nu, collision_pairs, nf, faces_out);
    return 0;
}
