
/* bbox: AABB from vertices, with degenerate/sliver detection. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAXV 32

typedef struct { float x, y, z; } V3;
typedef struct { V3 min, max; int valid; } Box;

static void box_empty(Box* b) {
    b->min.x = b->min.y = b->min.z = 1e30f;
    b->max.x = b->max.y = b->max.z = -1e30f;
    b->valid = 0;
}

static int box_grow(Box* b, V3 p) {
    if (b->min.x > p.x) b->min.x = p.x;
    if (b->min.y > p.y) b->min.y = p.y;
    if (b->min.z > p.z) b->min.z = p.z;
    if (b->max.x < p.x) b->max.x = p.x;
    if (b->max.y < p.y) b->max.y = p.y;
    if (b->max.z < p.z) b->max.z = p.z;
    b->valid = 1;
    return 0;
}

static V3 box_size(const Box* b) {
    V3 s = { b->max.x - b->min.x, b->max.y - b->min.y, b->max.z - b->min.z };
    return s;
}

static float box_diag(const Box* b) {
    V3 s = box_size(b);
    return sqrtf(s.x*s.x + s.y*s.y + s.z*s.z);
}

/* sliver: smallest extent less than ratio of diagonal */
static int box_is_sliver(const Box* b, float ratio) {
    if (!b->valid) return -1;
    if (ratio <= 0.0f || ratio >= 1.0f) return -1;
    V3 s = box_size(b);
    float m = s.x < s.y ? s.x : s.y;
    if (s.z < m) m = s.z;
    float d = box_diag(b);
    if (d <= 0.0f) return -1;
    return (m < ratio * d) ? 1 : 0;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Box b; box_empty(&b);
    check("empty_invalid", b.valid == 0);
    check("empty_sliver_err", box_is_sliver(&b, 0.1f) == -1);

    box_grow(&b, (V3){0,0,0});
    check("single_valid", b.valid == 1);
    check("single_zero_diag", fabsf(box_diag(&b)) < 1e-6f);

    box_grow(&b, (V3){1,1,1});
    V3 s = box_size(&b);
    check("cube_size", fabsf(s.x-1)<1e-6f && fabsf(s.y-1)<1e-6f && fabsf(s.z-1)<1e-6f);
    check("cube_not_sliver", box_is_sliver(&b, 0.1f) == 0);

    Box sl; box_empty(&sl);
    box_grow(&sl, (V3){0,0,0});
    box_grow(&sl, (V3){10,10,0.01f});
    check("sliver_detected", box_is_sliver(&sl, 0.1f) == 1);

    check("bad_ratio", box_is_sliver(&sl, 0.0f) == -1 && box_is_sliver(&sl, 1.0f) == -1);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
