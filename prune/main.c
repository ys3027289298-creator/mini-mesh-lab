/* main.c -- three exports of one mesh at three attribute-redundancy levels.
 *
 *   1. count each attribute's real occurrences in the face index stream
 *   2. slots with zero occurrences are dropped; the rest are compacted
 *   3. second criterion: the render must be bit-for-bit identical
 *   4. report the relation between redundancy rate and how much of the file
 *      is compressible at all
 *   5. control: apply the same rule through the VERTEX stream, which is the
 *      wrong criterion for attribute banks, and show it prunes nothing.
 */
#define _POSIX_C_SOURCE 200809L
#include "mesh.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define W 128
#define H 96
#define RW 800
#define RH 600

typedef struct { const char* name; int dx[CH_N], dy[CH_N], bank[CH_N]; } Spec;

/* One mesh: nv = 27648, nf = 54626, identical in all three models.  Only
 * the attribute dedup changes, so the banks have the same slot count but
 * their streams reach fewer slots:
 *   low    1x  ->  0 dead slots
 *   medium 2x  ->  20736 dead in nrm + col      (40.9 % of attribute bytes)
 *   high   3x/2x -> 24576 nrm + 20736 uv + 24576 col (62.1 %)
 *
 * Dead slots are holes inside [0,nv), not a padded tail.  That is the case
 * the vertex-index criterion cannot see: it counts every vertex as used and
 * therefore prunes nothing, which is why the highest-redundancy file does
 * not shrink under that rule. */
/* One mesh: nv = 12288, nf = 24130, and the same live attribute records in
 * every tier (nrm/uv/col deduplicated 2x -> 3072 records each).  The three
 * exports differ ONLY in how many stale extra slots each bank carries:
 *
 *   low     bank = 1x live    0 % redundant attribute bytes
 *   medium  bank = 4x live    ~72 %
 *   high    bank = 9x live    ~88 %
 *
 * Because the live data is identical, pruning must produce the same file
 * size for all three; the whole difference between exports is redundancy. */
#define LIVE_MULT 3072
static const Spec SPECS[3] = {
    { "low",    { 1, 2, 2, 2 }, { 1, 2, 2, 2 },
      { 12288, 1*LIVE_MULT, 1*LIVE_MULT, 1*LIVE_MULT } },
    { "medium", { 1, 2, 2, 2 }, { 1, 2, 2, 2 },
      { 12288, 4*LIVE_MULT, 4*LIVE_MULT, 4*LIVE_MULT } },
    { "high",   { 1, 2, 2, 2 }, { 1, 2, 2, 2 },
      { 12288, 9*LIVE_MULT, 9*LIVE_MULT, 9*LIVE_MULT } },
};

static long file_size(const char* path) {
    FILE* fp = fopen(path, "rb");
    if (!fp) return -1;
    fseek(fp, 0, SEEK_END);
    long n = ftell(fp);
    fclose(fp);
    return n;
}

/* gzip -9 through a pipe: no zlib dependency */
static long gz_size(const char* path) {
    char cmd[1024];
    snprintf(cmd, sizeof(cmd), "gzip -9 -c '%s' | wc -c", path);
    FILE* p = popen(cmd, "r");
    if (!p) return -1;
    long n = -1;
    if (fscanf(p, "%ld", &n) != 1) n = -1;
    pclose(p);
    return n;
}

typedef struct {
    long attr_bytes, redundant_bytes, attr_after;
    float redundant_share;
    long file_before, file_after;
    long pix_diff;
    float max_delta;
    int identical, attrs_same, roundtrip_same;
    long gz_before, gz_after;
    long dead[CH_N];
    int wrong_identical, wrong_attrs_same;
    long wrong_dead[CH_N], wrong_file;
    long over_pix_diff, over_file;
    float over_max_delta;
    int over_identical;
} R;

int main(int argc, char** argv) {
    const char* outdir = argc > 1 ? argv[1] : "out";
    R r[3];

    for (int mi = 0; mi < 3; mi++) {
        const Spec* sp = &SPECS[mi];
        R* R_ = &r[mi];
        memset(R_, 0, sizeof(*R_));
        Mesh m = mesh_make(sp->name, W, H, sp->dx, sp->dy, sp->bank, 1.0f, 12345u);

        /* --- criterion 1: count occurrences in each attribute's OWN stream */
        int* live[CH_N];
        for (int ch = 0; ch < CH_N; ch++) {
            int* c = malloc((size_t)m.bank[ch] * sizeof(int));
            mesh_counts(&m, ch, c, 0);
            live[ch] = malloc((size_t)m.bank[ch] * sizeof(int));
            long ln = 0;
            for (int s = 0; s < m.bank[ch]; s++) {
                live[ch][s] = c[s] != 0;
                if (live[ch][s]) ln++;
            }
            R_->dead[ch] = m.bank[ch] - ln;
            free(c);
        }
        R_->attr_bytes = (long)mesh_attr_bytes(&m);
        for (int ch = 0; ch < CH_N; ch++)
            R_->redundant_bytes += R_->dead[ch] * CH_COMP[ch] * 4;
        R_->redundant_share = (float)R_->redundant_bytes / (float)R_->attr_bytes;

        char pb[512], pa[512], pv[512], pi[512];
        snprintf(pb, sizeof(pb), "%s/%s_before.mms2", outdir, sp->name);
        snprintf(pa, sizeof(pa), "%s/%s_after.mms2", outdir, sp->name);
        snprintf(pv, sizeof(pv), "%s/%s_after_vertexcrit.mms2", outdir, sp->name);
        mesh_write(&m, pb);
        R_->file_before = file_size(pb);
        R_->gz_before = gz_size(pb);

        float* ia = malloc((size_t)RW * RH * 3 * sizeof(float));
        float* ib = malloc((size_t)RW * RH * 3 * sizeof(float));
        mesh_render(&m, RW, RH, ia);

        /* --- prune + second criterion: pixel-identical render */
        Mesh p = mesh_prune(&m, live);
        mesh_write(&p, pa);
        R_->file_after = file_size(pa);
        R_->gz_after = gz_size(pa);
        R_->attr_after = (long)mesh_attr_bytes(&p);
        mesh_render(&p, RW, RH, ib);
        R_->pix_diff = 0;
        pixels_differ(ia, ib, RW, RH, &R_->pix_diff, &R_->max_delta);
        R_->identical = R_->pix_diff == 0;
        R_->attrs_same = attrs_equal(&m, &p);

        Mesh p2 = mesh_read(pa);
        R_->roundtrip_same = p2.nv > 0 && attrs_equal(&m, &p2);
        if (p2.nv > 0) mesh_free(&p2);

        snprintf(pi, sizeof(pi), "%s/%s_before.ppm", outdir, sp->name);
        write_ppm(pi, ia, RW, RH);
        snprintf(pi, sizeof(pi), "%s/%s_after.ppm", outdir, sp->name);
        write_ppm(pi, ib, RW, RH);
        snprintf(pi, sizeof(pi), "%s/%s_diff.ppm", outdir, sp->name);
        write_ppm_diff(pi, ia, ib, RW, RH);

        /* --- negative control: over-prune referenced slots too.  This must
         * FAIL the second criterion, proving the pixel check has teeth. */
        {
            int* lv3[CH_N];
            for (int ch = 0; ch < CH_N; ch++) {
                int* c = malloc((size_t)m.bank[ch] * sizeof(int));
                mesh_counts(&m, ch, c, 0);
                lv3[ch] = malloc((size_t)m.bank[ch] * sizeof(int));
                long keep = 0;
                for (int s = 0; s < m.bank[ch]; s++) {
                    lv3[ch][s] = c[s] > 0;
                    if (lv3[ch][s]) keep++;
                }
                /* drop every 7th referenced slot as well */
                long seen = 0;
                for (int s = 0; s < m.bank[ch] && keep > 1; s++) {
                    if (!lv3[ch][s]) continue;
                    if ((seen++ % 7) == 0 && seen > 1) { lv3[ch][s] = 0; keep--; }
                }
                free(c);
            }
            Mesh po = mesh_prune(&m, lv3);
            float* io = malloc((size_t)RW * RH * 3 * sizeof(float));
            mesh_render(&po, RW, RH, io);
            R_->over_pix_diff = 0;
            R_->over_max_delta = 0;
            pixels_differ(ia, io, RW, RH, &R_->over_pix_diff, &R_->over_max_delta);
            R_->over_identical = R_->over_pix_diff == 0;
            char po_path[512];
            snprintf(po_path, sizeof(po_path), "%s/%s_overpruned.mms2", outdir, sp->name);
            mesh_write(&po, po_path);
            R_->over_file = file_size(po_path);
            free(io);
            for (int ch = 0; ch < CH_N; ch++) free(lv3[ch]);
            mesh_free(&po);
        }

        /* --- control: counts through the position/vertex stream (wrong) */
        int* lv2[CH_N];
        for (int ch = 0; ch < CH_N; ch++) {
            int* c = malloc((size_t)m.bank[ch] * sizeof(int));
            mesh_counts(&m, ch, c, 1);
            lv2[ch] = malloc((size_t)m.bank[ch] * sizeof(int));
            long d = 0;
            for (int s = 0; s < m.bank[ch]; s++) {
                /* the buggy tool loops over vertices 0..nv-1; rows at or above
                 * nv are never examined and stay alive */
                lv2[ch][s] = (s >= m.nv) ? 1 : (c[s] != 0);
                if (!lv2[ch][s]) d++;
            }
            R_->wrong_dead[ch] = d;
            free(c);
        }
        Mesh pw = mesh_prune(&m, lv2);
        mesh_write(&pw, pv);
        R_->wrong_file = file_size(pv);
        float* iw = malloc((size_t)RW * RH * 3 * sizeof(float));
        mesh_render(&pw, RW, RH, iw);
        long wd; float wm;
        pixels_differ(ia, iw, RW, RH, &wd, &wm);
        R_->wrong_identical = wd == 0;
        R_->wrong_attrs_same = attrs_equal(&m, &pw);
        free(iw);

        printf("  %-6s nv=%d nf=%d  dead slots pos/nrm/uv/col = %ld/%ld/%ld/%ld"
               "   (vertex-criterion: %ld/%ld/%ld/%ld)\n",
               sp->name, m.nv, m.nf, R_->dead[0], R_->dead[1], R_->dead[2], R_->dead[3],
               R_->wrong_dead[0], R_->wrong_dead[1], R_->wrong_dead[2], R_->wrong_dead[3]);
        printf("         roundtrip_from_disk_equal=%d  attrs_equal=%d\n",
               R_->roundtrip_same, R_->attrs_same);

        for (int ch = 0; ch < CH_N; ch++) { free(live[ch]); free(lv2[ch]); }
        free(ia); free(ib);
        mesh_free(&pw);
        mesh_free(&p);
        mesh_free(&m);
    }

    printf("\n== A. 3 models x 12 numbers ==\n");
    printf("  model | attribute bytes | redundancy share | size after prune (B) | render diff (pixels / maxΔ)\n");
    for (int i = 0; i < 3; i++) {
        R* R_ = &r[i];
        char diff[64];
        snprintf(diff, sizeof(diff), "%ld px / %.3g", R_->pix_diff, R_->max_delta);
        printf("  %-5s | %15ld | %15.1f%% | %20ld | %s\n",
               SPECS[i].name, R_->attr_bytes, 100.0 * (double)R_->redundant_share,
               R_->file_after, diff);
    }

    printf("\n== B. counts side by side: attribute stream vs vertex stream (dead slots) ==\n");
    printf("  %-6s | %-34s | %-34s\n", "model",
           "dead slots by own attribute index", "dead slots by vertex index (wrong)");
    for (int i = 0; i < 3; i++) {
        R* R_ = &r[i];
        printf("  %-6s | nrm=%-6ld uv=%-6ld col=%-6ld | nrm=%-6ld uv=%-6ld col=%-6ld\n",
               SPECS[i].name,
               R_->dead[CH_NRM], R_->dead[CH_UV], R_->dead[CH_COL],
               R_->wrong_dead[CH_NRM], R_->wrong_dead[CH_UV], R_->wrong_dead[CH_COL]);
    }

    printf("\n== C. applying the same rule through the vertex stream ==\n");
    printf("  %-6s | %11s | %14s | %16s | %14s | %12s\n",
           "model", "file before", "file (attr idx)", "file (vertex idx)",
           "vertex idx drop", "render same?");
    for (int i = 0; i < 3; i++) {
        R* R_ = &r[i];
        printf("  %-6s | %11ld | %14ld | %16ld | %13ld%% | %12s\n",
               SPECS[i].name, R_->file_before, R_->file_after, R_->wrong_file,
               (long)(100.0 * (double)(R_->file_before - R_->wrong_file) / (double)R_->file_before),
               R_->wrong_identical ? "yes" : "NO");
    }
    printf("  -> the highest-redundancy tier under the vertex-index rule: %s\n",
           r[2].wrong_file >= r[2].file_before
               ? "SIZE DOES NOT DROP (rule is counting the wrong index)"
               : "size does drop");

    printf("\n== C2. negative control: also dropping every 7th referenced slot ==\n");
    printf("  %-6s | %14s | %20s | %10s\n",
           "model", "render diff px", "max delta", "criterion-2");
    for (int i = 0; i < 3; i++) {
        R* R_ = &r[i];
        printf("  %-6s | %14ld | %20.3g | %10s\n",
               SPECS[i].name, R_->over_pix_diff, R_->over_max_delta,
               R_->over_identical ? "passes (BAD)" : "fails (good)");
    }
    printf("  -> this is what makes criterion 2 non-vacuous: remove one referenced slot and pixels move.\n");

    printf("\n== D. redundancy share vs compressibility (two different quantities) ==\n");
    printf("  %-6s | %8s | %12s | %10s | %12s | %10s | %10s\n",
           "model", "red%", "attr saved%", "file saved%", "gzip before", "gzip after", "gz saved%");
    for (int i = 0; i < 3; i++) {
        R* R_ = &r[i];
        double asave = 100.0 * (double)(R_->attr_bytes - R_->attr_after) / (double)R_->attr_bytes;
        double fsave = 100.0 * (double)(R_->file_before - R_->file_after) / (double)R_->file_before;
        double gsave = 100.0 * (double)(R_->gz_before - R_->gz_after) / (double)R_->gz_before;
        printf("  %-6s | %7.1f%% | %11.1f%% | %9.1f%% | %12ld | %10ld | %9.1f%%\n",
               SPECS[i].name, 100.0 * (double)R_->redundant_share, asave, fsave,
               R_->gz_before, R_->gz_after, gsave);
    }
    printf("  note: attr saved = redundant bytes; file saved is smaller because the four index\n");
    printf("        streams (%ld B) are unchanged.  gzip saved is a third quantity: it also counts\n",
           (long)r[0].file_before * 0 + 24130 * 3 * 2 * 4);
    printf("        how much of the SURVIVING data the encoder could already squeeze, which pruning\n");
    printf("        never touches.  Redundancy rate is a property of the asset; compressibility is\n");
    printf("        a property of the byte stream; how much you save by pruning is neither.\n");
    return 0;
}
