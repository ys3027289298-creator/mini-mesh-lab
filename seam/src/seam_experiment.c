/* seam experiment: tolerance curves, three seam densities, channel diagnostics.
 *
 * Outputs CSV/PPM under the directory given in argv[1] (default artifacts/seam). */
#include "seam.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define N_DENSITY 3
#define N_SCALE   3
#define N_SWEEP   29
#define ANGLE_GUARD_DEG 30.0f
#define SHADE_TOL 0.001f

static const int seam_counts[N_DENSITY] = { 2, 7, 14 };
static const char* density_tag[N_DENSITY] = { "sparse", "medium", "dense" };
static const float scales[N_SCALE] = { 0.5f, 1.0f, 2.0f };

typedef struct {
    float L, k_seam, k_crease;
    int   nseam_pairs, ncrease_pairs, uv_channels;
} ModelInfo;

static void path_join(char* out, size_t n, const char* dir, const char* name) {
    snprintf(out, n, "%s/%s", dir, name);
}

int main(int argc, char** argv) {
    const char* dir = (argc > 1) ? argv[1] : "artifacts/seam";
    char path[512];

    FILE* fsweep = fopen((path_join(path, sizeof(path), dir, "sweep.csv"), path), "w");
    if (!fsweep) { perror(path); return 1; }
    fprintf(fsweep, "density,scale,k,mode,unmerged_seam_pairs,mismerged_crease_pairs,"
                    "bright_pixels,mean_normal_error_deg\n");

    FILE* fsumm = fopen((path_join(path, sizeof(path), dir, "summary.csv"), path), "w");
    if (!fsumm) { perror(path); return 1; }
    fprintf(fsumm, "density,scale,L,k_seam,k_crease,k_star,seam_pairs,crease_pairs,"
                   "uv_channels,variant,unmerged_vertices,mismerged_crease_pairs,"
                   "bright_pixels,mean_normal_error_deg,normal_values,position_open_pairs\n");

    ModelInfo info[N_DENSITY][N_SCALE];

    for (int d = 0; d < N_DENSITY; d++) {
        for (int s = 0; s < N_SCALE; s++) {
            Mesh m;
            mesh_build_folded_strip(&m, seam_counts[d], scales[s]);
            PairSets ps;
            mesh_find_pairs(&m, &ps);
            float L = mesh_characteristic_length(&m);

            /* empirical thresholds: largest k needed to close every seam,
             * smallest k that breaks the first crease */
            float seam_min = 0.0f, crease_first = 1.0f;
            for (int i = 0; i < ps.nseam; i++) {
                float k = ps.seam[i].gap / L;
                if (k > seam_min) seam_min = k;
            }
            for (int i = 0; i < ps.ncrease; i++) {
                float k = ps.crease[i].gap / L;
                if (k < crease_first) crease_first = k;
            }
            float k_star = sqrtf(seam_min * crease_first);

            info[d][s].L = L;
            info[d][s].k_seam = seam_min;
            info[d][s].k_crease = crease_first;
            info[d][s].nseam_pairs = ps.nseam;
            info[d][s].ncrease_pairs = ps.ncrease;
            info[d][s].uv_channels = count_uv_channels(&m);

            /* log-spaced k sweep at unit scale, both weld modes */
            if (s == 1) {
                for (int w = 0; w < N_SWEEP; w++) {
                    float lgk = -5.0f + (float)w * (5.0f / (N_SWEEP - 1));
                    float k = powf(10.0f, lgk);
                    float eps = k * L;
                    V3 nrm[MAXV];
                    UF u;
                    weld_position(&m, eps, &u);
                    share_group_normals(&m, &u, nrm);
                    int bright = render_bright_pixels(&m, nrm, SHADE_TOL, NULL);
                    fprintf(fsweep, "%s,%.2f,%.6g,position,%d,%d,%d,%.4f\n",
                            density_tag[d], scales[s], k,
                            count_unmerged_seams(&ps, &u),
                            count_mismerged_creases(&ps, &u), bright,
                            mean_normal_error_deg(&m, nrm));
                    UF ug;
                    weld_guarded(&m, eps, ANGLE_GUARD_DEG, &ug);
                    V3 gnrm[MAXV];
                    share_group_normals(&m, &ug, gnrm);
                    int gbright = render_bright_pixels(&m, gnrm, SHADE_TOL, NULL);
                    fprintf(fsweep, "%s,%.2f,%.6g,guarded,%d,%d,%d,%.4f\n",
                            density_tag[d], scales[s], k,
                            count_unmerged_seams(&ps, &ug),
                            count_mismerged_creases(&ps, &ug), gbright,
                            mean_normal_error_deg(&m, gnrm));
                }
            }

            /* fixed tolerance at the geometric-mean operating point */
            float eps = k_star * L;
            UF up; weld_position(&m, eps, &up);
            UF ug; weld_guarded(&m, eps, ANGLE_GUARD_DEG, &ug);

            V3 per_record[MAXV], shared_pos[MAXV], shared_grd[MAXV];
            for (int i = 0; i < m.nv; i++) per_record[i] = m.nrm[i];
            share_group_normals(&m, &up, shared_pos);
            share_group_normals(&m, &ug, shared_grd);

            struct { const char* tag; UF* u; V3* nrm; } cases[3] = {
                { "position_match_only", &up, per_record },
                { "channel_open_position", &up, shared_pos },
                { "channel_open_guarded", &ug, shared_grd },
            };
            for (int c = 0; c < 3; c++) {
                int bright = render_bright_pixels(&m, cases[c].nrm, SHADE_TOL, NULL);
                int open_normal_pairs = count_open_normal_seams(&ps, cases[c].nrm);
                int normal_values = count_normal_values(&m, cases[c].u, cases[c].nrm);
                int position_open = count_unmerged_seams(&ps, cases[c].u);
                int mismerged = count_mismerged_creases(&ps, cases[c].u);
                float nerr = mean_normal_error_deg(&m, cases[c].nrm);
                fprintf(fsumm,
                        "%s,%.2f,%.6g,%.6g,%.6g,%.6g,%d,%d,%d,%s,%d,%d,%d,%.4f,%d,%d\n",
                        density_tag[d], scales[s], L, seam_min, crease_first, k_star,
                        ps.nseam, ps.ncrease, info[d][s].uv_channels,
                        cases[c].tag,
                        2 * open_normal_pairs,   /* unmerged vertices */
                        mismerged, bright, nerr, normal_values, position_open);
            }

            /* unit-scale renders for the visual record */
            if (s == 1) {
                char ppm[512];
                snprintf(ppm, sizeof(ppm), "%s/render_%s_split.ppm", dir, density_tag[d]);
                render_bright_pixels(&m, per_record, SHADE_TOL, ppm);
                snprintf(ppm, sizeof(ppm), "%s/render_%s_merged.ppm", dir, density_tag[d]);
                render_bright_pixels(&m, shared_pos, SHADE_TOL, ppm);
            }
        }
    }
    fclose(fsweep);
    fclose(fsumm);

    /* scale invariance check: k thresholds must be identical across scales */
    FILE* fscale = fopen((path_join(path, sizeof(path), dir, "scale.csv"), path), "w");
    fprintf(fscale, "density,scale,L,k_seam,k_crease\n");
    for (int d = 0; d < N_DENSITY; d++)
        for (int s = 0; s < N_SCALE; s++)
            fprintf(fscale, "%s,%.2f,%.6g,%.6g,%.6g\n", density_tag[d], scales[s],
                    info[d][s].L, info[d][s].k_seam, info[d][s].k_crease);
    fclose(fscale);

    printf("done: wrote sweep.csv summary.csv scale.csv and renders to %s\n", dir);
    return 0;
}
