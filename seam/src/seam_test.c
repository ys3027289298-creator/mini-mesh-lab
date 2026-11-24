/* seam self-test: geometry construction, weld window, crease guard, channels. */
#include "seam.h"
#include <stdio.h>
#include <math.h>

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Mesh m;
    mesh_build_folded_strip(&m, 2, 1.0f);

    PairSets ps;
    mesh_find_pairs(&m, &ps);
    check("seam_pairs_present", ps.nseam > 0);
    check("crease_pairs_present", ps.ncrease > 0);
    check("pairs_classified", ps.nseam == 2 * (SEAM_ROW + 1)
                              && ps.ncrease == (SEAM_PANELS - 1) * (SEAM_ROW + 1));

    float seam_gap = ps.seam[0].gap;
    float crease_gap = ps.crease[0].gap;
    check("gap_ordering", seam_gap < crease_gap);
    /* the seam twin is displaced tangentially on the surface, so its 3D
     * distance is slightly larger than the pure x shift F_SEAM*CELL */
    check("seam_gap_value", fabsf(seam_gap - SEAM_F_SEAM * SEAM_CELL) < 2.5e-3f);
    check("crease_gap_value", fabsf(crease_gap - SEAM_F_CREASE * SEAM_CELL) < 1e-4f);

    float L = mesh_characteristic_length(&m);
    check("characteristic_length", L > 0.0f);

    /* below the seam gap: nothing merges */
    UF u1;
    weld_position(&m, seam_gap * 0.5f, &u1);
    check("tight_eps_seams_open", count_unmerged_seams(&ps, &u1) == ps.nseam);
    check("tight_eps_creases_safe", count_mismerged_creases(&ps, &u1) == 0);

    /* inside the window: all seams merged, every crease still separate */
    UF u2;
    weld_position(&m, 0.5f * (seam_gap + crease_gap), &u2);
    check("window_seams_closed", count_unmerged_seams(&ps, &u2) == 0);
    check("window_creases_safe", count_mismerged_creases(&ps, &u2) == 0);

    /* above the crease gap: creases start collapsing */
    UF u3;
    weld_position(&m, crease_gap * 1.2f, &u3);
    check("loose_eps_creases_hit", count_mismerged_creases(&ps, &u3) > 0);

    /* guarded weld refuses the fold even above the crease gap */
    UF u4;
    weld_guarded(&m, crease_gap * 1.2f, 30.0f, &u4);
    check("guard_saves_creases", count_mismerged_creases(&ps, &u4) == 0);
    check("guard_still_joins_seams", count_unmerged_seams(&ps, &u4) == 0);

    /* UV topology: four panels + one extra island per internal seam */
    int channels = count_uv_channels(&m);
    check("uv_channels", channels == SEAM_PANELS + 2);

    /* shared group normals must be closer to ground truth than split records */
    V3 split[MAXV], shared[MAXV];
    mesh_ground_truth_normals(&m, split);
    /* split field = per-record normals, unchanged */
    for (int i = 0; i < m.nv; i++) split[i] = m.nrm[i];
    share_group_normals(&m, &u2, shared);
    float err_split = mean_normal_error_deg(&m, split);
    float err_shared = mean_normal_error_deg(&m, shared);
    check("normal_error_reduced", err_shared < err_split);

    int bright = render_bright_pixels(&m, shared, 0.001f, NULL);
    check("render_runs", bright >= 0);
    int bright_split = render_bright_pixels(&m, split, 0.001f, NULL);
    check("seams_visible_before_merge", bright_split > bright);

    /* empty / degenerate input safety */
    Mesh e; mesh_init(&e);
    check("empty_length", mesh_characteristic_length(&e) < 0.0f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
