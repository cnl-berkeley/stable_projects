#!/usr/bin/env python3
"""
dice/dice_prcusp_vs_other_sulci.py

Overlap (Dice similarity coefficient, DSC) between each of the eight PMC sulci
present in every hemisphere and the maximally face-selective region, and the
comparison of the posterior precuneal sulcus (prcus-p) with each other sulcus.

Reproduces
  Results, "The prcus-p and its immediate vicinity overlap with the maximally
    face-selective region in PMC": prcus-p vs each of the other seven sulci in
    the 0 mm, 0-5 mm and 0-10 mm windows (Wilcoxon signed-rank tests,
    Holm-corrected across the seven comparisons within each window and
    hemisphere); the sulcus with the highest mean DSC at every distance; the
    distance at which the prcus-p DSC peaks
  Methods, "Dice coefficient analysis...": the share of hemispheres in which
    the mcgs does not overlap the face-selective region across 0-10 mm
  Figure 5 caption: group-mean 0-mm DSC of the prcus-p
  Supplemental Material, "Dice similarity coefficient (DSC) results...": median
    area of the face-selective region and all prcus-p comparisons at the 2.5%
    and 10% thresholds

Method
  DSC = 2|S & F| / (|S| + |F|) on surface area, where F is the highest-mass
  connected component of the top-k% face-selective (faces - objects) PMC
  vertices and S is the sulcal label dilated by r mm (r = 0, 2.5, 5, 7.5, 10)
  of geodesic distance on the participant's native surface. For the 0-5 and
  0-10 mm windows each participant's DSC is the area under the DSC-by-distance
  curve (trapezoidal rule) divided by the window width.

Inputs (data/)
  dice_overlap.csv          DSC per participant, hemisphere, threshold, sulcus
                            and dilation distance (computed from FreeSurfer
                            surfaces with surface_producers/dice_dilated_native.py)
  face_patch_geometry.csv   area of the face-selective region

Outputs (outputs/dice/)
  dice_prcusp_vs_other_sulci.txt     all statistics, every threshold
  dice_prcusp_vs_other_sulci.csv     Wilcoxon and Holm p-values, every test
  dice_window_values.csv             per-participant window DSC (input to
                                     dice_group_tests.R)

Run from the repository root:  python3 dice/dice_prcusp_vs_other_sulci.py
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

OUTDIR = os.path.join("outputs", "dice")
FOCUS = "prcus-p"
THRESHOLDS = [5.0, 2.5, 10.0]                 # primary threshold first
LADDER = [0.0, 2.5, 5.0, 7.5, 10.0]
WINDOWS = [("0 mm", [0.0]), ("0-5 mm", [0.0, 2.5, 5.0]), ("0-10 mm", LADDER)]


def aggregate(d, knots):
    """Per-participant DSC for one window (trapezoidal area / window width)."""
    g = (d[d.radius_mm.isin(knots)]
         .pivot_table(index=["sub", "hemi", "sulcus"], columns="radius_mm", values="dice"))
    g = g[sorted(knots)]
    assert not g.isna().any().any(), "missing distances in the DSC table"
    if len(knots) == 1:
        v = g.iloc[:, 0]
    else:
        x = np.asarray(sorted(knots), float)
        v = pd.Series(np.trapezoid(g.values, x, axis=1) / (x[-1] - x[0]), index=g.index)
    return v.rename("stat").reset_index()


def holm(pairs):
    """Holm adjustment; pairs = [(p, name)] -> {name: adjusted p}."""
    ps, out, running = sorted(pairs), {}, 0.0
    for i, (p, nm) in enumerate(ps):
        running = max(running, p * (len(ps) - i))
        out[nm] = min(running, 1.0)
    return out


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    dice = pd.read_csv(os.path.join("data", "dice_overlap.csv"))
    geom = pd.read_csv(os.path.join("data", "face_patch_geometry.csv"))
    L, rows, windows_out = [], [], []

    for pct in THRESHOLDS:
        d = dice[dice.threshold_pct == pct]
        sulci = sorted(d.sulcus.unique())
        assert len(sulci) == 8 and d["sub"].nunique() == 47
        area = geom.loc[geom.threshold_pct == pct, "top_component_area_cm2"].median()
        L += ["=" * 78, f"TOP {pct:g}% THRESHOLD   (median area of the face-selective region "
              f"{area:.2f} cm2)", "=" * 78]

        for wname, knots in WINDOWS:
            a = aggregate(d, knots)
            windows_out.append(a.assign(threshold_pct=pct, window=wname))
            for hemi in ("lh", "rh"):
                h = a[a.hemi == hemi].pivot_table(index="sub", columns="sulcus", values="stat")
                pairs = []
                for s in sulci:
                    if s == FOCUS:
                        continue
                    diff = h[FOCUS] - h[s]
                    p = 1.0 if np.allclose(diff, 0) else float(wilcoxon(h[FOCUS], h[s]).pvalue)
                    pairs.append((p, s))
                hp = holm(pairs)
                means = h.mean().sort_values(ascending=False)
                first = int((h.idxmax(axis=1) == FOCUS).sum())
                L.append(f"\n{wname:7s} {hemi.upper()}: mean DSC " +
                         ", ".join(f"{s} {m:.3f}" for s, m in means.items()))
                L.append(f"   prcus-p vs each sulcus (Wilcoxon signed-rank; Holm across 7): " +
                         "; ".join(f"{s} p={p:.4g} holm={hp[s]:.4g}" for p, s in
                                   sorted(pairs, key=lambda t: hp[t[1]])))
                L.append(f"   significant after Holm: {sum(v < .05 for v in hp.values())}/7; "
                         f"max Holm p over significant = "
                         f"{max([v for v in hp.values() if v < .05], default=np.nan):.4f}; "
                         f"min Holm p over non-significant = "
                         f"{min([v for v in hp.values() if v >= .05], default=np.nan):.4f}; "
                         f"prcus-p highest in {first}/47 hemispheres")
                for p, s in pairs:
                    rows.append(dict(threshold_pct=pct, window=wname, hemi=hemi, sulcus=s,
                                     prcusp_mean=h[FOCUS].mean(), sulcus_mean=h[s].mean(),
                                     p_wilcoxon=p, p_holm=hp[s]))

        # mean DSC curves: which sulcus is highest at every distance, and where prcus-p peaks
        curve = d.groupby(["hemi", "radius_mm", "sulcus"])["dice"].mean().reset_index()
        L.append("")
        for hemi in ("lh", "rh"):
            c = curve[curve.hemi == hemi]
            top = c.loc[c.groupby("radius_mm")["dice"].idxmax()]
            pk = c[c.sulcus == FOCUS].set_index("radius_mm")["dice"]
            L.append(f"{hemi.upper()}: highest mean DSC at each distance: " +
                     ", ".join(f"{r:g} mm {s}" for r, s in zip(top.radius_mm, top.sulcus)) +
                     f"; prcus-p peaks at {pk.idxmax():g} mm ({pk.max():.3f})")
        z = aggregate(d, LADDER)
        z = z[z.sulcus == "mcgs"]
        L.append(f"mcgs 0-10 mm DSC = 0 in {100 * (z.stat == 0).mean():.1f}% of hemispheres")
        L.append("")

    pd.DataFrame(rows).to_csv(os.path.join(OUTDIR, "dice_prcusp_vs_other_sulci.csv"), index=False)
    pd.concat(windows_out).to_csv(os.path.join(OUTDIR, "dice_window_values.csv"), index=False)
    txt = "\n".join(L)
    with open(os.path.join(OUTDIR, "dice_prcusp_vs_other_sulci.txt"), "w") as fh:
        fh.write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
