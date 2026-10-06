#!/usr/bin/env python3
"""
surface_producers/patch_geometry.py

Geometry of the face-selective region in PMC for each participant and
hemisphere, at the top 2.5%, 5% and 10% thresholds.

  territory   six Destrieux parcels covering PMC (G_precuneus,
              G_cingul-Post-dorsal, G_cingul-Post-ventral, S_subparietal,
              S_parieto_occipital, S_cingul-Marginalis)
  patch       top X% of territory vertices by faces - objects t
  components  connected components of the patch on the white-surface mesh;
              mass of a component = sum over its vertices of area x t

Output columns (one row per participant, hemisphere and threshold) include the
number of territory and patch vertices, the number of components, the mass
fraction of the highest-mass component and its number of vertices. These are
the source of data/face_patch_geometry.csv. The functions territory_mask_destrieux,
vertex_areas and face_patch are also used by centroid_distance.py and
dice_dilated_native.py.

Requires per-participant FreeSurfer surfaces (not shared); see surface_io.py.
Usage:  python patch_geometry.py --id-map MAP.csv [--subjects ...] [--out patch_geometry.csv]
"""
import argparse, os, sys, time
import numpy as np
import pandas as pd
import nibabel as nib
from scipy.sparse.csgraph import connected_components

import surface_io as ge

SUBJECTS_DIR = ge.SUBJECTS_DIR
HEMIS = ge.HEMIS

# PMC territory: six Destrieux parcels, matched by name.
DESTRIEUX_PMC = ["G_cingul-Post-dorsal", "G_cingul-Post-ventral", "G_precuneus",
                 "S_subparietal", "S_parieto_occipital", "S_cingul-Marginalis"]

LADDER = [2.5, 5.0, 10.0]             # thresholds reported in the paper
MASS_FRAC_CRITERION = 0.50


def territory_mask_destrieux(sub, hemi, n):
    """Boolean PMC-territory mask from Destrieux aparc.a2009s.annot (the six PMC parcels)."""
    p = os.path.join(SUBJECTS_DIR, str(sub), "label", f"{hemi}.aparc.a2009s.annot")
    labels, _, names = nib.freesurfer.read_annot(p)
    names = [x.decode() if isinstance(x, bytes) else x for x in names]
    missing = [w for w in DESTRIEUX_PMC if w not in names]
    if missing:
        raise ValueError(f"sub {sub} {hemi}: Destrieux parcels missing from annot: {missing}")
    out = np.zeros(n, bool)
    for w in DESTRIEUX_PMC:
        out |= (labels == names.index(w))
    return out


def vertex_areas(verts, faces):
    """Barycentric per-vertex surface area (each triangle's area split 1/3 to each vertex)."""
    v0, v1, v2 = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    tri = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)
    a = np.zeros(verts.shape[0])
    for k in range(3):
        np.add.at(a, faces[:, k], tri / 3.0)
    return a


def face_patch(fmap, territory, pct):
    """Top-pct% of territory vertices by signed faces-objects t (percentile WITHIN territory)."""
    terr = np.where(territory)[0]
    if terr.size == 0:
        return np.array([], dtype=int), np.nan
    vals = fmap[terr]
    cut = float(np.percentile(vals, 100.0 - pct))
    return np.sort(terr[vals >= cut]), cut


def component_stats(graph, idx, areas, fmap):
    """Connected components of `idx` on the mesh graph, with mass and area of each."""
    idx = np.asarray(idx)
    if idx.size == 0:
        return dict(n_comp=0, top_mass_frac=np.nan, top_area_frac=np.nan,
                    concordant=np.nan, top_mass_n=0, top_mass_area=np.nan,
                    total_mass=np.nan, min_t=np.nan, n_neg_t=0)
    sub = graph[idx][:, idx]
    n_comp, lab = connected_components(sub, directed=False)
    mass = np.zeros(n_comp)
    area = np.zeros(n_comp)
    cnt = np.zeros(n_comp, int)
    w = areas[idx] * fmap[idx]
    for c in range(n_comp):
        m = lab == c
        mass[c] = w[m].sum()
        area[c] = areas[idx][m].sum()
        cnt[c] = int(m.sum())
    i_mass = int(np.argmax(mass))
    i_area = int(np.argmax(area))
    return dict(
        n_comp=int(n_comp),
        top_mass_frac=float(mass[i_mass] / mass.sum()) if mass.sum() > 0 else np.nan,
        top_area_frac=float(area[i_area] / area.sum()),
        concordant=bool(i_mass == i_area),
        top_mass_n=int(cnt[i_mass]),
        top_mass_area=float(area[i_mass]),
        total_mass=float(mass.sum()),
        min_t=float(fmap[idx].min()),
        n_neg_t=int((fmap[idx] < 0).sum()),
    )


def run(subjects, out_csv, id_map):
    rows = []
    for sub in subjects:
        t0 = time.time()
        for hemi in HEMIS:
            if not os.path.isdir(os.path.join(SUBJECTS_DIR, str(sub))):
                continue
            try:
                verts, faces, curv, graph = ge.load_surface(sub, hemi)
                fmap = ge.load_face_map(sub, hemi)
                terr = territory_mask_destrieux(sub, hemi, verts.shape[0])
                areas = vertex_areas(verts, faces)
            except Exception as ex:
                print(f"  [skip] sub {sub} {hemi}: {ex}", file=sys.stderr)
                continue
            if fmap.shape[0] != verts.shape[0]:
                print(f"  [skip] sub {sub} {hemi}: face map {fmap.shape} != surface {verts.shape[0]}",
                      file=sys.stderr)
                continue
            for pct in LADDER:
                idx, cut = face_patch(fmap, terr, pct)
                st = component_stats(graph, idx, areas, fmap)
                rows.append(dict(sub=id_map[sub][0], group=id_map[sub][1], hemi=hemi, pct=pct, n_terr=int(terr.sum()),
                                 n_patch=len(idx), t_cut=cut, **st))
            ge._cache.pop((sub, hemi), None)
        print(f"  done sub {sub}  ({time.time()-t0:.1f}s)", file=sys.stderr)
    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"\nwrote {out_csv}  ({len(df)} rows)")
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--id-map", required=True, help="identifier map CSV (see surface_io.py)")
    ap.add_argument("--subjects", nargs="+", help="FreeSurfer subject names (default: all in the map)")
    ap.add_argument("--out", default="patch_geometry.csv")
    a = ap.parse_args()
    id_map = ge.read_id_map(a.id_map)
    run(a.subjects or sorted(id_map), a.out, id_map)
