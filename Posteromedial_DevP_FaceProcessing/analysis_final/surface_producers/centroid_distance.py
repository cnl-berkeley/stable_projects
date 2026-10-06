#!/usr/bin/env python3
"""
surface_producers/centroid_distance.py

Geodesic distance, on each participant's native white surface, from the
centroid of the maximally face-selective PMC region to the centroid of each
sulcus (source of data/centroid_distance.csv; column geodesic_mm).

  face-selective region  highest-mass connected component of the top X% of
                         PMC-territory vertices by faces - objects t
                         (patch_geometry.face_patch / top_mass_component)
  region centroid        area-weighted centroid, snapped to the nearest vertex
  sulcal centroid        area-weighted label centroid from the CNL scalpel
                         package (label_centroid, white surface), snapped to the
                         nearest vertex (https://github.com/b-parker/CNL_scalpel;
                         pass its src directory with --scalpel-src)
  distance               Dijkstra on the white-surface mesh graph

Requires per-participant FreeSurfer surfaces; see surface_io.py.
Usage:  python centroid_distance.py --id-map MAP.csv --face-pct 5 2.5 10 [--subjects ...] [--scalpel-src PATH]
"""
import argparse, os, sys, time
import numpy as np
import pandas as pd
import nibabel as nib
from scipy.sparse.csgraph import dijkstra, connected_components

import surface_io as ge
import patch_geometry as pg

SUBJECTS_DIR = ge.SUBJECTS_DIR
HEMIS = ge.HEMIS

# Sulci and their label stems (prculs = prculs-d; pmc_x = sspls-v; pmc_2 = ifrms).
REPORTED = {"prcus-p": "prcus1", "prcus-i": "prcus2", "prcus-a": "prcus3",
            "sspls-v": "pmc_x", "ifrms": "pmc_2", "mcgs": "MCGS", "pos": "POS",
            "prculs-d": "prculs"}


def _scalpel_subject(sub, hemi, scalpel_src):
    """ScalpelSubject on the WHITE surface (never the 'inflated' default)."""
    if scalpel_src and scalpel_src not in sys.path:
        sys.path.insert(0, scalpel_src)
    from scalpel.subject import ScalpelSubject
    return ScalpelSubject(str(sub), hemi, SUBJECTS_DIR, surface_type="white")


def sulcal_centroid(scalpel_sub, verts, stem):
    """CNL scalpel label_centroid -> (vertex index, RAS). None if the label is absent/degenerate."""
    try:
        scalpel_sub.load_label(stem)
    except Exception:
        return None, None
    if stem not in scalpel_sub.labels or len(scalpel_sub.labels[stem].vertex_indexes) == 0:
        return None, None
    try:
        ras = scalpel_sub.label_centroid(stem, load=False)
    except Exception as ex:
        print(f"    [centroid fail] {stem}: {ex}", file=sys.stderr)
        return None, None
    ras = np.asarray(ras).reshape(3)
    return int(np.argmin(np.linalg.norm(verts - ras, axis=1))), ras


def top_mass_component(graph, idx, areas, fmap):
    """(component vertex indices, concordant_with_area) for the HIGHEST-MASS component."""
    idx = np.asarray(idx)
    if idx.size == 0:
        return np.array([], dtype=int), np.nan
    n_comp, lab = connected_components(graph[idx][:, idx], directed=False)
    w = areas[idx] * fmap[idx]
    mass = np.array([w[lab == c].sum() for c in range(n_comp)])
    area = np.array([areas[idx][lab == c].sum() for c in range(n_comp)])
    i_mass, i_area = int(np.argmax(mass)), int(np.argmax(area))
    return idx[lab == i_mass], bool(i_mass == i_area)


def area_weighted_centroid_vertex(verts, areas, idx, restrict_to_set=False):
    """Area-weighted 3D centroid of `idx`, snapped to the nearest vertex.

    Snap matches CNL scalpel (nearest vertex over the WHOLE surface) so the face centroid and the
    sulcal centroid are defined the same way. `restrict_to_set=True` gives the in-set variant as a
    QC comparison -- not the reported definition."""
    idx = np.asarray(idx)
    if idx.size == 0:
        return None, None
    w = areas[idx]
    if w.sum() <= 0:
        w = np.ones_like(w)
    c = (verts[idx] * w[:, None]).sum(axis=0) / w.sum()
    if restrict_to_set:
        return int(idx[np.argmin(np.linalg.norm(verts[idx] - c, axis=1))]), c
    return int(np.argmin(np.linalg.norm(verts - c, axis=1))), c


def _flush(rows, out_csv):
    """Append one participant's rows, so an interrupted run can be resumed with --resume."""
    if not rows:
        return
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    hdr = not os.path.exists(out_csv)
    pd.DataFrame(rows).to_csv(out_csv, mode="a", header=hdr, index=False)


def run(subjects, pcts, out_csv, scalpel_src, id_map, resume=False):
    if resume and os.path.exists(out_csv):
        have = set(pd.read_csv(out_csv, usecols=["sub"])["sub"].unique().tolist())
        skip = [s for s in subjects if id_map[s][0] in have]
        subjects = [s for s in subjects if id_map[s][0] not in have]
        print(f"resume: {len(skip)} subjects already in {os.path.basename(out_csv)}, "
              f"{len(subjects)} to go", file=sys.stderr)
    elif os.path.exists(out_csv):
        os.remove(out_csv)      # not resuming -> start clean, never append to a stale run
    for sub in subjects:
        t0 = time.time()
        rows = []
        for hemi in HEMIS:
            if not os.path.isdir(os.path.join(SUBJECTS_DIR, str(sub))):
                continue
            try:
                verts, faces, curv, graph = ge.load_surface(sub, hemi)
                fmap = ge.load_face_map(sub, hemi)
                terr = pg.territory_mask_destrieux(sub, hemi, verts.shape[0])
                areas = pg.vertex_areas(verts, faces)
                ssub = _scalpel_subject(sub, hemi, scalpel_src)
            except Exception as ex:
                print(f"  [skip] sub {sub} {hemi}: {ex}", file=sys.stderr)
                continue

            # --- face-selective region, one per threshold ---
            targets = {}
            for pct in pcts:
                idx, cut = pg.face_patch(fmap, terr, pct)
                comp, concord = top_mass_component(graph, idx, areas, fmap)
                fcv, _ = area_weighted_centroid_vertex(verts, areas, comp)
                # the whole-surface snap can land outside the component; recorded per row
                in_comp = (fcv in set(comp.tolist())) if fcv is not None else None
                targets[pct] = dict(patch=idx, comp=comp, fcv=fcv, concord=concord, t_cut=cut,
                                    in_comp=in_comp)

            # --- sulcal centroids (scalpel, white) ---
            cents = {}
            for name, stem in REPORTED.items():
                S = ge.label_vertices(sub, hemi, stem)
                if S is None or len(S) == 0:
                    cents[name] = (None, None, None)
                    continue
                cv, _ = sulcal_centroid(ssub, verts, stem)
                cents[name] = (S, cv, (cv in set(S.tolist())) if cv is not None else None)

            # --- distances: one Dijkstra per sulcal centroid, reused across thresholds ---
            for name, (S, scv, in_label) in cents.items():
                dmap = (dijkstra(graph, directed=False, indices=int(scv))
                        if scv is not None else None)
                for pct in pcts:
                    tg = targets[pct]
                    geo = float(dmap[tg["fcv"]]) if (dmap is not None and tg["fcv"] is not None) else np.nan
                    # secondary: selectivity-weighted mean geodesic distance to ALL suprathreshold verts
                    if dmap is not None and len(tg["patch"]):
                        w = fmap[tg["patch"]] * areas[tg["patch"]]
                        d = dmap[tg["patch"]]
                        ok = np.isfinite(d) & (w > 0)
                        wmean = float((d[ok] * w[ok]).sum() / w[ok].sum()) if ok.any() else np.nan
                    else:
                        wmean = np.nan
                    rows.append(dict(
                        sub=id_map[sub][0], group=id_map[sub][1], hemi=hemi, sulcus=name, face_pct=pct,
                        present=S is not None, n_sulcus=(0 if S is None else len(S)),
                        t_cut=tg["t_cut"], n_patch=len(tg["patch"]), n_top_comp=len(tg["comp"]),
                        mass_area_concordant=tg["concord"],
                        sulc_centroid_vert=scv, sulc_centroid_in_label=in_label,
                        face_centroid_vert=tg["fcv"], face_centroid_in_comp=tg["in_comp"],
                        geodesic_mm=geo,
                        weighted_geodesic_mm=wmean,
                        euclid_mm=(float(np.linalg.norm(verts[scv] - verts[tg["fcv"]]))
                                   if (scv is not None and tg["fcv"] is not None) else np.nan),
                    ))
            ge._cache.pop((sub, hemi), None)
        _flush(rows, out_csv)
        print(f"  done sub {sub} ({time.time()-t0:.1f}s, {len(rows)} rows flushed)", file=sys.stderr)

    df = pd.read_csv(out_csv)
    print(f"\nwrote {out_csv}  ({len(df)} rows)")
    bad = df[(df.present) & (df.sulc_centroid_in_label == False)]
    if len(bad):
        u = bad[["sub", "hemi", "sulcus"]].drop_duplicates()
        print(f"[QC] {len(u)} sulcal centroids landed OUTSIDE their own label "
              f"(scalpel snaps over the whole surface): {u.sulcus.value_counts().to_dict()}")
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--id-map", required=True, help="identifier map CSV (see surface_io.py)")
    ap.add_argument("--subjects", nargs="+", help="FreeSurfer subject names (default: all in the map)")
    ap.add_argument("--face-pct", nargs="+", type=float, required=True,
                    help="threshold(s), in %% (the paper uses 5 2.5 10)")
    ap.add_argument("--scalpel-src", default=os.environ.get("SCALPEL_SRC", ""),
                    help="path to CNL_scalpel/src (or put it on PYTHONPATH)")
    ap.add_argument("--out", default="centroid_distance.csv")
    ap.add_argument("--resume", action="store_true",
                    help="keep an existing --out and skip participants already in it")
    a = ap.parse_args()
    id_map = ge.read_id_map(a.id_map)
    subjects = a.subjects or sorted(id_map)
    print(f"subjects n={len(subjects)}  pcts={a.face_pct}  sulci={list(REPORTED)}", file=sys.stderr)
    run(subjects, a.face_pct, a.out, a.scalpel_src, id_map, resume=a.resume)
