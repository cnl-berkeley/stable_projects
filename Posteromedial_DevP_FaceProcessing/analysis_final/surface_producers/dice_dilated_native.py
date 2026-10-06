#!/usr/bin/env python3
"""
surface_producers/dice_dilated_native.py

Dice similarity coefficient between each of the eight consistently present PMC
sulci and the maximally face-selective region, with the sulcal label dilated by
0-20 mm in 2.5-mm steps, on each participant's native white surface (source of
data/dice_overlap.csv; the paper uses 0-10 mm).

  DSC = 2|S_r & F| / (|S_r| + |F|) on surface area, where F is the highest-mass
  connected component of the top X% of PMC-territory vertices by faces - objects
  t and S_r is the sulcal label dilated by r mm of geodesic distance (Dijkstra
  on the white-surface mesh graph). Vertex areas are barycentric.

Requires per-participant FreeSurfer surfaces; see surface_io.py.
Usage:  python dice_dilated_native.py --id-map MAP.csv [--pct 5] [--out FILE] [--subjects ...] [--resume]
"""
import argparse, os, sys, time
import numpy as np
import pandas as pd
from scipy.sparse.csgraph import dijkstra

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import surface_io as ge
import patch_geometry as pg
from centroid_distance import top_mass_component

PCT = 5.0
RADII = [0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0]   # the paper uses 0-10 mm
# the eight PMC sulci present in every hemisphere, and their label stems
SULCI = {"prcus-p": "prcus1", "prcus-i": "prcus2", "prcus-a": "prcus3", "prculs-d": "prculs",
         "pos": "POS", "sbps": "sbps", "mcgs": "MCGS", "ifrms": "pmc_2"}
OUT = "dice_dilated_native.csv"


def run(subjects, out, resume, id_map, pct=PCT):
    if resume and os.path.exists(out):
        done = set(pd.read_csv(out, usecols=["sub"])["sub"])
        subjects = [s for s in subjects if id_map[s][0] not in done]
    elif os.path.exists(out):
        os.remove(out)
    for sub in subjects:
        t0, rows = time.time(), []
        for hemi in ("lh", "rh"):
            verts, faces, _curv, graph = ge.load_surface(sub, hemi)
            fmap = ge.load_face_map(sub, hemi)
            assert fmap.shape[0] == verts.shape[0], (sub, hemi, fmap.shape, verts.shape)
            terr = pg.territory_mask_destrieux(sub, hemi, verts.shape[0])
            areas = pg.vertex_areas(verts, faces)
            idx, _cut = pg.face_patch(fmap, terr, pct)
            comp, _ = top_mass_component(graph, idx, areas, fmap)
            F = np.zeros(verts.shape[0], bool); F[comp] = True
            aF = areas[F].sum()
            for name, stem in SULCI.items():
                S = ge.label_vertices(sub, hemi, stem)
                if S is None or len(S) == 0:
                    continue
                dd = dijkstra(graph, directed=False, indices=S, min_only=True)
                for r in RADII:
                    band = dd <= r
                    aS, aI = areas[band].sum(), areas[band & F].sum()
                    rows.append(dict(sub=id_map[sub][0], hemi=hemi, group=id_map[sub][1], sulcus=name,
                                     radius_mm=r, band_cm2=aS / 100, patch_cm2=aF / 100,
                                     capture=aI / aF, precision=aI / aS if aS > 0 else np.nan,
                                     dice=2 * aI / (aS + aF), n_vert=verts.shape[0]))
            ge._cache.pop((sub, hemi), None)
        pd.DataFrame(rows).to_csv(out, mode="a", header=not os.path.exists(out), index=False)
        print(f"  sub {sub}: {len(rows)} rows, {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--id-map", required=True, help="identifier map CSV (see surface_io.py)")
    ap.add_argument("--subjects", nargs="+", help="FreeSurfer subject names (default: all in the map)")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--pct", type=float, default=PCT,
                    help="threshold, top pct%% of the PMC territory (the paper uses 5, 2.5 and 10)")
    a = ap.parse_args()
    id_map = ge.read_id_map(a.id_map)
    run(a.subjects or sorted(id_map), a.out, a.resume, id_map, a.pct)
