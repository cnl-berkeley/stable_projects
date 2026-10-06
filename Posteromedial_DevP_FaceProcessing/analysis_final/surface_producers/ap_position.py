#!/usr/bin/env python3
"""
surface_producers/ap_position.py

Anterior-posterior position of the prcus-p and prculs-d labels (source of
data/sulcal_ap_position.csv; column ap_position).

For every vertex v of the white surface,
    ap(v) = geo(v -> pos) / (geo(v -> pos) + geo(v -> mcgs)),
where geo is the geodesic distance (Dijkstra on the mesh graph) to the nearest
vertex of the parieto-occipital sulcus (pos) or of the marginal ramus of the
cingulate sulcus (mcgs) label; ap = 0 at the pos and 1 at the mcgs. Each
sulcus's position is the mean of ap over all its label vertices (anchor
vertices excluded).

Requires per-participant FreeSurfer surfaces; see surface_io.py.
Usage:  python ap_position.py --id-map MAP.csv [--out ap_position.csv]
"""
import argparse, os, sys

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import dijkstra

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import surface_io as ge

ANCHOR_POST, ANCHOR_ANT = "POS", "MCGS"        # ap = 0 at pos, 1 at mcgs
SULCI = {"prcus-p": "prcus1", "prculs-d": "prculs_d"}


def multisource_geodesic(graph, sources, n):
    """Geodesic distance from EVERY vertex to the NEAREST vertex in `sources`.

    Implemented with a virtual node n (0-weight edges to each source) so this is a single
    Dijkstra call, equivalent to the minimum over sources of the per-source distance."""
    sources = np.asarray(sources)
    rows = np.r_[np.full(sources.size, n), sources]
    cols = np.r_[sources, np.full(sources.size, n)]
    ext = sparse.vstack([sparse.hstack([graph, sparse.csr_matrix((n, 1))]),
                         sparse.csr_matrix((1, n + 1))]).tolil()
    ext[rows, cols] = 1e-12            # 0 would be dropped as "no edge" by sparse storage
    d = dijkstra(ext.tocsr(), directed=False, indices=n)
    return d[:n]


def hemi_rows(sub, hemi, id_map):
    verts, faces, curv, graph = ge.load_surface(sub, hemi)
    n = verts.shape[0]
    pos = ge.label_vertices(sub, hemi, ANCHOR_POST)
    mcg = ge.label_vertices(sub, hemi, ANCHOR_ANT)
    if pos is None or mcg is None or len(pos) == 0 or len(mcg) == 0:
        raise ValueError(f"{sub} {hemi}: anchor missing")
    d_pos = multisource_geodesic(graph, pos, n)
    d_mcg = multisource_geodesic(graph, mcg, n)
    with np.errstate(invalid="ignore", divide="ignore"):
        coord = d_pos / (d_pos + d_mcg)
    anchors = np.union1d(pos, mcg)
    rows = []
    for name, stem in SULCI.items():
        L = ge.label_vertices(sub, hemi, stem)
        if L is None or len(L) == 0:
            raise ValueError(f"{sub} {hemi}: {name} label missing (expected present in all hemispheres)")
        use = np.setdiff1d(L, anchors)
        use = use[np.isfinite(coord[use])]
        rows.append(dict(sub=id_map[sub][0], hemi=hemi, group=id_map[sub][1], sulcus=name, n_label=len(L),
                         n_used=len(use), mean_ap=float(coord[use].mean())))
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--id-map", required=True, help="identifier map CSV (see surface_io.py)")
    p.add_argument("--out", default="ap_position.csv")
    a = p.parse_args()
    id_map = ge.read_id_map(a.id_map)
    rows = []
    for sub in sorted(id_map):
        for hemi in ge.HEMIS:
            rows += hemi_rows(sub, hemi, id_map)
            ge._cache.pop((sub, hemi), None)
    pd.DataFrame(rows).to_csv(a.out, index=False)


if __name__ == "__main__":
    main()
