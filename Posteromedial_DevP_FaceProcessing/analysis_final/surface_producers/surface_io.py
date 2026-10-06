#!/usr/bin/env python3
"""
surface_producers/surface_io.py

Reading of per-participant FreeSurfer data shared by the surface producers:
white-surface mesh and curvature, the geodesic mesh graph, the native
faces - objects face-selectivity map, and manually defined sulcal labels.

Expected layout under $SUBJECTS_DIR (FreeSurfer subjects directory; not shared):
  <sub>/surf/{lh,rh}.white, {lh,rh}.curv, {lh,rh}.sulc
  <sub>/label/{lh,rh}.<sulcus stem>.label, {lh,rh}.aparc.a2009s.annot
  <sub>/face_map_{lh,rh}.nii    faces - objects t map on the native surface
                                (mean across the five localizer runs)
Participant identifiers: the shared data use the labels of Figures S2-S5
(NT#/DP#). FreeSurfer subject directories use the original study codes, so the
producers take an identifier map (not shared): a CSV with columns `subject`
(FreeSurfer subject directory name), `sub` (NT#/DP# label) and `group`, passed
with --id-map. Output rows carry the NT#/DP# label.

Sulcal label stems: prcus1 = prcus-p, prcus2 = prcus-i, prcus3 = prcus-a,
prculs / prculs_d = prculs-d, pmc_2 = ifrms, pmc_x = sspls-v, MCGS = mcgs,
POS = pos, sbps = spls.
"""
import os

import nibabel as nib
import numpy as np
from scipy.sparse import csr_matrix

SUBJECTS_DIR = os.environ.get("SUBJECTS_DIR", "subjects")
HEMIS = ["lh", "rh"]
_cache = {}


def load_surface(sub, hemi):
    """verts, faces, curv, geodesic-graph (cached per sub/hemi). Read-only."""
    key = (sub, hemi)
    if key in _cache:
        return _cache[key]
    sd = os.path.join(SUBJECTS_DIR, str(sub))
    verts, faces = nib.freesurfer.read_geometry(os.path.join(sd, "surf", f"{hemi}.white"))
    curv = nib.freesurfer.read_morph_data(os.path.join(sd, "surf", f"{hemi}.curv"))
    graph = _mesh_graph(verts, faces)
    _cache[key] = (verts, faces, curv, graph)
    return _cache[key]

def _mesh_graph(verts, faces):
    """Sparse undirected graph; edge weight = euclidean length (for geodesic dijkstra)."""
    e = np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    e = np.sort(e, axis=1)
    e = np.unique(e, axis=0)
    d = np.linalg.norm(verts[e[:, 0]] - verts[e[:, 1]], axis=1)
    n = verts.shape[0]
    g = csr_matrix((np.concatenate([d, d]),
                    (np.concatenate([e[:, 0], e[:, 1]]),
                     np.concatenate([e[:, 1], e[:, 0]]))), shape=(n, n))
    return g

def load_face_map(sub, hemi):
    p = os.path.join(SUBJECTS_DIR, str(sub), f"face_map_{hemi}.nii")
    return np.asarray(nib.load(p).dataobj).squeeze()

def label_vertices(sub, hemi, stem):
    """Vertex indices for a standalone {h}.<stem>.label, or None if absent."""
    p = os.path.join(SUBJECTS_DIR, str(sub), "label", f"{hemi}.{stem}.label")
    if not os.path.exists(p):
        return None
    return np.unique(nib.freesurfer.read_label(p))


def read_id_map(path):
    """{FreeSurfer subject name: (participant label, group)} from the identifier map CSV."""
    import pandas as pd
    m = pd.read_csv(path, dtype=str)
    return {r.subject: (r.sub, r.group) for r in m.itertuples()}
