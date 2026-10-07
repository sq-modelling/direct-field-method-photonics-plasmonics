#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Qiang Sun
# SPDX-License-Identifier: BSD-3-Clause
"""Gap-graded curved six-node (Q6) surface meshes of the Fig. 3 'LgRd' rods, written in the
Prtl_*.msh import format of the release code (MeshRead = 1, mesh type 'Q').

Geometry: exactly the release 'LgRd' Rankine ovoid (rankine.py), unzoomed local frame, axis z,
tips at z = +-1, half-width b.  The release input still applies zoom (1.5 or 1.0), the pitch rotation
(0, 90, 0) that maps local z to global x, and the translation, so the imported rods coincide with the
generated ones (only the node distribution differs).

Mesh: body-of-revolution ring mesh.  Corner nodes lie on rings s = s_k of the meridian arc length
(a single node at each tip).  Target corner-edge length (physical units, a = 1):
    h(s) = min(h_side, h_tip + grad * d_gap(s), h_far + grad * d_far(s))
where d_gap is the meridian distance to a gap-facing tip and d_far to a free tip; the azimuthal node
count of a ring is round(2 pi rho / h) and the ring spacing is 0.866 h * aspect(s), with
aspect = 1 on the rounded ends and 'aspect' on the nearly cylindrical side (anisotropic side elements
like the release Icshdrl rods).  Neighbouring rings are joined by an angular zipper; every corner edge
gets a mid-side node at the parameter midpoint (s, phi), mapped exactly onto the surface.

usage (library): build_rod(b, zoom, gap_top, gap_bot, h_tip, h_side, grad, aspect, h_far) -> X, E
       (command line): meshgen_rod.py OUTDIR h_tip h_side grad aspect [h_far]
       writes OUTDIR/Prtl_0001.msh (rod 1: gap at local z=+1), Prtl_0002.msh (rod 2: gap at z=-1),
       Prtl_0003.msh (rod 3: gaps at both tips) and OUTDIR/meshinfo.txt.

The meshes supplied in ../ka5, ../ka0001 and ../ka50 (Fig. 3) are reproduced byte for byte by
       python3 meshgen_rod.py OUTDIR 0.005 0.08 0.3 1.0     (ka5, ka0001: 11 878 nodes)
       python3 meshgen_rod.py OUTDIR 0.005 0.04 0.3 1.0     (ka50: 41 798 nodes)
i.e. element size g = 0.005a at the gap-facing tips, growing at rate 0.3 to 0.08a or 0.04a; meshinfo_*.txt
are the corresponding check outputs.  Requires NumPy and SciPy (about 25 s for the ka50 set).
"""
import os, sys
import numpy as np
from rankine import Rankine

_RCACHE = {}


def rankine(b):
    if b not in _RCACHE:
        _RCACHE[b] = Rankine(b)
    return _RCACHE[b]


def surf(R, s, phi, zoom):
    """physical meridian arc length s (from local z=+1 tip) and azimuth -> unzoomed local xyz"""
    r, z = R.point_at_s(np.atleast_1d(s) / zoom)
    phi = np.atleast_1d(phi)
    return np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=1)


def build_rod(b, zoom, gap_top, gap_bot, h_tip, h_side, grad, aspect=1.0, h_far=None):
    R = rankine(b)
    S = R.S * zoom                                   # physical meridian length tip to tip
    h_far = h_side if h_far is None else h_far
    r0 = b * zoom

    def hfun(s):
        s = np.asarray(s, float)
        h = np.full_like(s, h_side)
        dt, db = s, S - s
        h = np.minimum(h, (h_tip if gap_top else h_far) + grad * dt)
        h = np.minimum(h, (h_tip if gap_bot else h_far) + grad * db)
        return h

    # rho(s) table (physical)
    sg = np.linspace(0, S, 20001)
    rg, zg = R.point_at_s(sg / zoom)
    rg = rg * zoom

    def asp(s):
        d = np.minimum(s, S - s)
        w = np.clip((d - 0.35) / 0.35, 0, 1)          # 1 on the rounded ends, 'aspect' on the side
        w = w * w * (3 - 2 * w)
        return 1 + (aspect - 1) * w

    # ring positions: integrate dk/ds = 1/(0.866 h aspect), then equalise to an integer count
    dk = 1.0 / (0.866 * hfun(sg) * asp(sg))
    k = np.concatenate([[0], np.cumsum(0.5 * (dk[1:] + dk[:-1]) * np.diff(sg))])
    M = max(4, int(round(k[-1])))
    s_r = np.interp(np.linspace(0, k[-1], M + 1), k, sg)
    s_r[0], s_r[-1] = 0.0, S
    # azimuthal counts
    rho_r = np.interp(s_r, sg, rg)
    n_r = np.maximum(5, np.round(2 * np.pi * rho_r / hfun(s_r)).astype(int))
    n_r[0] = n_r[-1] = 1
    # corner nodes: tip, rings, tip; staggered angular offsets
    nodes_s, nodes_p, ring_ids = [], [], []
    off = 0.0
    for kk in range(M + 1):
        n = n_r[kk]
        if n == 1:
            ids = [len(nodes_s)]; nodes_s.append(s_r[kk]); nodes_p.append(np.nan)
        else:
            ph = off + 2 * np.pi * np.arange(n) / n
            ids = list(range(len(nodes_s), len(nodes_s) + n))
            nodes_s += [s_r[kk]] * n; nodes_p += list(ph)
            off = off + np.pi / n                       # stagger the next ring by half a spacing
        ring_ids.append(ids)
    nodes_s = np.array(nodes_s); nodes_p = np.array(nodes_p)
    tris = []
    for kk in range(M):
        A, B = ring_ids[kk], ring_ids[kk + 1]
        if len(A) == 1:                                 # fan from the top tip
            for j in range(len(B)):
                tris.append((A[0], B[j], B[(j + 1) % len(B)]))
            continue
        if len(B) == 1:                                 # fan to the bottom tip
            for i in range(len(A)):
                tris.append((B[0], A[(i + 1) % len(A)], A[i]))
            continue
        na, nb = len(A), len(B)
        pa = nodes_p[A][0] + np.concatenate([[0.0], np.cumsum(np.diff(nodes_p[A]) % (2 * np.pi))])
        pa = np.append(pa, pa[0] + 2 * np.pi)              # pa[na] = A[0] again
        # ring B re-indexed to start at the node closest in angle to A[0], angles unwrapped near pa[0]
        dphi = (nodes_p[B] - pa[0] + np.pi) % (2 * np.pi) - np.pi
        j0 = int(np.argmin(np.abs(dphi)))
        Bi = [B[(j0 + j) % nb] for j in range(nb)]
        pb = pa[0] + dphi[j0] + np.concatenate([[0.0], np.cumsum(np.diff(nodes_p[Bi]) % (2 * np.pi))])
        pb = np.append(pb, pb[0] + 2 * np.pi)
        i = j = 0
        while i < na or j < nb:
            if j == nb or (i < na and pa[i + 1] <= pb[j + 1]):
                tris.append((A[i % na], Bi[j % nb], A[(i + 1) % na])); i += 1
            else:
                tris.append((A[i % na], Bi[j % nb], Bi[(j + 1) % nb])); j += 1
    tris = np.array(tris)
    Xc = surf(R, nodes_s, np.nan_to_num(nodes_p), zoom)
    # orient outward (convex body: centroid . n > 0)
    P = Xc[tris]
    nrm = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    flip = np.einsum('ij,ij->i', nrm, P.mean(1)) < 0
    tris[flip] = tris[flip][:, [0, 2, 1]]
    # mid-side nodes
    edge_mid = {}
    mid_s, mid_p = [], []
    E = np.zeros((len(tris), 6), int)
    nc = len(nodes_s)
    for e, (a, b_, c) in enumerate(tris):
        mids = []
        for (u, v) in ((a, b_), (b_, c), (c, a)):
            key = (min(u, v), max(u, v))
            if key not in edge_mid:
                su, sv = nodes_s[u], nodes_s[v]
                pu, pv = nodes_p[u], nodes_p[v]
                if np.isnan(pu):
                    pm = pv
                elif np.isnan(pv):
                    pm = pu
                else:
                    pm = pu + (((pv - pu) + np.pi) % (2 * np.pi) - np.pi) / 2
                edge_mid[key] = nc + len(mid_s)
                mid_s.append(0.5 * (su + sv)); mid_p.append(pm)
            mids.append(edge_mid[key])
        E[e] = [a, b_, c, mids[0], mids[1], mids[2]]
    Xm = surf(R, np.array(mid_s), np.array(mid_p), zoom)
    X = np.vstack([Xc, Xm])
    info = dict(M=M, n_r=n_r, s_r=s_r, ncorner=nc, nnodes=len(X), nelem=len(E))
    return X, E, info


def _q6_rule():
    """Standard isoparametric six-node triangle (mid-side nodes at the parametric midpoints; node order
    1 = (0,0), 2 = (1,0), 3 = (0,1), 4 = mid 1-2, 5 = mid 2-3, 6 = mid 1-3), used only for the area and
    volume checks: shape functions and derivatives at a 28-point rule (each element split into four
    parametric sub-triangles with the 7-point degree-5 rule).  Returns N, dN/dxi, dN/deta (28, 6), weights."""
    a1, b1 = 0.059715871789770, 0.470142064105115
    a2, b2 = 0.797426985353087, 0.101286507323456
    w0, w1, w2 = 0.225, 0.132394152788506, 0.125939180544827
    bary = [(1/3, 1/3, 1/3, w0)] + [(a1, b1, b1, w1), (b1, a1, b1, w1), (b1, b1, a1, w1)] + \
           [(a2, b2, b2, w2), (b2, a2, b2, w2), (b2, b2, a2, w2)]
    subs = [((0, 0), (.5, 0), (0, .5)), ((.5, 0), (1, 0), (.5, .5)), ((0, .5), (.5, .5), (0, 1)),
            ((.5, 0), (.5, .5), (0, .5))]
    pts, wts = [], []
    for (p0, p1, p2) in subs:
        for (l0, l1, l2, w) in bary:
            pts.append((l0*p0[0] + l1*p1[0] + l2*p2[0], l0*p0[1] + l1*p1[1] + l2*p2[1]))
            wts.append(w * 0.5 * 0.25)
    pts = np.array(pts); wts = np.array(wts)
    xi, et = pts[:, 0], pts[:, 1]; L = 1 - xi - et
    N = np.stack([L*(2*L-1), xi*(2*xi-1), et*(2*et-1), 4*L*xi, 4*xi*et, 4*L*et], axis=1)
    dNx = np.stack([-(4*L-1), 4*xi-1, 0*xi, 4*(L-xi), 4*et, -4*et], axis=1)
    dNe = np.stack([-(4*L-1), 0*xi, 4*et-1, -4*xi, 4*xi, 4*(L-et)], axis=1)
    return N, dNx, dNe, wts


def validate(X, E, b, zoom):
    """closedness, orientation consistency, Euler characteristic, area/volume vs exact, surface error"""
    from collections import Counter
    R = rankine(b)
    ed = Counter()
    for t in E[:, :3]:
        for u, v in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            ed[(u, v)] += 1
    closed = all(ed.get((v, u), 0) == 1 and c == 1 for (u, v), c in ed.items())
    nc = len(np.unique(E[:, :3])); ne = len(ed) // 2; nf = len(E)
    euler = nc - ne + nf
    # mid nodes shared by exactly two elements
    cm = Counter(E[:, 3:].ravel().tolist())
    midok = all(v == 2 for v in cm.values())
    # exact surface error of all nodes (unzoomed): F/|grad F| ~ distance
    rr = np.hypot(X[:, 0], X[:, 1]); zz = X[:, 2]
    dist = []
    for r_, z_ in zip(rr, zz):
        pr, pz = R.project(r_, z_)
        dist.append(np.hypot(pr - r_, pz - z_))
    # area / volume with the standard isoparametric Q6 (7-pt rule on 4 subtriangles) in physical units
    NN, dNx, dNe, WTS = _q6_rule()
    P = X[E] * zoom
    r = np.einsum('gj,mjd->mgd', NN, P)
    rx = np.einsum('gj,mjd->mgd', dNx, P); re = np.einsum('gj,mjd->mgd', dNe, P)
    nv = np.cross(rx, re); J = np.linalg.norm(nv, axis=2); n = nv / J[..., None]
    w = WTS[None, :] * J
    area = w.sum(); vol = np.sum(w * np.sum(r * n, axis=2)) / 3
    # exact area/volume of the Rankine body of revolution
    s = np.linspace(0, R.S, 200001); rs, zs = R.point_at_s(s)
    area_ex = 2 * np.pi * np.trapezoid(rs, s) * zoom**2
    vol_ex = np.pi * np.trapezoid(rs**2, -zs) * zoom**3
    # triangle quality (corner triangles, physical)
    C = X[E[:, :3]] * zoom
    L = np.stack([np.linalg.norm(C[:, 1] - C[:, 0], axis=1), np.linalg.norm(C[:, 2] - C[:, 1], axis=1),
                  np.linalg.norm(C[:, 0] - C[:, 2], axis=1)], 1)
    ang = []
    for i in range(3):
        a_, b2, c_ = L[:, i], L[:, (i + 1) % 3], L[:, (i + 2) % 3]
        ang.append(np.degrees(np.arccos(np.clip((b2**2 + c_**2 - a_**2) / (2 * b2 * c_), -1, 1))))
    ang = np.array(ang)
    return dict(closed=closed, mid_shared_twice=midok, euler=euler, outward=vol > 0, area=area, area_exact=area_ex,
                vol=vol, vol_exact=vol_ex, max_surf_dist=max(dist) * zoom, min_angle=ang.min(),
                max_edge=L.max(), min_edge=L.min())


def write_msh(fn, X, E, title):
    with open(fn, 'w') as f:
        f.write(' $MeshFormat\n $quadratic\n $EndMeshFormat\n $PhysicalNames\n $\n "%s"\n $EndPhysicalNames\n $Nodes\n'
                % title)
        f.write('%12d\n' % len(X))
        for i, p in enumerate(X, 1):
            f.write('%12d %26.17E %26.17E %26.17E\n' % (i, p[0], p[1], p[2]))
        f.write(' $EndNodes\n $Elements\n%12d\n' % len(E))
        for k, e in enumerate(E, 1):
            f.write('%9d%9d%9d%9d%9d' % (k, k, k, k, k) + ''.join('%9d' % (x + 1) for x in e) + '\n')
        f.write(' $EndElements\n')


TRANSL = (-2.505, 2.505, 0.0)   # release x offsets of particles 1, 2, 3


def to_global(X, i):
    """local unzoomed rod mesh of particle i (1..3) -> final frame of the release input
    (zoom, pitch rotation (0, 90, 0): x' = z, y' = y, z' = -x, then translation)"""
    b, zoom = RODS[i - 1][:2]
    Y = np.empty_like(X)
    Y[:, 0] = zoom * X[:, 2] + TRANSL[i - 1]; Y[:, 1] = zoom * X[:, 1]; Y[:, 2] = -zoom * X[:, 0]
    return Y


RODS = [  # (b, zoom, gap at local z=+1, gap at local z=-1) for release particles 1, 2, 3
    (0.2, 1.5, True, False),
    (0.2, 1.5, False, True),
    (0.3, 1.0, True, True)]


def make_set(outdir, h_tip, h_side, grad, aspect=1.0, h_far=None, quiet=False):
    os.makedirs(outdir, exist_ok=True)
    lines = []
    for i, (b, zoom, gt, gb) in enumerate(RODS, 1):
        if i == 2:      # rod 2 = mirror image (local z -> -z) of rod 1, so the mesh pair keeps the x -> -x symmetry
            X = X1.copy(); X[:, 2] *= -1; E = E1[:, [1, 0, 2, 3, 5, 4]].copy(); info = info1
        else:
            X, E, info = build_rod(b, zoom, gt, gb, h_tip, h_side, grad, aspect, h_far)
        if i == 1:
            X1, E1, info1 = X, E, info
        v = validate(X, E, b, zoom)
        write_msh(os.path.join(outdir, 'Prtl_%04d.msh' % i), X, E,
                  'three-rods gap-graded LgRd b=%g zoom=%g h_tip=%g h_side=%g grad=%g aspect=%g' %
                  (b, zoom, h_tip, h_side, grad, aspect))
        s = ('rod %d: N=%d elems=%d rings=%d closed=%s mids2=%s euler=%d outward=%s area=%.6f (exact %.6f) '
             'vol=%.6f (exact %.6f) max|node-surface|=%.1e min angle=%.1f deg edge %.4f..%.4f a'
             % (i, info['nnodes'], info['nelem'], info['M'], v['closed'], v['mid_shared_twice'], v['euler'],
                v['outward'], v['area'], v['area_exact'], v['vol'], v['vol_exact'], v['max_surf_dist'],
                v['min_angle'], v['min_edge'], v['max_edge']))
        lines.append(s)
        if not quiet:
            print(s)
    open(os.path.join(outdir, 'meshinfo.txt'), 'w').write(
        'h_tip=%g h_side=%g grad=%g aspect=%g h_far=%s\n' % (h_tip, h_side, grad, aspect, h_far) + '\n'.join(lines) + '\n')
    return lines


if __name__ == '__main__':
    out = sys.argv[1]
    h_tip, h_side, grad, aspect = map(float, sys.argv[2:6])
    h_far = float(sys.argv[6]) if len(sys.argv) > 6 else None
    make_set(out, h_tip, h_side, grad, aspect, h_far)
