#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Qiang Sun
# SPDX-License-Identifier: BSD-3-Clause
"""Fields and radar cross section of perfect conductors from the surface current J = n_out x H.

At very low frequencies the electric-field formulation for perfect conductors leaves the net charge of
each conductor only weakly determined (manuscript, Secs. 2.3 and 7.3), and the electric field written by
`field_solver` (surface fields, planar samples and RCS cuts) should then not be used.  For the conductors
of Fig. 3 the magnetic-field solution, which `field_solver` always computes for a plane wave, does not
have this problem.  This script evaluates the scattered field from the surface current density of that
solution, as for the case k_out a = 1e-3 of Fig. 3:

    J   = n_out x H_tot                       (n_out: outward unit normal of the conductor),
    rho = div_s J / (i k)                     (surface continuity equation),
    E_sc(x) = i k int G J dS - grad int G rho dS,         G = exp(ikR)/(4 pi R),
    far field  F = (i k / 4 pi) (I - r r) int J exp(-i k r.y) dS,   sigma = 4 pi |F|^2

(manuscript Eq. (39a); units with a = 1, eps = mu = 1, E0 = 1, omega = k, time dependence exp(-i omega t)).
The surface divergence is taken element by element on the curved six-node elements of the solver,
|A| div_s(n x H) = H_eta . r_xi - H_xi . r_eta, which needs no line charges because the tangential H
is continuous across element edges.  Integrals use the solver's 16-point triangle rule; for targets
closer than four element diameters the element is subdivided adaptively until every sub-triangle is
more than three of its own diameters away (at most 10 levels).

The solver writes only the real parts of the complex phasors.  The imaginary parts are obtained from a
second run of the same case in which the incident field is multiplied by i (incident order 3 with phase
pi/2 in Input_Phys_EM.dat): that run writes Re(i F) = -Im F, so that F = F_run0 - i F_run90.

Usage (from the repository root; see README.md, Figure 3):
  python3 postprocessing/pec_current_route.py prepare CASE_DIR RUN90_DIR
      copies the case directory CASE_DIR to RUN90_DIR with the incident field multiplied by i and the
      planar samples switched off;
  python3 postprocessing/pec_current_route.py compute RUN0_DIR RUN90_DIR OUT_DIR [--no-slice]
      reads Rslt_SurfCal3_EM.dat of both runs (and Rslt_Dmn3_plot.dat and RsltPos.dat of RUN0 for the
      planar samples) and writes to OUT_DIR:
        rcs_far_zx.csv, rcs_far_xy.csv, rcs_far_yz.csv  theta (deg, 1-degree steps, the angle convention of
            Rslt_RCS_*.dat: zx: x = cos t, z = sin t; xy: x = cos t, y = sin t; yz: y = cos t, z = sin t),
            sigma/a^2 in the far-field limit, and 10 log10(sigma/a^2);
        slice_<plane>.csv  x, y, z, pos (0 outside the conductors), and Re of the total field
            E_x, E_y, E_z (E0 units, 0 inside the conductors), at the planar samples of RUN0;
        summary.txt  charge and dipole moments of each conductor and RCS values.
Requirements: Python 3 with NumPy and SciPy.  Only an outer medium of refractive index 1 is supported.
The planar samples of Fig. 3 (24 000 points) take a few minutes.

The real part of the field is the in-phase response, which is accurate at low frequency; the imaginary
part of the field computed in this way contains the discretisation error of the static current divided
by k and is not used for Fig. 3 (README.md, Figure 3).
"""
import os, re, sys, shutil
import numpy as np
from scipy.spatial import cKDTree

# 16-point degree-8 symmetric triangle rule of the solver (Pre_GQ_FDM.f90); weights sum to 1/2
TRI16_XI = np.array([0.33333333333333333333, 0.45929258829272315603, 0.08141482341455368794,
                     0.45929258829272315603, 0.17056930775176020662, 0.65886138449647958676,
                     0.17056930775176020662, 0.05054722831703097546, 0.89890554336593804908,
                     0.05054722831703097546, 0.00839477740995760534, 0.72849239295540428124,
                     0.26311282963463811342, 0.72849239295540428124, 0.26311282963463811342,
                     0.00839477740995760534])
TRI16_ETA = np.array([0.33333333333333333333, 0.08141482341455368794, 0.45929258829272315603,
                      0.45929258829272315603, 0.65886138449647958676, 0.17056930775176020662,
                      0.17056930775176020662, 0.89890554336593804908, 0.05054722831703097546,
                      0.05054722831703097546, 0.72849239295540428124, 0.00839477740995760534,
                      0.72849239295540428124, 0.26311282963463811342, 0.00839477740995760534,
                      0.26311282963463811342])
TRI16_W = np.array([0.07215780383889358413] + [0.04754581713364231240] * 3 + [0.05160868526735912514] * 3
                   + [0.01622924881159904016] * 3 + [0.01361515708721749713] * 6)
PWE_PHASE0 = '"pwe", "p", 0, 0.0d0,'
PWE_PHASE90 = '"pwe", "p", 3, 1.5707963267948966d0,'


# ----------------------------------------------------------------------------------------------
# input and output files of field_solver
def _floats(rows):
    return np.array(' '.join(rows).replace(',', ' ').split(), float).reshape(len(rows), -1)


def _zones(path):
    """Tecplot FEPOINT zones of a solver output file: list of (title, node array, triangle array)"""
    lines = open(path).read().splitlines()
    out, i = [], 0
    while i < len(lines):
        if lines[i].lstrip().startswith('Zone'):
            h = lines[i]
            n = int(re.search(r'n=\s*(\d+)', h).group(1)); e = int(re.search(r'e=\s*(\d+)', h).group(1))
            title = h.split('"')[1].strip() if '"' in h else h
            A = _floats(lines[i + 1:i + 1 + n])
            T = _floats(lines[i + 1 + n:i + 1 + n + e]).astype(int) - 1 if e else np.zeros((0, 3), int)
            out.append((title, A, T)); i += 1 + n + e
        else:
            i += 1
    return out


def read_surface(path):
    """Rslt_SurfCal3_EM.dat -> nodes (N, 3), data (N, 23), Q6 elements (M, 6), zone index per node and element.
    Each Q6 element is written as four linear triangles (t1, t2, t3 at the corners, t4 = mid-side nodes)."""
    xs, data, els, prtl, eprtl, off = [], [], [], [], [], 0
    for p, (title, A, T) in enumerate(_zones(path), 1):
        assert len(T) % 4 == 0, 'unexpected connectivity in ' + path
        t1, t2, t3, t4 = T[0::4], T[1::4], T[2::4], T[3::4]
        q6 = np.column_stack([t1[:, 0], t2[:, 1], t3[:, 2], t1[:, 1], t2[:, 2], t1[:, 2]])
        assert np.array_equal(t4, q6[:, 3:6]), 'unexpected Q6 sub-triangles in ' + path
        data.append(A); els.append(q6 + off)
        prtl.append(np.full(len(A), p)); eprtl.append(np.full(len(q6), p)); off += len(A)
    D = np.vstack(data)
    return dict(x=D[:, 0:3], data=D, el=np.vstack(els), prtl=np.concatenate(prtl), eprtl=np.concatenate(eprtl))


def read_incident(run):
    """k, polarisation and direction of the plane wave in Input_Phys_EM.dat (record "pwe")"""
    for line in open(os.path.join(run, 'Input_Phys_EM.dat')):
        if line.startswith('"pwe"'):
            f = [t.strip().strip('"') for t in line.split(',')]
            num = lambda t: float(t.lower().replace('d', 'e'))
            assert f[4].lower() == 'k', 'the wavenumber must be given with "k"'
            k = num(f[5]); pol = np.array([num(t) for t in f[6:9]]); d = np.array([num(t) for t in f[9:12]])
            return k, pol / np.linalg.norm(pol), d / np.linalg.norm(d), (int(f[2]), num(f[3]))
    raise ValueError('no "pwe" record in Input_Phys_EM.dat')


def check_outer_index(run):
    lines = [l for l in open(os.path.join(run, 'Input_Phys_EM.dat')) if l.strip()]
    for i, l in enumerate(lines):
        if l.startswith('!the outer domain medium constants'):
            n = float(lines[i + 1].split(',')[0].lower().replace('d', 'e'))
            assert n == 1.0, 'only an outer medium of refractive index 1 is supported'
            return


# ----------------------------------------------------------------------------------------------
# curved six-node element of the solver (Geom_NormVec.f90): mid-side nodes at parameters alpha, beta, gamma
def q6_params(X):
    d = lambda i, j: np.linalg.norm(X[:, i] - X[:, j], axis=1)
    return (1.0 / (1.0 + d(3, 1) / d(3, 0)), 1.0 / (1.0 + d(5, 2) / d(5, 0)), 1.0 / (1.0 + d(4, 1) / d(4, 2)))


def q6_shape(xi, eta, a, b, g):
    """xi, eta (M, P); a, b, g (M, 1) -> N, dN/dxi, dN/deta (M, P, 6)"""
    N2 = xi * (xi - a + (a - g) / (1 - g) * eta) / (1 - a)
    N3 = eta * (eta - b + (b + g - 1) / g * xi) / (1 - b)
    N4 = xi * (1 - xi - eta) / (a * (1 - a))
    N5 = xi * eta / (g * (1 - g))
    N6 = eta * (1 - xi - eta) / (b * (1 - b))
    N1 = 1 - N2 - N3 - N4 - N5 - N6
    dx2 = (xi - a + (a - g) / (1 - g) * eta + xi) / (1 - a)
    dx3 = eta * (b + g - 1) / ((1 - b) * g)
    dx4 = (1 - 2 * xi - eta) / (a * (1 - a))
    dx5 = eta / (g * (1 - g))
    dx6 = -eta / (b * (1 - b))
    dx1 = -(dx2 + dx3 + dx4 + dx5 + dx6)
    de2 = xi * (a - g) / ((1 - a) * (1 - g))
    de3 = (eta - b + (b + g - 1) / g * xi + eta) / (1 - b)
    de4 = -xi / (a * (1 - a))
    de5 = xi / (g * (1 - g))
    de6 = (1 - xi - 2 * eta) / (b * (1 - b))
    de1 = -(de2 + de3 + de4 + de5 + de6)
    return (np.stack([N1, N2, N3, N4, N5, N6], -1), np.stack([dx1, dx2, dx3, dx4, dx5, dx6], -1),
            np.stack([de1, de2, de3, de4, de5, de6], -1))


class Surface:
    """Q6 surface with complex nodal H; J and rho at the 16-point rule"""

    def __init__(self, surf, H, k):
        self.x, self.el, self.prtl, self.eprtl = surf['x'], surf['el'], surf['prtl'], surf['eprtl']
        self.k, self.H = k, H
        self.X = self.x[self.el]
        a, b, g = q6_params(self.X)
        self.abg = (a[:, None], b[:, None], g[:, None])
        # outward orientation of each element, from the centroid of its conductor's nodes
        # (valid for conductors that are star-shaped about that point, such as the rods of Fig. 3)
        c = np.array([self.x[self.prtl == p].mean(0) for p in range(1, self.prtl.max() + 1)])
        N, Nx, Ne = q6_shape(np.full((len(a), 1), 1 / 3), np.full((len(a), 1), 1 / 3), *self.abg)
        rc = np.einsum('mpl,mlc->mpc', N, self.X)[:, 0]
        A = np.cross(np.einsum('mpl,mlc->mpc', Nx, self.X)[:, 0], np.einsum('mpl,mlc->mpc', Ne, self.X)[:, 0])
        self.s = np.sign(np.sum(A * (rc - c[self.eprtl - 1]), 1))
        assert np.all(self.s != 0)
        self.centroid = rc
        self.diam = np.max(np.linalg.norm(self.X - rc[:, None, :], axis=2), axis=1) * 2
        self.quad = self.eval_points(np.broadcast_to(TRI16_XI, (len(a), 16)),
                                     np.broadcast_to(TRI16_ETA, (len(a), 16)), np.broadcast_to(TRI16_W, (len(a), 16)))
        for p in range(1, self.prtl.max() + 1):     # outward check: volume = (1/3) oint r.n dS > 0
            m = self.eprtl == p
            vol = np.sum(self.quad['wA'][m] * np.sum((self.quad['r'][m] - c[p - 1]) * self.quad['n'][m], 2)) / 3
            assert vol > 0, 'inconsistent normal orientation on conductor %d' % p

    def eval_points(self, xi, eta, w, elems=None, geom_only=False):
        idx = np.arange(len(self.el)) if elems is None else elems
        a, b, g = (t[idx] for t in self.abg)
        X = self.X[idx]
        N, Nx, Ne = q6_shape(xi, eta, a, b, g)
        r = np.einsum('mpl,mlc->mpc', N, X)
        rx = np.einsum('mpl,mlc->mpc', Nx, X); re_ = np.einsum('mpl,mlc->mpc', Ne, X)
        A = np.cross(rx, re_) * self.s[idx][:, None, None]
        jac = np.linalg.norm(A, axis=2)
        out = dict(r=r, n=A / jac[..., None], wA=w * jac)
        if not geom_only:
            He = self.H[self.el[idx]]
            Hq = np.einsum('mpl,mlc->mpc', N, He)
            Hx = np.einsum('mpl,mlc->mpc', Nx, He); Hy = np.einsum('mpl,mlc->mpc', Ne, He)
            out['J'] = np.cross(out['n'], Hq)
            out['rho'] = self.s[idx][:, None] * (np.sum(Hy * rx, 2) - np.sum(Hx * re_, 2)) / jac / (1j * self.k)
        return out

    def moments(self):
        q, res = self.quad, {}
        for p in range(1, self.prtl.max() + 1):
            m = self.eprtl == p
            w, rho = q['wA'][m], q['rho'][m]
            res[p] = dict(Q=np.sum(w * rho), Qabs=np.sum(w * np.abs(rho)),
                          p=np.einsum('mp,mp,mpc->c', w, rho, q['r'][m]),
                          m=0.5 * np.einsum('mp,mpc->c', w, np.cross(q['r'][m], q['J'][m])))
        return res

    def farfield(self, dirs, chunk=64):
        q = self.quad
        Y = q['r'].reshape(-1, 3); Jw = (q['J'] * q['wA'][..., None]).reshape(-1, 3)
        F = np.empty((len(dirs), 3), complex)
        for s in range(0, len(dirs), chunk):
            d = dirs[s:s + chunk]
            Jt = np.exp(-1j * self.k * (d @ Y.T)) @ Jw
            Jt -= d * np.sum(d * Jt, 1)[:, None]
            F[s:s + chunk] = 1j * self.k / (4 * np.pi) * Jt
        return F

    @staticmethod
    def _contrib(D, Jw, qw, k):
        """E_sc from point sources: D = target - source (T, S, 3), current moments Jw (T, S, 3), charges qw (T, S)"""
        R = np.linalg.norm(D, axis=-1)
        e = np.exp(1j * k * R) / (4 * np.pi * R)
        return 1j * k * np.einsum('ps,psc->pc', e, Jw) - np.einsum('ps,psc->pc', e * (1j * k * R - 1) / R ** 2 * qw, D)

    def field(self, P, near_c=4.0, sub_c=3.0, lmax=10, chunk=64):
        """complex scattered E at exterior targets P (T, 3)"""
        q, k = self.quad, self.k
        Y = q['r'].reshape(-1, 3)
        Jw = (q['J'] * q['wA'][..., None]).reshape(-1, 3); qw = (q['rho'] * q['wA']).reshape(-1)
        E = np.zeros((len(P), 3), complex)
        for s in range(0, len(P), chunk):
            p = P[s:s + chunk]
            E[s:s + chunk] = self._contrib(p[:, None, :] - Y[None], np.broadcast_to(Jw, (len(p),) + Jw.shape),
                                           np.broadcast_to(qw, (len(p), len(qw))), k)
        # near target-element pairs: replace the 16-point contribution by an adaptive one
        lists = cKDTree(self.centroid).query_ball_point(P, r=near_c * self.diam.max())
        ti, ei = [], []
        for t, L in enumerate(lists):
            if L:
                L = np.array(L)
                L = L[np.linalg.norm(self.centroid[L] - P[t], axis=1) < near_c * self.diam[L]]
                ti.append(np.full(len(L), t)); ei.append(L)
        if not ti:
            return E
        ti = np.concatenate(ti); ei = np.concatenate(ei)
        np.add.at(E, ti, -self._contrib(P[ti][:, None, :] - q['r'][ei], q['J'][ei] * q['wA'][ei][..., None],
                                        q['rho'][ei] * q['wA'][ei], k))
        V = np.tile(np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]), (len(ti), 1, 1))
        pid = np.arange(len(ti))
        for lev in range(lmax + 1):
            vc = V.mean(1)
            gv = self.eval_points(np.concatenate([V[:, :, 0], vc[:, None, 0]], 1),
                                  np.concatenate([V[:, :, 1], vc[:, None, 1]], 1), np.zeros((len(pid), 4)),
                                  elems=ei[pid], geom_only=True)
            rv = gv['r']
            size = np.max(np.linalg.norm(rv[:, :3] - rv[:, 3:4], axis=2), 1) * 2
            acc = (np.linalg.norm(P[ti[pid]] - rv[:, 3], axis=1) > sub_c * size) | (lev == lmax)
            if acc.any():
                Va, pa = V[acc], pid[acc]
                e1 = Va[:, 1] - Va[:, 0]; e2 = Va[:, 2] - Va[:, 0]
                det = np.abs(e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0])
                xi = Va[:, 0, 0][:, None] + TRI16_XI[None] * e1[:, 0][:, None] + TRI16_ETA[None] * e2[:, 0][:, None]
                eta = Va[:, 0, 1][:, None] + TRI16_XI[None] * e1[:, 1][:, None] + TRI16_ETA[None] * e2[:, 1][:, None]
                w = det[:, None] * TRI16_W[None]
                for s in range(0, len(pa), 20000):
                    sl = slice(s, s + 20000)
                    g = self.eval_points(xi[sl], eta[sl], w[sl], elems=ei[pa[sl]])
                    np.add.at(E, ti[pa[sl]], self._contrib(P[ti[pa[sl]]][:, None, :] - g['r'],
                                                           g['J'] * g['wA'][..., None], g['rho'] * g['wA'], k))
            if (~acc).any():
                Vr, pr = V[~acc], pid[~acc]
                m01 = 0.5 * (Vr[:, 0] + Vr[:, 1]); m12 = 0.5 * (Vr[:, 1] + Vr[:, 2]); m20 = 0.5 * (Vr[:, 2] + Vr[:, 0])
                V = np.concatenate([np.stack([Vr[:, 0], m01, m20], 1), np.stack([m01, Vr[:, 1], m12], 1),
                                    np.stack([m20, m12, Vr[:, 2]], 1), np.stack([m01, m12, m20], 1)])
                pid = np.concatenate([pr, pr, pr, pr])
            else:
                break
        return E


def rcs_dirs(cut, n=361):
    t = np.linspace(0, 2 * np.pi, n); c, s, z = np.cos(t), np.sin(t), np.zeros(n)
    return np.rad2deg(t), {'zx': np.column_stack([c, z, s]), 'xy': np.column_stack([c, s, z]),
                           'yz': np.column_stack([z, c, s])}[cut]


# ----------------------------------------------------------------------------------------------
def prepare(case, run90):
    if os.path.exists(run90):
        sys.exit('%s exists; give a new directory' % run90)
    shutil.copytree(case, run90)
    p = os.path.join(run90, 'Input_Phys_EM.dat')
    lines = open(p).read().splitlines(keepends=True)
    n_pwe = n_slice = 0
    for i, l in enumerate(lines):
        if l.startswith(PWE_PHASE0):
            lines[i] = l.replace(PWE_PHASE0, PWE_PHASE90, 1); n_pwe += 1
        elif l.startswith('!If turn on field slice calculation'):
            j = i + 1
            if lines[j].startswith('1,'):
                lines[j] = '0,' + lines[j][2:]; n_slice += 1
    if n_pwe != 1:
        sys.exit('expected one record starting with %s in %s' % (PWE_PHASE0, p))
    open(p, 'w').write(''.join(lines))
    print('%s: incident field multiplied by i (%s)%s' % (run90, PWE_PHASE90, '; planar samples off' if n_slice else ''))


def compute(run0, run90, out, slice_=True):
    os.makedirs(out, exist_ok=True)
    check_outer_index(run0)
    k, pol, kdir, (order0, phase0) = read_incident(run0)
    k9, pol9, kdir9, (order9, phase9) = read_incident(run90)
    assert (order0, phase0) == (0, 0.0) and (order9, abs(phase9 - np.pi / 2) < 1e-12) == (3, True), \
        'RUN0 must have incident phase 0 and RUN90 incident phase pi/2 (use "prepare")'
    assert k9 == k and np.allclose(pol9, pol) and np.allclose(kdir9, kdir), 'the two runs differ'
    s0 = read_surface(os.path.join(run0, 'Rslt_SurfCal3_EM.dat'))
    s9 = read_surface(os.path.join(run90, 'Rslt_SurfCal3_EM.dat'))
    assert np.array_equal(s0['x'], s9['x']) and np.array_equal(s0['el'], s9['el']), 'the two runs have different meshes'
    H = s0['data'][:, 13:16] - 1j * s9['data'][:, 13:16]          # columns exH3x, exH3y, exH3z
    S = Surface(s0, H, k)
    lines = ['k = %g (k_out a for a = 1), %d nodes, %d elements, %d conductors'
             % (k, len(S.x), len(S.el), S.prtl.max())]
    mom = S.moments()
    ptot = sum(m['p'] for m in mom.values()); mtot = sum(m['m'] for m in mom.values())
    for p, m in mom.items():
        lines.append('conductor %d: |Q| / oint|rho| dS = %.1e' % (p, abs(m['Q']) / m['Qabs']))
    f = lambda v: '(' + ', '.join('%.6g%+.3gi' % (c.real, c.imag) for c in v) + ')'
    lines.append('electric dipole moment p = int r rho dS = ' + f(ptot))
    lines.append('magnetic dipole moment m = (1/2) int r x J dS = ' + f(mtot))
    for cut in ('zx', 'xy', 'yz'):
        th, d = rcs_dirs(cut)
        sig = 4 * np.pi * np.sum(np.abs(S.farfield(d)) ** 2, 1)
        np.savetxt(os.path.join(out, 'rcs_far_%s.csv' % cut), np.column_stack([th, sig, 10 * np.log10(sig)]),
                   delimiter=',', fmt=['%.1f', '%.10e', '%.6f'], header='theta_deg,sigma_over_a2,dB', comments='')
        i = int(np.argmax(sig))
        lines.append('RCS %s (far field): max %.4e a^2 (%.2f dB) at %.0f deg; at 90 deg %.4e a^2, at 270 deg %.4e a^2'
                     % (cut, sig[i], 10 * np.log10(sig[i]), th[i], sig[90], sig[270]))
    dmn = os.path.join(run0, 'Rslt_Dmn3_plot.dat')
    if slice_ and os.path.isfile(dmn):
        pos = {t.replace('pos_', '').strip(): A for t, A, _ in _zones(os.path.join(run0, 'RsltPos.dat'))}
        for title, A, _ in _zones(dmn):
            if 'plane' not in title:
                continue
            plane = title.split('plane')[-1].strip(' -')
            ps = pos['plane - ' + plane]
            assert np.allclose(ps[:, 0:3], A[:, 0:3]), 'RsltPos.dat does not match ' + dmn
            P = A[:, 0:3]; ext = ps[:, 3] == 0
            E = np.zeros((len(P), 3), complex)
            E[ext] = S.field(P[ext]) + pol[None, :] * np.exp(1j * k * P[ext] @ kdir)[:, None]
            np.savetxt(os.path.join(out, 'slice_%s.csv' % plane),
                       np.column_stack([P, ps[:, 3], E.real]), delimiter=',',
                       fmt=['%.8f'] * 3 + ['%d'] + ['%.10e'] * 3,
                       header='x,y,z,pos,Re_Ex_tot,Re_Ey_tot,Re_Ez_tot', comments='')
            lines.append('planar samples %s: %d points (%d outside the conductors), Re E_x^tot from %.4f to %.4f'
                         % (plane, len(P), ext.sum(), E[ext, 0].real.min(), E[ext, 0].real.max()))
    open(os.path.join(out, 'summary.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    a = sys.argv[1:]
    if len(a) == 3 and a[0] == 'prepare':
        prepare(a[1], a[2])
    elif len(a) in (4, 5) and a[0] == 'compute':
        compute(a[1], a[2], a[3], slice_='--no-slice' not in a)
    else:
        sys.exit(__doc__)
