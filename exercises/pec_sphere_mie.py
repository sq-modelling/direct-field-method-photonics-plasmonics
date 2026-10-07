#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Qiang Sun
# SPDX-License-Identifier: BSD-3-Clause
"""Exercise 1: compare the DFM solution for a PEC sphere with Mie theory.

This script reproduces the two numerical checks of Fig. 2 of the manuscript
from the output of ``code/dfm_field_solver``.  It needs only NumPy and SciPy;
Matplotlib is used only when ``--plot`` is given.

Sub-commands
------------
surface RUN
    Read the total surface field ``Rslt_SurfCal3_EM.dat`` of a PEC-sphere run
    (RUN is the run directory or the file itself), evaluate the exact Mie
    solution at the same nodes and print the relative L2 error of |E_tot|
    over all nodes, as in Fig. 2(a).  On a PEC surface the total field is
    purely normal, so |E_tot| = |E_n|.

point RUN [RUN ...]
    For one or more runs of ``exercises/pec_sphere_convergence_ka4p5`` (one
    run per mesh level), read |E_sc| at an observation point (default
    (0, 0, 1.2a)) from the planar samples in ``Rslt_Dmn2_plot.dat``, compare
    it with Mie theory and fit the convergence rate, as in Fig. 2(b).

selftest
    Check the Mie implementation itself: the series for the incident wave
    must reproduce exp(i k z) x, and the tangential total field must vanish
    on the sphere.

k_out a, the polarisation and the propagation direction are read from
``Input_Phys_EM.dat`` in the run directory, and the centre and radius of the
sphere from the surface nodes.  For a result file without its input files,
give ``--ka`` (and, if they differ from x and z, ``--pol`` and ``--dir``).

Conventions are those of the solver: time dependence exp(-i omega t); incident
field E_inc = E0 p exp(i k_out khat.r), where the solver normalises the
polarisation vector p to unit length, so that all fields are in units of E0
(code/dfm_field_solver/EM_SurfCal_Input.f90).  The
Mie series follows Bohren and Huffman, "Absorption and Scattering of Light by
Small Particles" (Wiley, 1983), ch. 4, in the perfect-conductor limit
a_n = psi_n'(x)/xi_n'(x), b_n = psi_n(x)/xi_n(x), x = k_out a.

Examples (from the repository root)
-----------------------------------
    python3 exercises/pec_sphere_mie.py surface ../dfm-runs/sphere
    python3 exercises/pec_sphere_mie.py surface ../dfm-runs/sphere --plot eplane.png
    python3 exercises/pec_sphere_mie.py point ../dfm-runs/sphere-convergence/level*

For the sphere case of Fig. 1(a) (inputs/fig01_pec_sphere_cube/sphere,
k_out a = 15, 2562 nodes), ``surface`` reports a relative L2 error of 1.15 %
and a largest nodal error of 0.0365 E0.
"""

import argparse
import os
import re
import sys

import numpy as np
from scipy.special import spherical_jn, spherical_yn

SURFACE_FILE = "Rslt_SurfCal3_EM.dat"
SLICE_FILE = "Rslt_Dmn2_plot.dat"
PHYS_FILE = "Input_Phys_EM.dat"


# ---------------------------------------------------------------------------
# Reading the solver files
# ---------------------------------------------------------------------------

def _float(token):
    """Convert a Fortran real, e.g. 1.5d0 or 0.12-100 (G format), to float."""
    t = token.strip().strip("'\"").replace("d", "e").replace("D", "E")
    try:
        return float(t)
    except ValueError:
        m = re.fullmatch(r"([+-]?\d*\.\d*)([+-]\d+)", t)
        if m:
            return float(m.group(1) + "E" + m.group(2))
        raise


def read_zones(path):
    """Return the node data of every Tecplot FEPOINT zone in a solver file.

    Each item is (zone_title, array of shape (n_nodes, n_columns)); the
    element connectivity that follows the nodes is skipped.
    """
    zones = []
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.lstrip().lower().startswith("zone"):
            title = re.search(r'T\s*=\s*"([^"]*)"', line)
            nn = re.search(r"\bn\s*=\s*(\d+)", line)
            if nn is None:
                raise ValueError("cannot read the node count in %s: %s" % (path, line))
            n = int(nn.group(1))
            rows = [[_float(v) for v in re.split(r"[,\s]+", lines[i + 1 + k].strip()) if v]
                    for k in range(n)]
            zones.append((title.group(1).strip() if title else "", np.array(rows)))
            i += n + 1
        else:
            i += 1
    if not zones:
        raise ValueError("no Tecplot zone found in %s" % path)
    return zones


def read_phys_input(path):
    """Read the incident plane wave and boundary types from Input_Phys_EM.dat.

    The file is record oriented (see code/dfm_field_solver/EM_SurfCal_Input.f90):
    record 11 holds the excitation, record 14 the exterior medium, and record
    21 + 10(p - 1) the boundary type of surface p.
    """
    with open(path) as f:
        rec = f.read().splitlines()

    def tokens(line):
        return [t.strip().strip("'\"") for t in re.split(r"[,\s]+", line.strip()) if t.strip()]

    t = tokens(rec[10])
    exc, order, feature, mode, value = t[0], int(t[2]), _float(t[3]), t[4].lower(), _float(t[5])
    pol = np.array([_float(v) for v in t[6:9]])
    kdir = np.array([_float(v) for v in t[9:12]])
    if exc.lower() != "pwe":
        raise ValueError("only plane-wave excitation 'pwe' is supported, found %r" % exc)
    if mode == "k":
        k0 = value
    elif mode == "l":
        k0 = 2.0 * np.pi / value
    else:
        raise ValueError("wavenumber given as %r; use --ka for this case" % mode)
    m = [_float(v) for v in tokens(rec[13])[:4]]
    eps, mu = complex(m[0], m[1]) ** 2, complex(m[2], m[3]) ** 2
    k_out = k0 * np.sqrt(eps * mu)
    if abs(k_out.imag) > 1e-12 * abs(k_out):
        raise ValueError("the exterior medium is lossy; this exercise assumes a lossless host")
    # The solver applies the phase 'feature' to y, z or x for order 1, 2 or 3 (0 -> 1).
    phased = {1: 1, 2: 2, 3: 0}.get(order if order in (1, 2, 3) else 1)
    if abs(feature) > 1e-12 and abs(pol[phased]) > 1e-12:
        raise ValueError("elliptic polarisation is not supported by this exercise")
    bcs = []
    p = 0
    while 20 + 10 * p < len(rec) and rec[18 + 10 * p].lstrip().startswith("~"):
        bcs.append(tokens(rec[20 + 10 * p])[0].upper())
        p += 1
    return dict(k_out=k_out.real, pol=pol, kdir=kdir, bcs=bcs)


def find_file(run, name):
    """Return (path of the result file, directory searched for the input files)."""
    if os.path.isdir(run):
        return os.path.join(run, name), run
    return run, os.path.dirname(os.path.abspath(run))


def incident_wave(args, run_dir, radius):
    """k_out a, polarisation and direction, from the arguments or the input file."""
    if args.ka is not None:
        return args.ka, np.array(args.pol, float), np.array(args.dir, float), "command line"
    path = os.path.join(run_dir, PHYS_FILE) if run_dir else None
    if not path or not os.path.isfile(path):
        sys.exit("error: %s not found; give k_out a with --ka" % (path or PHYS_FILE))
    ph = read_phys_input(path)
    if len(ph["bcs"]) != 1 or ph["bcs"][0] != "PEC":
        print("warning: boundary types %s; this exercise is for a single PEC sphere" % ph["bcs"])
    return ph["k_out"] * radius, ph["pol"], ph["kdir"], path


def local_frame(pol, kdir):
    """Orthonormal frame (e1, e2, e3) with e1 along E_inc and e3 along k_inc."""
    e1 = np.asarray(pol, float) / np.linalg.norm(pol)
    e3 = np.asarray(kdir, float) / np.linalg.norm(kdir)
    if abs(e1 @ e3) > 1e-8:
        sys.exit("error: the polarisation is not perpendicular to the propagation direction")
    return np.array([e1, np.cross(e3, e1), e3])


def sphere_fit(xyz):
    """Centre and radius of the sphere from its surface nodes."""
    c = xyz.mean(axis=0)
    r = np.linalg.norm(xyz - c, axis=1)
    return c, r.mean(), (r.max() - r.min()) / r.mean()


# ---------------------------------------------------------------------------
# Mie solution for a perfectly conducting sphere (Bohren & Huffman, ch. 4)
# ---------------------------------------------------------------------------

def mie_pec_field(x, xyz, part="total"):
    """Complex electric field of a unit plane wave on a PEC sphere.

    x    : size parameter k_out a.
    xyz  : (P, 3) points in the local frame, in units of a, with |xyz| >= 1;
           the incident field is (1, 0, 0) exp(i x z).
    part : 'incident', 'scattered' or 'total'.
    Returns the (P, 3) Cartesian field in the same frame.
    """
    xyz = np.atleast_2d(xyz)
    r = np.linalg.norm(xyz, axis=1)
    theta = np.arccos(np.clip(xyz[:, 2] / r, -1.0, 1.0))
    phi = np.arctan2(xyz[:, 1], xyz[:, 0])
    rho = x * r

    nmax = int(np.ceil(x * r.max() + 4.05 * (x * r.max()) ** (1.0 / 3.0) + 10))
    n = np.arange(1, nmax + 1)[:, None]
    En = 1j ** n * (2 * n + 1) / (n * (n + 1))

    def psi(z, d=False):
        return spherical_jn(n, z) + z * spherical_jn(n, z, True) if d else z * spherical_jn(n, z)

    def xi(z, d=False):
        h, dh = spherical_jn(n, z) + 1j * spherical_yn(n, z), \
            spherical_jn(n, z, True) + 1j * spherical_yn(n, z, True)
        return h + z * dh if d else z * h

    an = psi(x, True) / xi(x, True)
    bn = psi(x) / xi(x)

    mu = np.cos(theta)                       # angular functions pi_n, tau_n
    pi = np.zeros((nmax, mu.size))
    pi[0] = 1.0
    if nmax > 1:
        pi[1] = 3.0 * mu
    for k in range(2, nmax):
        m = k + 1
        pi[k] = (2 * m - 1) / (m - 1) * mu * pi[k - 1] - m / (m - 1) * pi[k - 2]
    pim1 = np.vstack([np.zeros_like(mu), pi[:-1]])
    tau = n * mu * pi - (n + 1) * pim1

    j, dj = spherical_jn(n, rho), spherical_jn(n, rho, True)
    h = j + 1j * spherical_yn(n, rho)
    dh = dj + 1j * spherical_yn(n, rho, True)
    Dj, Dh = j / rho + dj, h / rho + dh      # (rho z_n)'/rho
    nn1 = n * (n + 1)

    er = np.zeros(mu.size, complex)
    et = np.zeros(mu.size, complex)
    ep = np.zeros(mu.size, complex)
    if part in ("incident", "total"):       # E_i = sum E_n (M_o1n - i N_e1n), z_n = j_n
        er += np.sum(En * (-1j) * nn1 * pi * j / rho, axis=0)
        et += np.sum(En * (pi * j - 1j * tau * Dj), axis=0)
        ep += np.sum(En * (-tau * j + 1j * pi * Dj), axis=0)
    if part in ("scattered", "total"):      # E_s = sum E_n (i a_n N_e1n - b_n M_o1n), z_n = h_n
        er += np.sum(En * 1j * an * nn1 * pi * h / rho, axis=0)
        et += np.sum(En * (1j * an * tau * Dh - bn * pi * h), axis=0)
        ep += np.sum(En * (-1j * an * pi * Dh + bn * tau * h), axis=0)
    cp, sp, ct, st = np.cos(phi), np.sin(phi), np.cos(theta), np.sin(theta)
    er, et, ep = er * cp * st, et * cp, ep * sp
    return np.stack([er * st * cp + et * ct * cp - ep * sp,
                     er * st * sp + et * ct * sp + ep * cp,
                     er * ct - et * st], axis=1)


# ---------------------------------------------------------------------------
# Sub-commands
# ---------------------------------------------------------------------------

def cmd_surface(args):
    path, run_dir = find_file(args.run, SURFACE_FILE)
    if not os.path.isfile(path):
        sys.exit("error: %s not found" % path)
    zones = read_zones(path)
    if len(zones) > 1:
        print("warning: %s holds %d zones (appended runs?); using the last one" % (path, len(zones)))
    data = zones[-1][1]
    xyz, re_e, mag = data[:, 0:3], data[:, 3:6], data[:, 6]
    centre, a, spread = sphere_fit(xyz)
    ka, pol, kdir, source = incident_wave(args, run_dir, a)
    R = local_frame(pol, kdir)
    loc = (xyz - centre) @ R.T / a
    loc /= np.linalg.norm(loc, axis=1)[:, None]          # evaluate exactly on r = a

    e_mie = mie_pec_field(ka, loc) @ R                   # back to the solver frame
    mag_mie = np.linalg.norm(np.abs(e_mie), axis=1)
    err = mag - mag_mie
    rel_l2 = np.linalg.norm(err) / np.linalg.norm(mag_mie)
    rel_l2_re = np.linalg.norm(re_e - e_mie.real) / np.linalg.norm(e_mie.real)
    eplane = np.abs(loc[:, 1]) < args.band

    print("DFM result      : %s (%d nodes)" % (path, len(mag)))
    print("sphere          : centre (%.3g, %.3g, %.3g), radius a = %.6g (node radii spread %.1e a)"
          % (*centre, a, spread))
    print("incident wave   : k_out a = %.6g, polarisation %s, direction %s  [%s]"
          % (ka, np.round(R[0], 6), np.round(R[2], 6), source))
    print("|E_tot|/E0 on S : relative L2 error over all %d nodes = %.2f %%" % (len(mag), 100 * rel_l2))
    print("                  largest nodal error = %.4f E0 (max |E_tot| = %.4f E0)"
          % (np.abs(err).max(), mag_mie.max()))
    print("Re(E_tot)       : relative L2 error of the real-part vectors = %.2f %%" % (100 * rel_l2_re))
    print("E-plane nodes   : %d with |y'| < %.3g a" % (eplane.sum(), args.band))

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        t = np.linspace(0.0, np.pi, 721)
        curve = [np.linalg.norm(np.abs(mie_pec_field(ka, np.c_[s * np.sin(t), 0 * t, np.cos(t)])), axis=1)
                 for s in (1.0, -1.0)]
        ang = np.degrees(np.arccos(np.clip(loc[eplane, 2], -1, 1))) * np.where(loc[eplane, 0] >= 0, 1, -1)
        fig, ax = plt.subplots(figsize=(5.5, 3.2))
        ax.plot(np.degrees(t), curve[0], "k-", lw=1, label="Mie")
        ax.plot(-np.degrees(t), curve[1], "k-", lw=1)
        ax.plot(ang, mag[eplane], "o", mfc="none", ms=4, label="DFM, %d nodes" % len(mag))
        ax.set_xlim(-180, 180)
        ax.set_xlabel(r"polar angle $\theta$ in the $E$-plane (deg)")
        ax.set_ylabel(r"$|E^{\mathrm{tot}}|/E_0$ on $S$")
        ax.set_title(r"$k_{\mathrm{out}}a = %g$, relative $L_2$ error %.2f%%" % (ka, 100 * rel_l2))
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(args.plot, dpi=200)
        print("plot written to : %s" % args.plot)


def cmd_point(args):
    target = np.array(args.xyz, float)
    rows = []
    for run in args.runs:
        missing = [f for f in (SURFACE_FILE, SLICE_FILE) if not os.path.isfile(os.path.join(run, f))]
        if missing:
            print("warning: skipping %s (no %s; unfinished or failed run?)" % (run, " or ".join(missing)))
            continue
        surf = read_zones(os.path.join(run, SURFACE_FILE))[-1][1]
        centre, a, _ = sphere_fit(surf[:, 0:3])
        ka, pol, kdir, _ = incident_wave(args, run, a)
        R = local_frame(pol, kdir)
        planes = [z for t, z in read_zones(os.path.join(run, SLICE_FILE)) if "plane" in t]
        if not planes:
            sys.exit("error: no planar samples in %s; enable the slice in %s"
                     % (os.path.join(run, SLICE_FILE), PHYS_FILE))
        pts = np.vstack(planes)
        want = centre + a * (target @ R)                 # target is given in the local frame
        d = np.linalg.norm(pts[:, 0:3] - want, axis=1)
        k = int(np.argmin(d))
        if d[k] > 1e-6 * a:
            sys.exit("error: %s has no sample at %s (nearest %s)" % (run, want, pts[k, 0:3]))
        e_dfm = pts[k, 6]                                # |E2| = |E_sc|
        e_mie = np.linalg.norm(np.abs(mie_pec_field(ka, target[None, :], "scattered")))
        elapsed = None
        tfile = os.path.join(run, "elapsed_s.txt")
        if os.path.isfile(tfile):
            with open(tfile) as f:
                elapsed = f.read().strip()
        rows.append((run, len(surf), ka, e_dfm, e_mie, abs(e_dfm - e_mie) / e_mie, elapsed))

    print("observation point (local frame, units of a): (%g, %g, %g)" % tuple(target))
    print("%-34s %7s %8s %16s %16s %10s %8s" % ("run", "N", "k_out a", "|E_sc| DFM", "|E_sc| Mie",
                                              "rel. err.", "wall s"))
    for run, n_nodes, ka, ed, em, rel, el in sorted(rows, key=lambda r: r[1]):
        print("%-34s %7d %8.4g %16.11f %16.11f %10.3e %8s"
              % (os.path.basename(os.path.normpath(run)), n_nodes, ka, ed, em, rel, el or "-"))
    if len(rows) >= 2:
        N = np.array([r[1] for r in rows], float)
        err = np.array([r[5] for r in rows])
        slope = np.polyfit(np.log(N), np.log(err), 1)[0]
        print("least-squares slope of log(error) against log(N): %.2f  (N^-2, i.e. h^4, gives -2)" % slope)


def cmd_selftest(args):
    rng = np.random.default_rng(1)
    for x in (0.5, 4.5, 15.0):
        v = rng.normal(size=(400, 3))
        v /= np.linalg.norm(v, axis=1)[:, None]
        p = v * rng.uniform(1.0, 3.0, size=(400, 1))
        e_inc = mie_pec_field(x, p, "incident")
        exact = np.zeros_like(e_inc)
        exact[:, 0] = np.exp(1j * x * p[:, 2])
        e_tot = mie_pec_field(x, v, "total")
        tang = e_tot - np.sum(e_tot * v, axis=1)[:, None] * v
        print("k a = %5.1f: max |E_inc(series) - exp(ikz) x| = %.1e;  max |E_tan| on S = %.1e"
              % (x, np.abs(e_inc - exact).max(), np.linalg.norm(tang, axis=1).max()))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--ka", type=float, help="k_out a; overrides Input_Phys_EM.dat")
    common.add_argument("--pol", type=float, nargs=3, default=(1, 0, 0),
                        help="incident polarisation, used with --ka (default 1 0 0)")
    common.add_argument("--dir", type=float, nargs=3, default=(0, 0, 1),
                        help="propagation direction, used with --ka (default 0 0 1)")
    s = sub.add_parser("surface", parents=[common], help="surface field against Mie (Fig. 2(a))")
    s.add_argument("run", help="run directory or %s" % SURFACE_FILE)
    s.add_argument("--band", type=float, default=0.03,
                   help="half-width of the E-plane band, in units of a (default 0.03)")
    s.add_argument("--plot", metavar="PNG", help="also plot |E_tot| in the E-plane (needs Matplotlib)")
    p = sub.add_parser("point", parents=[common], help="|E_sc| at a point against Mie (Fig. 2(b))")
    p.add_argument("runs", nargs="+", help="run directories, e.g. one per mesh level")
    p.add_argument("--xyz", type=float, nargs=3, default=(0.0, 0.0, 1.2),
                   help="observation point in units of a, in the frame (E_inc, k x E_inc, k_inc); "
                        "default 0 0 1.2")
    sub.add_parser("selftest", help="check the Mie implementation")
    args = ap.parse_args()
    {"surface": cmd_surface, "point": cmd_point, "selftest": cmd_selftest}[args.cmd](args)


if __name__ == "__main__":
    main()
