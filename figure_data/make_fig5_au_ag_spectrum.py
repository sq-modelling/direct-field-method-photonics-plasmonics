#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Qiang Sun
# SPDX-License-Identifier: BSD-3-Clause
"""Manuscript Fig. 5: absorption spectra of eccentric Au-core/Ag-shell spheres and of
homogeneous Ag and Au spheres in water.

(a) Lines: Direct Field Method, replotted unchanged from the archived run
    fig5_Rslt_SPR_XS_archived.dat (Ag shell of outer radius a_Ag = 50 nm; Au core of
    radius a_Au = 10, 15, ..., 45 nm, displaced along the polarization so that the
    thinnest part of the Ag shell is h = 5 nm; a_Au = 45 nm is therefore concentric).
    Open circles: coated-sphere Mie theory (Aden-Kerker; Bohren-Huffman 1983, Sec. 8.1)
    for the concentric case a_Au = 45 nm, computed here with the same n, k tables.
(b) Mie theory for homogeneous Ag and Au spheres of radius 50 nm, computed here.

Horizontal axis: free-space wavelength lambda_0.  Host medium: water, n = 1.33.
The Au and Ag optical constants are read from the case inputs of the spectrum
program, ../inputs/fig04_au_ag_spectrum/core_shell/Input_nk_{Ag,Au}.dat, which are
the tables used by the archived run.  Usage, from this directory:

    python3 make_fig5_au_ag_spectrum.py

writes output/Fig_spr_coreshell.pdf, output/Fig_spr_coreshell.png and the
DFM/Mie comparison for a_Au = 45 nm, output/fig5a_mie_check.csv (compare with the
archived copy fig5a_mie_check.csv).  Requires NumPy, SciPy and Matplotlib.
See README.md in this directory for the provenance of the data files.
"""
from pathlib import Path

import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess, glob, os
import matplotlib.font_manager as fm


def _register_termes():
    """Make TeX Gyre Termes available to Matplotlib if a TeX installation provides it."""
    try:
        p = subprocess.run(['kpsewhich', 'texgyretermes-regular.otf'], capture_output=True, text=True,
                           timeout=20).stdout.strip()
    except Exception:
        p = ''
    for f in glob.glob(os.path.join(os.path.dirname(p), 'texgyretermes-*.otf')) if p else []:
        try:
            fm.fontManager.addfont(f)
        except Exception:
            pass


_register_termes()
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from scipy.special import spherical_jn, spherical_yn

HERE = Path(__file__).resolve().parent
NK_DIR = HERE.parent / 'inputs' / 'fig04_au_ag_spectrum' / 'core_shell'
OUT = HERE / 'output'
OUT.mkdir(exist_ok=True)

# ------------------------------------------------------------------ Mie (homogeneous and coated sphere)
def psi(n, z): return z*spherical_jn(n, z)
def dpsi(n, z): return spherical_jn(n, z) + z*spherical_jn(n, z, derivative=True)
def chi(n, z): return -z*spherical_yn(n, z)
def dchi(n, z): return -(spherical_yn(n, z) + z*spherical_yn(n, z, derivative=True))
def xi(n, z): return psi(n, z) - 1j*chi(n, z)
def dxi(n, z): return dpsi(n, z) - 1j*dchi(n, z)

def coated_ab(m1, m2, x, y):
    """Bohren & Huffman (1983) Sec. 8.1: core index m1 (size x), shell m2 (size y), relative indices."""
    n = np.arange(1, int(y + 4*y**(1/3) + 4) + 1)
    An = (m2*psi(n, m2*x)*dpsi(n, m1*x) - m1*dpsi(n, m2*x)*psi(n, m1*x)) / \
         (m2*chi(n, m2*x)*dpsi(n, m1*x) - m1*dchi(n, m2*x)*psi(n, m1*x))
    Bn = (m2*psi(n, m1*x)*dpsi(n, m2*x) - m1*psi(n, m2*x)*dpsi(n, m1*x)) / \
         (m2*dchi(n, m2*x)*psi(n, m1*x) - m1*dpsi(n, m1*x)*chi(n, m2*x))
    a = (psi(n, y)*(dpsi(n, m2*y) - An*dchi(n, m2*y)) - m2*dpsi(n, y)*(psi(n, m2*y) - An*chi(n, m2*y))) / \
        (xi(n, y)*(dpsi(n, m2*y) - An*dchi(n, m2*y)) - m2*dxi(n, y)*(psi(n, m2*y) - An*chi(n, m2*y)))
    b = (m2*psi(n, y)*(dpsi(n, m2*y) - Bn*dchi(n, m2*y)) - dpsi(n, y)*(psi(n, m2*y) - Bn*chi(n, m2*y))) / \
        (m2*xi(n, y)*(dpsi(n, m2*y) - Bn*dchi(n, m2*y)) - dxi(n, y)*(psi(n, m2*y) - Bn*chi(n, m2*y)))
    return n, a, b

def sigma_abs(n, a, b, k):
    ext = 2*np.pi/k**2*np.sum((2*n + 1)*np.real(a + b))
    sca = 2*np.pi/k**2*np.sum((2*n + 1)*(abs(a)**2 + abs(b)**2))
    return (ext - sca)*1e-6                                  # nm^2 -> um^2

def load_nk(f):
    d = np.loadtxt(f, skiprows=1, delimiter=','); return d[:, 0]*1000, d[:, 1] + 1j*d[:, 2]
lAg, nAg = load_nk(NK_DIR/'Input_nk_Ag.dat'); lAu, nAu = load_nk(NK_DIR/'Input_nk_Au.dat')
def nk(l, L, N): return np.interp(l, L, N.real) + 1j*np.interp(l, L, N.imag)
nmed = 1.33

def mie_coated(lam, rc, rs):
    k = 2*np.pi*nmed/lam
    n, a, b = coated_ab(nk(lam, lAu, nAu)/nmed, nk(lam, lAg, nAg)/nmed, k*rc, k*rs)
    return sigma_abs(n, a, b, k)

def mie_sphere(lam, L, N, r):
    # homogeneous sphere = coated sphere with core index equal to shell index
    k = 2*np.pi*nmed/lam; m = nk(lam, L, N)/nmed
    n, a, b = coated_ab(m, m, 0.5*k*r, k*r)
    return sigma_abs(n, a, b, k)

# ------------------------------------------------------------------ archived DFM data
# Zones of the archived file: Au-core radius 5, 10, ..., 45 nm; columns: wavelength (nm),
# incident-flux residual, scattering, extinction and absorption cross sections (um^2).
zones = {}; cur = None
for line in open(HERE/'fig5_Rslt_SPR_XS_archived.dat'):
    s = line.strip()
    if s.startswith('Zone'):
        cur = float(s.split('radius')[1].split('nm')[0]); zones[cur] = []
    elif cur is not None and s[:1].isdigit():
        zones[cur].append([float(v) for v in s.replace(',', ' ').split()])
zones = {k: np.array(v) for k, v in zones.items()}
radii = [10, 15, 20, 25, 30, 35, 40, 45]

lam_m = np.arange(300, 801, 20.0)
mie45 = np.array([mie_coated(l, 45.0, 50.0) for l in lam_m])
d45 = zones[45.0]; dfm45 = np.interp(lam_m, d45[:, 0], d45[:, 4])
rel = (dfm45 - mie45)/mie45
print('a_Au = 45 nm, DFM vs coated Mie: max |rel| 300-600 nm %.2f%%, 620-800 nm %.2f%%' %
      (100*abs(rel[lam_m <= 600]).max(), 100*abs(rel[lam_m > 600]).max()))
lam_f = np.arange(300, 800.01, 2.0)
sAg = np.array([mie_sphere(l, lAg, nAg, 50.0) for l in lam_f])
sAu = np.array([mie_sphere(l, lAu, nAu, 50.0) for l in lam_f])
print('Mie Ag peak %.4f um^2 at %d nm; Au peak %.4f um^2 at %d nm' %
      (sAg.max(), lam_f[sAg.argmax()], sAu.max(), lam_f[sAu.argmax()]))
np.savetxt(OUT/'fig5a_mie_check.csv', np.c_[lam_m, dfm45, mie45, rel], delimiter=',', comments='', fmt='%.6e',
           header='lambda0_nm,sigma_abs_DFM_aAu45nm_um2,sigma_abs_coated_Mie_um2,rel_err')

# ------------------------------------------------------------------ plotting
plt.rcParams.update({'font.family': 'serif',
                     'font.serif': ['TeX Gyre Termes', 'Liberation Serif', 'DejaVu Serif'],
                     'mathtext.fontset': 'stix', 'font.size': 9.5, 'axes.linewidth': 0.9,
                     'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
                     'xtick.major.size': 4, 'ytick.major.size': 4, 'xtick.minor.size': 2, 'ytick.minor.size': 2,
                     'legend.fontsize': 8.5})
ink = '#0b0b0b'; mid = '#52514e'; agc = '#c9c9c9'; auc = '#f28e1c'
cmap = plt.get_cmap('viridis')
cols = [cmap(0.82*i/(len(radii) - 1)) for i in range(len(radii))]

fig, axs = plt.subplots(1, 2, figsize=(6.5, 2.9), sharey=True)
fig.subplots_adjust(left=0.095, right=0.985, bottom=0.16, top=0.975, wspace=0.13)

def frame(ax):
    ax.set_xlim(300, 800); ax.set_ylim(0, 0.025)
    ax.xaxis.set_major_locator(MultipleLocator(100)); ax.xaxis.set_minor_locator(MultipleLocator(20))
    ax.yaxis.set_major_locator(MultipleLocator(0.005)); ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.set_xlabel(r'Free-space wavelength, $\lambda_0$ (nm)')

def arrow(ax, p, q, **kw):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle=kw.pop('style', '-|>'), mutation_scale=kw.pop('ms', 8),
                                 shrinkA=0, shrinkB=0, color=kw.pop('color', ink), lw=kw.pop('lw', 0.9), **kw))

# (a)
ax = axs[0]; frame(ax)
for r, c in zip(radii, cols):
    d = zones[float(r)]
    ax.plot(d[:, 0], d[:, 4], color=c, lw=1.2)
ax.plot(lam_m, mie45, 'o', ms=3.6, mfc='white', mec=ink, mew=0.8, zorder=5,
        label='coated-sphere Mie\n' + r'($a_{\rm Au}=45$ nm)')
arrow(ax, (485, 0.0068), (528, 0.0205), color=ink, lw=0.8, style='->', ms=9)
ax.text(478, 0.0232, r'$a_{\rm Au} = 10, 15, \ldots, 45$ nm', ha='center', va='center', fontsize=8.5)
ax.text(478, 0.0214, r'($a_{\rm Ag} = 50$ nm, $h = 5$ nm)', ha='center', va='center', fontsize=8.5, color=mid)
ax.legend(loc='upper right', frameon=False, fontsize=8, handletextpad=0.3, borderaxespad=0.25, handlelength=1.2)
ax.set_ylabel(r'$\sigma_{\rm abs}$ (' + 'µ' + r'm$^2$)')
ax.text(0.025, 0.97, '(a)', transform=ax.transAxes, ha='left', va='top', fontsize=11, fontweight='bold')
# inset: eccentric core-shell sphere
ia = ax.inset_axes([0.62, 0.50, 0.30, 0.36]); ia.set_xlim(-1.25, 1.75); ia.set_ylim(-1.5, 1.2)
ia.set_aspect('equal'); ia.axis('off')
ia.add_patch(Circle((0, 0), 1.0, fc=agc, ec='none'))
rc = 0.55; xc = 1.0 - 0.15 - rc                 # thinnest shell part h (exaggerated in the sketch)
ia.add_patch(Circle((xc, 0), rc, fc=auc, ec='none'))
arrow(ia, (0, 0), (-0.62, 0.62), ms=5, lw=0.7, style='->')
ia.text(-0.58, 0.83, r'$a_{\rm Ag}$', fontsize=8, ha='center', va='center')
arrow(ia, (xc, 0), (xc, -rc), ms=5, lw=0.7, style='->')
ia.text(xc + 0.05, -0.30, r'$a_{\rm Au}$', fontsize=8, ha='left', va='center')
ia.text(xc - 0.18, 0.18, 'Au', fontsize=8, ha='center', va='center')
ia.text(-0.55, -0.25, 'Ag', fontsize=8, ha='center', va='center')
ia.plot([0.925, 1.22], [0.0, 0.42], color=ink, lw=0.5)
ia.text(1.25, 0.50, r'$h$', fontsize=8, ha='left', va='center')
arrow(ia, (-0.95, -1.35), (-0.95, -0.95), ms=5, lw=0.8)
arrow(ia, (-0.95, -1.35), (-0.55, -1.35), ms=5, lw=0.8)
ia.text(-1.05, -0.88, r'$\mathbfit{k}^{\rm inc}$', fontsize=8, ha='right', va='center')
ia.text(-0.48, -1.35, r'$\mathbfit{E}^{\rm inc}$', fontsize=8, ha='left', va='center')

# (b)
ax = axs[1]; frame(ax)
ax.plot(lam_f, sAg, color='#6d6d6d', lw=1.3, ls=(0, (7, 3)), label=r'Ag, $a = 50$ nm')
ax.plot(lam_f, sAu, color=auc, lw=1.3, ls='-', label=r'Au, $a = 50$ nm')
ax.legend(loc='upper right', frameon=False, title='Mie theory', title_fontsize=8.5)
ax.text(0.025, 0.97, '(b)', transform=ax.transAxes, ha='left', va='top', fontsize=11, fontweight='bold')
ib = ax.inset_axes([0.74, 0.30, 0.24, 0.32]); ib.set_xlim(-1.25, 1.5); ib.set_ylim(-1.5, 1.2)
ib.set_aspect('equal'); ib.axis('off')
ib.add_patch(Circle((0, 0), 1.0, fc=agc, ec='none'))
arrow(ib, (0, 0), (-0.62, 0.62), ms=5, lw=0.7, style='->')
ib.text(-0.25, 0.55, r'$a$', fontsize=8, ha='center', va='center')
arrow(ib, (-0.95, -1.35), (-0.95, -0.95), ms=5, lw=0.8)
arrow(ib, (-0.95, -1.35), (-0.55, -1.35), ms=5, lw=0.8)
ib.text(-1.05, -0.88, r'$\mathbfit{k}^{\rm inc}$', fontsize=8, ha='right', va='center')
ib.text(-0.48, -1.35, r'$\mathbfit{E}^{\rm inc}$', fontsize=8, ha='left', va='center')

fig.savefig(OUT/'Fig_spr_coreshell.pdf'); fig.savefig(OUT/'Fig_spr_coreshell.png', dpi=600)
print('saved output/Fig_spr_coreshell.pdf/.png')
