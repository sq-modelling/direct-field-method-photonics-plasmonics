#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Qiang Sun
# SPDX-License-Identifier: BSD-3-Clause
"""Manuscript Fig. 3: scattering by three nearly touching PEC rods, k_out a = 1e-3, 5 and 50.

(a)-(c) Re(E_x^tot/E0) in the plane y = 0 (120 x 200 samples, |x| <= 3a, -2.5a <= z <= 7.5a) on two colour scales
        ((a): 0 to 2, 1 = incident amplitude; (b), (c): -2 to 2, values beyond the end arrows).  Data: fig3_slice_y0_ka{0001,5,50}.csv.gz (columns x/a, z/a, pos,
        Re(E_x^tot/E0); pos > 0 inside the rods, where the value is 0 and the cell is masked).
(d)-(f) Bistatic radar cross section in the x-z plane, 1-degree steps (theta: x = r cos theta,
        z = r sin theta, so 90 deg is forward), as 10 log10(sigma/a^2) over 60 dB per panel.
        Data: fig3_rcs_zx.csv, columns sigma_over_a2_ka0001_far (far-field limit),
        sigma_over_a2_ka5_r1000 and sigma_over_a2_ka50_r1000 (r = 1000a); the two further columns are
        quoted in the manuscript text but not plotted.

k_out a = 5 and 50 are outputs of code/dfm_field_solver for inputs/fig02_pec_rods/{ka5,ka50} (E solve);
k_out a = 1e-3 is computed from the surface current of the magnetic-field solution for
inputs/fig02_pec_rods/ka0001 with ../postprocessing/pec_current_route.py.  See README.md in this
directory for the provenance.  Usage, from this directory:

    python3 make_fig3_three_rods.py

writes output/Fig_three_rods.pdf, output/Fig_three_rods.png and output/fig3_summary.csv (compare with the
copy fig3_summary.csv).  Requires NumPy and Matplotlib.
"""
from pathlib import Path
import glob, os, subprocess
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import matplotlib.font_manager as fm
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import FancyBboxPatch


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
plt.rcParams.update({'font.family': 'serif',
                     'font.serif': ['TeX Gyre Termes', 'Liberation Serif', 'DejaVu Serif'],
                     'mathtext.fontset': 'stix', 'font.size': 9.5, 'axes.linewidth': 0.8,
                     'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True})
ink = '#20262b'; rodfill = '#d9dee3'; red = '#a51c30'; grey = '#e1e5e8'
RODS = ((-2.505, 3.0), (0.0, 2.0), (2.505, 3.0)); R = 0.3      # centre, length (units of a); radius
VMAX = 2.0                                                       # colour scale of (b), (c): -VMAX .. +VMAX
A_NORM = TwoSlopeNorm(vmin=0.0, vcenter=1.0, vmax=2.0)           # colour scale of (a): 1 = incident amplitude
CASES = (  # file tag, field title, RCS title, dB range of the polar plot (centre, outer circle; 60 dB in each)
    ('ka0001', r'(a) $k_{\rm out}a = 10^{-3}$', r'(d) $k_{\rm out}a = 10^{-3}$, far field', (-140.0, -80.0)),
    ('ka5', r'(b) $k_{\rm out}a = 5$', r'(e) $k_{\rm out}a = 5$, $r = 1000a$', (-30.0, 30.0)),
    ('ka50', r'(c) $k_{\rm out}a = 50$', r'(f) $k_{\rm out}a = 50$, $r = 1000a$', (-15.0, 45.0)))
DB_STEP = 20.0                                                   # rings every 20 dB from the centre value


HERE = Path(__file__).resolve().parent
OUT = HERE / 'output'
OUT.mkdir(exist_ok=True)
RCS = np.genfromtxt(HERE / 'fig3_rcs_zx.csv', delimiter=',', names=True)


def load(case):
    a = np.loadtxt(HERE / ('fig3_slice_y0_%s.csv.gz' % case), delimiter=',', skiprows=1)
    x = np.unique(a[:, 0]); z = np.unique(a[:, 1])
    assert len(x) * len(z) == len(a) and np.array_equal(a[:len(x), 0], x)
    ex = np.where(a[:, 2] > 0, np.nan, a[:, 3]).reshape(len(z), len(x))    # interior of the rods masked
    col = {'ka0001': 'sigma_over_a2_ka0001_far', 'ka5': 'sigma_over_a2_ka5_r1000',
           'ka50': 'sigma_over_a2_ka50_r1000'}[case]
    return x, z, ex, np.deg2rad(RCS['theta_deg']), RCS[col]


def capsule(ax, c, L, z0=0.0, r=R, lw=0.8, zorder=6):
    ax.add_patch(FancyBboxPatch((c - L/2, z0 - r), L, 2*r, boxstyle='round,pad=0,rounding_size=%g' % r,
                                fc=rodfill, ec=ink, lw=lw, zorder=zorder))


cmap = plt.get_cmap('RdBu_r').copy(); cmap.set_bad(rodfill)


def field_panel(ax, x, z, ex, title, ylab, norm):
    ax.set_facecolor(grey)                  # outside the stored slice |x| <= 3a
    X, Z = np.meshgrid(x, z)
    im = ax.pcolormesh(X, Z, np.ma.masked_invalid(ex), cmap=cmap,
                       norm=norm, shading='nearest', rasterized=True)
    for c, L in RODS:
        capsule(ax, c, L)
    ax.set_xlim(-4.2, 4.2); ax.set_ylim(-2.5, 7.5); ax.set_aspect('equal')
    ax.set_xticks([-4, -2, 0, 2, 4]); ax.set_yticks([-2, 0, 2, 4, 6])
    ax.set_xlabel(r'$x/a$', labelpad=2)
    if ylab:
        ax.set_ylabel(r'$z/a$', labelpad=2)
    else:
        ax.set_yticklabels([])
    ax.set_title(title, loc='left', pad=3, fontsize=10.5)
    # incident-wave indicator (lower left)
    o = (-2.55, -1.72); st = dict(arrowstyle='-|>', color=ink, lw=1.1, mutation_scale=10, shrinkA=0, shrinkB=0)
    halo = [pe.Stroke(linewidth=3, foreground='white'), pe.Normal()]
    for tip in [(o[0], -0.62), (-1.35, o[1])]:
        a = ax.annotate('', xy=tip, xytext=o, arrowprops=st, zorder=9); a.arrow_patch.set_path_effects(halo)
    for (xx, yy, s, ha, va) in [(-2.42, -1.05, r'$\mathbfit{k}^{\rm inc}$', 'left', 'center'),
                                (-1.95, -1.95, r'$\mathbfit{E}^{\rm inc}$', 'center', 'top')]:
        t = ax.text(xx, yy, s, ha=ha, va=va, fontsize=9, color=ink, zorder=10); t.set_path_effects(halo)
    return im


def polar_panel(ax, th, s, title, rng):
    """Polar plot of 10 log10(sigma/a^2), clipped to rng = (centre value, outer-circle value)."""
    floor, ceil = rng
    db = 10*np.log10(np.maximum(s, 1e-300))
    ax.plot(th, np.clip(db, floor, ceil) - floor, color=red, lw=0.9)
    ax.set_theta_zero_location('E'); ax.set_theta_direction(1)
    ax.set_thetagrids([0, 90, 180, 270], [r'$+x$', r'$+z$', r'$-x$', r'$-z$'], fontsize=9)
    ax.tick_params(axis='x', pad=-1)
    ax.set_rlim(0, ceil - floor)
    ax.set_rticks(np.arange(DB_STEP, ceil - floor - 1e-9, DB_STEP)); ax.set_yticklabels([])
    ax.grid(color='#aab2b8', lw=0.4)
    ax.set_title(title, va='bottom', pad=9, fontsize=10)


def db_scale(fig, cx, y, rad, rng):
    """Radial dB scale below a polar plot: from its centre (x = cx) to its outer circle (x = cx + rad), inches."""
    floor, ceil = rng
    sx = fig.add_axes(box(cx, y, rad, 0.001)); sx.set_xlim(0, ceil - floor)
    for side in ('left', 'right', 'top'):
        sx.spines[side].set_visible(False)
    sx.spines['bottom'].set_linewidth(0.7); sx.set_yticks([])
    vals = np.arange(floor, ceil + 1e-9, DB_STEP)
    sx.set_xticks(vals - floor); sx.set_xticklabels(['$%d$' % v for v in vals])
    sx.tick_params(axis='x', which='both', direction='in', length=3.5, width=0.7, labelsize=8, pad=2, top=False)
    fig.text((cx - 0.17) / W, (y + 0.01) / H, r'$10\,\log_{10}(\sigma/a^2)$ (dB)', ha='right', va='center', fontsize=8.5,
             color=ink)


data = [load(c[0]) for c in CASES]
KOUTA = {'ka0001': '0.001', 'ka5': '5', 'ka50': '50'}
rows = ['k_out_a,rcs_max_dB,rcs_max_theta_deg,rcs_forward_dB,rcs_backward_dB,Re_Ex_min,Re_Ex_max']
for (tag, *_), (x, z, ex, th, s) in zip(CASES, data):
    assert np.array_equal(x, data[1][0]) and np.array_equal(z, data[1][1])
    db = 10*np.log10(s)
    print('%-8s RCS max %7.2f dB at %5.1f deg; forward %7.2f dB, backscatter %7.2f dB; Re Ex in [%.2f, %.2f]'
          % (tag, db.max(), np.rad2deg(th[db.argmax()]), db[90], db[270], np.nanmin(ex), np.nanmax(ex)))
    rows.append('%s,%.4f,%.1f,%.4f,%.4f,%.4f,%.4f' % (KOUTA[tag], db.max(), np.rad2deg(th[db.argmax()]), db[90], db[270],
                                                     np.nanmin(ex), np.nanmax(ex)))
(OUT / 'fig3_summary.csv').write_text('\n'.join(rows) + '\n')

W, H = 7.1, 7.12
fig = plt.figure(figsize=(W, H))
def box(l, b, w, h):                        # axes rectangle given in inches
    return [l / W, b / H, w / W, h / H]

# dimensioned schematic (gaps enlarged for clarity); 9.8 x 1.67 data units at 0.48 in per unit
g = fig.add_axes(box(1.13, 6.26, 4.70, 0.80)); g.axis('off')
disp = ((-2.62, 3.0), (0.0, 2.0), (2.62, 3.0))
for i, (c, L) in zip((1, 3, 2), disp):     # particle numbers of Input_Geom.dat (and of the slice column pos)
    capsule(g, c, L, lw=0.9, zorder=2)
    g.text(c, 0, str(i), ha='center', va='center', fontsize=10, fontweight='bold', color=ink, zorder=3)
def dim(ax, p, q, lab, off):
    ax.annotate('', xy=q, xytext=p, arrowprops=dict(arrowstyle='|-|', color=ink, lw=0.8, mutation_scale=5,
                shrinkA=0, shrinkB=0))
    ax.text((p[0]+q[0])/2 + off[0], (p[1]+q[1])/2 + off[1], lab, ha='center', va='center', fontsize=9.5, color=ink)
dim(g, (-4.12, 0.45), (-1.12, 0.45), r'$3a$', (0, 0.14)); dim(g, (-1.0, 0.45), (1.0, 0.45), r'$2a$', (0, 0.14))
dim(g, (1.12, 0.45), (4.12, 0.45), r'$3a$', (0, 0.14)); dim(g, (4.32, -0.3), (4.32, 0.3), r'$0.6a$', (0.44, 0))
for gx in (-1.06, 1.06):
    g.annotate('', xy=(gx, -0.31), xytext=(gx, -0.6), arrowprops=dict(arrowstyle='-|>', color=red, lw=0.8,
               mutation_scale=7))
g.text(0, -0.66, r'two gaps, $g = 0.005a$ (enlarged for clarity)', ha='center', va='top', fontsize=9, color='#5c2027')
g.set_xlim(-4.9, 4.9); g.set_ylim(-0.95, 0.72); g.set_aspect('equal')

# field maps: 8.4a x 10a at 2.0 in width -> 2.381 in height (equal aspect)
fw = 2.0; fh = fw * 10.0 / 8.4; fb = 3.15; gap = 0.27; left = 0.50
lefts = [left + j * (fw + gap) for j in range(3)]
ims = []
for j, ((tag, ftitle, *_), (x, z, ex, th, s)) in enumerate(zip(CASES, data)):
    ax = fig.add_axes(box(lefts[j], fb, fw, fh))
    norm = A_NORM if j == 0 else TwoSlopeNorm(vmin=-VMAX, vcenter=0, vmax=VMAX)
    ims.append(field_panel(ax, x, z, ex, ftitle, j == 0, norm))
cb_y = 5.90
def colourbar(im, l, w, ticks, extend):
    cax = fig.add_axes(box(l, cb_y, w, 0.085))
    cb = fig.colorbar(im, cax=cax, orientation='horizontal', ticks=ticks, extend=extend, extendfrac=0.06)
    cb.ax.xaxis.set_ticks_position('top'); cb.ax.tick_params(labelsize=9, width=0.7, direction='out', pad=1.5)
    cb.outline.set_linewidth(0.6)
    fig.text((l - 0.12) / W, (cb_y + 0.0425) / H, r'$\mathrm{Re}(E_x^{\rm tot}/E_0)$', ha='right', va='center',
             fontsize=10)
colourbar(ims[0], 1.40, 1.05, [0, 1, 2], 'max')          # (a)
colourbar(ims[1], 4.05, 2.55, [-2, -1, 0, 1, 2], 'both')    # (b), (c)

# far-field polar plots, centred under the field maps, each with its radial dB scale below
pd = 1.78; pb = 0.50
for j, ((tag, ftitle, ptitle, rng), (x, z, ex, th, s)) in enumerate(zip(CASES, data)):
    cx = lefts[j] + fw/2
    axp = fig.add_axes(box(cx - pd/2, pb, pd, pd), projection='polar')
    polar_panel(axp, th, s, ptitle, rng)
    db_scale(fig, cx, pb - 0.24, pd/2, rng)
fig.savefig(OUT / 'Fig_three_rods.pdf', dpi=600)   # dpi of the rasterized field maps (10 px per cell)
fig.savefig(OUT / 'Fig_three_rods.png', dpi=600)
print('wrote output/Fig_three_rods.pdf, output/Fig_three_rods.png, output/fig3_summary.csv; font:',
      fm.findfont(fm.FontProperties(family='serif')))
