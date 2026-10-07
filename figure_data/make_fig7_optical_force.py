#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Qiang Sun
# SPDX-License-Identifier: BSD-3-Clause
"""Manuscript Fig. 7: optical force and torque on two dielectric rods (spherocylinders,
refractive index 1.5) in water, free-space wavelength 840 nm.

Data: fig7_Rslt_FrcTrq_archived.dat, the archived output Rslt_FrcTrq.dat of the case
inputs/fig06_optical_force/pill (73 orientations, 0-360 deg in 5-deg steps; the figure
shows 0-180 deg).  Columns: theta (deg), then F_x, F_y, F_z, N_x, N_y, N_z of rod 1
and of rod 2, in units of eps0|E0|^2 a^2 (force) and eps0|E0|^2 a^3 (torque),
a = 0.5 um.

(a) Vector schematic: spherocylinders of length 2a and width 2b, b/a = 0.5, centres at
    x = -/+1.01a; rod 1 turned by theta about +y, i.e. clockwise in this x-z view
    (code/optical_force/Main_EM.f90: anglecal_y(1) = theta, rotation in Geom_Mesh.f90).
(b), (c) The archived forces F_x and F_z, replotted unchanged.
(d) The torque N_y about each rod's own centre.  The archived torques are taken about
    the global origin: Input_Geom.dat gives the 'centre of mass' of BOTH rods as
    (0, 0, 0), while the rods are translated to x = -/+0.505 um = -/+1.01a.  Hence
        N_own = N_O - r_c x F,   r_c = (x_c, 0, 0),   N_y,own = N_y,O + x_c F_z,
    with x_c = -1.01a for rod 1 and +1.01a for rod 2.

Usage, from this directory:

    python3 make_fig7_optical_force.py

writes output/Fig_optical_force.pdf, output/Fig_optical_force.png and the plotted
values, output/fig7_plotted_data.csv (compare with the archived copy
fig7_plotted_data.csv), which also lists the torques about the origin.
Requires NumPy and Matplotlib.  See README.md in this directory for the provenance
of the data file.
"""
from pathlib import Path

import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Arc, FancyArrowPatch
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

HERE = Path(__file__).resolve().parent
OUT = HERE / 'output'
OUT.mkdir(exist_ok=True)

# ------------------------------------------------------------------ data
d = np.loadtxt(HERE/'fig7_Rslt_FrcTrq_archived.dat', skiprows=1, delimiter=',', usecols=range(13))
th, Fx1, Fy1, Fz1, Nx1, Ny1, Nz1, Fx2, Fy2, Fz2, Nx2, Ny2, Nz2 = d.T
xc1, xc2 = -1.01, +1.01
N1 = Ny1 + xc1*Fz1            # about rod 1's own centre
N2 = Ny2 + xc2*Fz2            # about rod 2's own centre
sel = th <= 180 + 1e-9
t = th[sel]
np.savetxt(OUT/'fig7_plotted_data.csv',
           np.c_[t, Fx1[sel], Fx2[sel], Fz1[sel], Fz2[sel], N1[sel], N2[sel], Ny1[sel], Ny2[sel]],
           delimiter=',', comments='', fmt='%.6e',
           header='theta_deg,Fx_1,Fx_2,Fz_1,Fz_2,Ny_1_own_centre,Ny_2_own_centre,'
                  'Ny_1_about_origin,Ny_2_about_origin')

def zeros(x, y):
    out = []
    for i in range(len(x)-1):
        if np.sign(y[i]) != np.sign(y[i+1]):
            out.append((x[i] - y[i]*(x[i+1]-x[i])/(y[i+1]-y[i]), 'down' if y[i+1] < y[i] else 'up'))
    return out
print('rod 1 own-centre N_y: %.4f..%.4f; zero crossings:' % (N1[sel].min(), N1[sel].max()),
      ['%.1f deg (%s)' % z for z in zeros(t, N1[sel])])
print('rod 2 own-centre N_y: %.4f..%.4f' % (N2[sel].min(), N2[sel].max()))
print('F_x1: %.4f (%.0f deg) .. %.4f (%.0f deg)' % (Fx1[sel].max(), t[Fx1[sel].argmax()], Fx1[sel].min(), t[Fx1[sel].argmin()]))
print('F_z1: min %.4f (%.0f deg); F_z2: %.4f (%.0f) .. %.4f (%.0f)' % (Fz1[sel].min(), t[Fz1[sel].argmin()],
      Fz2[sel].min(), t[Fz2[sel].argmin()], Fz2[sel].max(), t[Fz2[sel].argmax()]))
# head-tail symmetry of the rod: theta and theta+180 are the same geometry (different mesh orientation)
print('max |q(theta)-q(theta+180)|: F_x1 %.1e, F_z1 %.1e, N_1 %.1e' % (
      np.abs(Fx1[:37]-Fx1[36:73]).max(), np.abs(Fz1[:37]-Fz1[36:73]).max(), np.abs(N1[:37]-N1[36:73]).max()))

# ------------------------------------------------------------------ style
plt.rcParams.update({'font.family': 'serif',
                     'font.serif': ['TeX Gyre Termes', 'Liberation Serif', 'DejaVu Serif'],
                     'mathtext.fontset': 'stix', 'font.size': 9.5, 'axes.linewidth': 0.9,
                     'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
                     'xtick.major.size': 4, 'ytick.major.size': 4, 'xtick.minor.size': 2, 'ytick.minor.size': 2,
                     'legend.fontsize': 8.5, 'legend.handlelength': 2.6})
ink = '#0b0b0b'; mid = '#52514e'; fillc = '#d9d9d9'; edgec = '#3f3f3f'; zc = '#9a9a9a'
ls1 = (0, (6, 3))            # particle 1: dashed
ls2 = '-'                    # particle 2: solid

fig = plt.figure(figsize=(6.5, 5.7))
gs = fig.add_gridspec(2, 2, left=0.095, right=0.985, bottom=0.08, top=0.985, wspace=0.32, hspace=0.32)
axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0]); axD = fig.add_subplot(gs[1, 1])

# ------------------------------------------------------------------ (a) schematic
def stadium(cx, cy, theta_deg, a=1.0, b=0.5, n=80):
    """Outline of a spherocylinder of tip-to-tip length 2a and width 2b in the x-z plane (x right,
    z up), turned by theta about +y (right-handed about +y = clockwise in this view)."""
    L = a - b
    t1 = np.linspace(-np.pi/2, np.pi/2, n); t2 = np.linspace(np.pi/2, 3*np.pi/2, n)
    pts = np.vstack([np.c_[L + b*np.cos(t1), b*np.sin(t1)], np.c_[-L + b*np.cos(t2), b*np.sin(t2)]])
    q = -np.radians(theta_deg)
    R = np.array([[np.cos(q), -np.sin(q)], [np.sin(q), np.cos(q)]])
    return pts @ R.T + [cx, cy]

def rod(ax, cx, cy, theta, s=1.0, lw=1.1, label=None, fs=8, dx=0.0, dy=0.0):
    ax.add_patch(Polygon(stadium(cx, cy, theta, a=s, b=0.5*s), closed=True, fc=fillc, ec=edgec, lw=lw, zorder=2))
    if label:
        ax.text(cx+dx, cy+dy, label, ha='center', va='center', fontsize=fs, color='#262626', zorder=4)

ax = axA
ax.set_position([0.005, 0.545, 0.475, 0.445])
ax.set_xlim(-3.6, 3.6); ax.set_ylim(-3.05, 2.85); ax.set_aspect('equal'); ax.axis('off')
s = 0.36
for (cx, ang) in [(-2.72, 0), (-0.9, 45), (0.92, 90), (2.72, 135)]:
    cy = 2.2
    rod(ax, cx - 1.01*s, cy, ang, s=s, lw=0.8, label='1', fs=8)
    rod(ax, cx + 1.01*s, cy, 0, s=s, lw=0.8, label='2', fs=8)
    ax.text(cx, cy - 0.60, r'$\theta = %d^\circ$' % ang, ha='center', va='center', fontsize=9)

# large schematic (rod 1 drawn at theta = 45 deg)
yc = -0.38; c1 = np.array([-1.65, yc]); c2 = np.array([0.45, yc]); T = 45
rod(ax, *c1, T); rod(ax, *c2, 0)
ax.plot([c1[0], c2[0]], [yc, yc], color=ink, lw=0.8, ls=(0, (7, 2, 1.5, 2)), zorder=3)          # line of centres
u = np.array([np.cos(np.radians(T)), -np.sin(np.radians(T))])
ax.plot([c1[0], c1[0] + u[0]], [yc, yc + u[1]], color=ink, lw=0.8, ls=(0, (3, 2)), zorder=3)      # rod-1 axis
ax.add_patch(Arc(c1, 0.9, 0.9, theta1=-T, theta2=0, color=ink, lw=0.8, zorder=3))
ax.text(c1[0] + 0.72*np.cos(np.radians(T/2)), yc - 0.72*np.sin(np.radians(T/2)), r'$\theta$',
        ha='center', va='center', fontsize=10.5, zorder=4)
ax.text(c1[0] - 0.30, yc + 0.30, '1', ha='center', va='center', fontsize=9.5, color='#262626', zorder=4)
ax.text(c2[0], yc + 0.25, '2', ha='center', va='center', fontsize=9.5, color='#262626', zorder=4)
# 2b: width of rod 2, dimensioned outside its right cap
xb = c2[0] + 1.30
for yy in (yc - 0.5, yc + 0.5):
    ax.plot([c2[0] + 0.5, xb + 0.12], [yy, yy], color=ink, lw=0.55, ls=(0, (4, 1.5, 1, 1.5)), zorder=1)
ax.add_patch(FancyArrowPatch((xb, yc - 0.5), (xb, yc + 0.5), arrowstyle='<->', mutation_scale=7,
                             lw=0.8, color=ink, shrinkA=0, shrinkB=0, zorder=4))
ax.text(xb + 0.12, yc, r'$2b$', ha='left', va='center', fontsize=9.5, zorder=4)
# 2a: tip-to-tip length of rod 2
x0, x1 = c2[0] - 1.0, c2[0] + 1.0; ya = yc - 0.85
for xx in (x0, x1):
    ax.plot([xx, xx], [yc - 0.02, ya - 0.12], color=ink, lw=0.55, ls=(0, (4, 1.5, 1, 1.5)), zorder=1)
ax.add_patch(FancyArrowPatch((x0, ya), (x1, ya), arrowstyle='<->', mutation_scale=7, lw=0.8, color=ink,
                             shrinkA=0, shrinkB=0, zorder=4))
ax.text(c2[0], ya - 0.27, r'$2a$', ha='center', va='center', fontsize=9.5)
ax.text(x1 + 0.35, ya, r'$b/a = 0.5$', ha='left', va='center', fontsize=9.5)
# incident plane wave
o = np.array([-3.35, -2.78])
ax.add_patch(FancyArrowPatch(o, o + [0, 0.72], arrowstyle='-|>', mutation_scale=10, lw=1.1, color=ink,
                             shrinkA=0, shrinkB=0))
ax.add_patch(FancyArrowPatch(o, o + [0.72, 0], arrowstyle='-|>', mutation_scale=10, lw=1.1, color=ink,
                             shrinkA=0, shrinkB=0))
ax.text(o[0] + 0.15, o[1] + 0.68, r'$\mathbfit{k}^{\mathrm{inc}} = (0, 0, k_{\mathrm{out}})$', ha='left', va='center', fontsize=9)
ax.text(o[0] + 0.85, o[1], r'$\mathbfit{E}^{\mathrm{inc}} = (E_0, 0, 0)$', ha='left', va='center', fontsize=9)
ax.text(-3.6, 2.85, '(a)', ha='left', va='top', fontsize=11, fontweight='bold')

# ------------------------------------------------------------------ (b)-(d)
def panel(ax, y1, y2, ylab, lab, ylim, yticks, legloc, zero=False, note=None):
    if zero:
        ax.axhline(0, color=zc, lw=0.6, ls=(0, (1, 1.5)), zorder=0)
    ax.plot(t, y1, color=ink, lw=1.3, ls=ls1, label='Rod 1 (rotated)')
    ax.plot(t, y2, color=ink, lw=1.3, ls=ls2, label='Rod 2 (fixed)')
    ax.set_xlim(0, 180); ax.set_xticks([0, 45, 90, 135, 180]); ax.xaxis.set_minor_locator(MultipleLocator(15))
    ax.set_ylim(*ylim); ax.set_yticks(yticks); ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.set_xlabel(r'Relative orientation, $\theta$ (deg)')
    ax.set_ylabel(ylab)
    ax.legend(loc=legloc, frameon=False)
    ax.text(0.03, 0.97, lab, transform=ax.transAxes, ha='left', va='top', fontsize=11, fontweight='bold')
    if note:
        ax.text(0.97, 0.97, note, transform=ax.transAxes, ha='right', va='top', fontsize=8.5, color=mid)

panel(axB, Fx1[sel], Fx2[sel], r'$F_x/[\varepsilon_0|E_0|^2a^2]$', '(b)', (-0.07, 0.07),
      [-0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06], 'lower left', zero=True)
panel(axC, Fz1[sel], Fz2[sel], r'$F_z/[\varepsilon_0|E_0|^2a^2]$', '(c)', (0.066, 0.082),
      [0.066, 0.070, 0.074, 0.078, 0.082], 'lower center')
panel(axD, N1[sel], N2[sel], r'$N_y/[\varepsilon_0|E_0|^2a^3]$', '(d)', (-0.03, 0.03),
      [-0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03], 'lower left', zero=True,
      note="about each rod's own center")
def fmt(nd):
    return matplotlib.ticker.FuncFormatter(lambda v, p: ('%.*f' % (nd, 0.0 if abs(v) < 1e-12 else v)).replace('-', '−'))
axB.yaxis.set_major_formatter(fmt(2)); axC.yaxis.set_major_formatter(fmt(3)); axD.yaxis.set_major_formatter(fmt(2))

fig.savefig(OUT/'Fig_optical_force.pdf')
fig.savefig(OUT/'Fig_optical_force.png', dpi=600)
print('saved output/Fig_optical_force.pdf/.png')
