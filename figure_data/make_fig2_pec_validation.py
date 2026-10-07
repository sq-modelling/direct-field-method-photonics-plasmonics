#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Qiang Sun
# SPDX-License-Identifier: BSD-3-Clause
"""Manuscript Fig. 2: validation of the Direct Field Method for a PEC sphere.

(a) Total surface field |E_tot|/E0 on a PEC sphere of radius a, k_out a = 15,
    2562-node Q6 mesh (the case of Fig. 1(a)), along the great circle in the
    E-plane (y = 0): DFM nodes with |y| < 0.03a against the exact Mie solution.
    Data: fig2a_pec_sphere_ka15_surface_field.csv (x, y, z and |E_tot| of all
    2562 surface nodes, from Rslt_SurfCal3_EM.dat of the archived run).
(b) Mesh convergence at k_out a = 4.5: relative error of |E_sc| at the
    near-field point (0, 0, 1.2a) against the number of surface nodes N
    (mesh levels 4-14), with the run times of three meshes.
    Data: fig2b_pec_sphere_ka4p5_convergence.csv (archived DFM values).

The Mie reference values are computed here.  Usage, from this directory:

    python3 make_fig2_pec_validation.py

writes output/Fig_validation_PEC.pdf, output/Fig_validation_PEC.png and the
plotted values of panel (b), output/fig2b_plotted_data.csv (compare with the
archived copy fig2b_plotted_data.csv).  Requires NumPy, SciPy and Matplotlib.
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
from scipy.special import spherical_jn, spherical_yn

HERE = Path(__file__).resolve().parent
OUT = HERE / 'output'
OUT.mkdir(exist_ok=True)

def mie_pec_coeffs(x):
    nmax=int(x+4*x**(1/3)+15); n=np.arange(1,nmax+1)
    psi=lambda z: z*spherical_jn(n,z)
    dpsi=lambda z: spherical_jn(n,z)+z*spherical_jn(n,z,derivative=True)
    h=lambda z: spherical_jn(n,z)+1j*spherical_yn(n,z)
    dh=lambda z: spherical_jn(n,z,derivative=True)+1j*spherical_yn(n,z,derivative=True)
    xi=lambda z: z*h(z); dxi=lambda z: h(z)+z*dh(z)
    return n, dpsi(x)/dxi(x), psi(x)/xi(x), (1j**n)*(2*n+1)/(n*(n+1)), h, dxi, xi

def pi_n(n, mu):
    nmax=n[-1]; p=np.zeros((nmax,)+np.shape(mu)); p[0]=1
    if nmax>1: p[1]=3*mu
    for k in range(2,nmax):
        m=k+1; p[k]=(2*m-1)/(m-1)*mu*p[k-1]-m/(m-1)*p[k-2]
    return p

def Er_surface(x, theta, phi):
    n,an,bn,En,h,dxi,xi=mie_pec_coeffs(x)
    coef=En*n*(n+1)*(-1j*spherical_jn(n,x)/x+1j*an*h(x)/x)
    return np.cos(phi)*np.sin(theta)*np.einsum('n,n...->...',coef,pi_n(n,np.cos(theta)))

def Esc_axis(x, r):
    n,an,bn,En,h,dxi,xi=mie_pec_coeffs(x); rho=x*r
    return np.sum(En*(n*(n+1)/2)*(1j*an*dxi(rho)-bn*xi(rho))/rho)

# ---- (a) data: x/a, y/a, z/a, |E_tot|/E0 at the 2562 surface nodes
D=np.loadtxt(HERE/'fig2a_pec_sphere_ka15_surface_field.csv',delimiter=',',skiprows=1)
X,Y,Z,Emag=D[:,0],D[:,1],D[:,2],D[:,3]
r=np.sqrt(X**2+Y**2+Z**2); th=np.arccos(np.clip(Z/r,-1,1)); ph=np.arctan2(Y,X)
Er_all=Er_surface(15.0,th,ph)
err_all=np.abs(Emag-np.abs(Er_all))
relL2=np.linalg.norm(err_all)/np.linalg.norm(np.abs(Er_all))
sel=np.abs(Y)<0.03
sgn=np.where(X[sel]>=0,1,-1)
ang=np.degrees(th[sel])*sgn            # signed polar angle in the E-plane
tt=np.linspace(0,np.pi,1801)
Emie_pos=np.abs(Er_surface(15.0,tt,0.0)); Emie_neg=np.abs(Er_surface(15.0,tt,np.pi))
# ---- (b) data: archived DFM values of |E_sc| at (0, 0, 1.2a) and run times
C=np.loadtxt(HERE/'fig2b_pec_sphere_ka4p5_convergence.csv',delimiter=',',skiprows=1)
N=C[:,1].astype(int)
Enum=C[:,3]
cpu=C[:,4].astype(int)
Eref=abs(Esc_axis(4.5,1.2))
rel=np.abs(Enum-Eref)/Eref
p=np.polyfit(np.log(N),np.log(rel),1)[0]
print(f"(a) nodes in E-plane: {sel.sum()}, relative L2 error over all {len(Emag)} nodes: {100*relL2:.2f}%, max abs err {err_all.max():.4f}")
print(f"(b) Mie |E_sc|(0,0,1.2a)={Eref:.8f}; rel err {rel[0]:.2e} -> {rel[-1]:.2e}; slope {p:.2f}")
np.savetxt(OUT/'fig2b_plotted_data.csv',np.c_[N,Enum,rel,cpu],delimiter=',',header='N_nodes,abs_Esc_DFM,rel_err_vs_Mie,run_time_s',comments='',fmt=['%d','%.11f','%.4e','%d'])

# ---- plotting
plt.rcParams.update({'font.family':'serif','font.serif':['TeX Gyre Termes','Liberation Serif','DejaVu Serif'],
                     'mathtext.fontset':'stix','font.size':9,'axes.linewidth':0.8,
                     'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True})
ink='#0b0b0b'; c_dfm='#2a78d6'; c_ref='#0b0b0b'; grid='#dddddd'
fig,axs=plt.subplots(1,2,figsize=(6.5,2.7),gridspec_kw={'width_ratios':[1.35,1]})
ax=axs[0]
ax.plot(np.degrees(tt),Emie_pos,color=c_ref,lw=1.0,label='Mie (exact)')
ax.plot(-np.degrees(tt),Emie_neg,color=c_ref,lw=1.0)
ax.plot(ang,Emag[sel],'o',ms=3.2,mfc='none',mec=c_dfm,mew=0.9,label='DFM, 2562 nodes\n(relative $L_2$ error %.2f%%)'%(100*relL2))
ax.set_xlim(-180,180); ax.set_xticks([-180,-120,-60,0,60,120,180])
ax.set_xlabel(r'polar angle $\theta$ in the $E$-plane (deg)')
ax.set_ylabel(r'$|\mathbfit{E}^{\mathrm{tot}}|/E_0$ on $S$')
ax.set_ylim(0,2.45)
ax.legend(frameon=False,loc='upper center',fontsize=8)
ax.text(0.99,0.97,r'(a) $k_{\mathrm{out}}a=15$',transform=ax.transAxes,ha='right',va='top')
# relative L2 error shown in the legend entry, which avoids overlapping the data

ax=axs[1]
ax.loglog(N,rel,'o-',ms=3.5,color=c_dfm,lw=1.0,mfc='white',label='DFM vs Mie')
NN=np.array([600,9000]); ax.loglog(NN,rel[0]*(NN/N[0])**-2,'--',color='#888888',lw=0.9,label=r'$\propto N^{-2}\ (\propto h^{4})$')
for i,off in [(0,(7,-3)),(4,(7,3)),(10,(-30,-11))]:
    ax.annotate('<1 s' if i==0 else '%d s'%cpu[i],  # time of a complete run; the 642-node run took <1 s
                (N[i],rel[i]),textcoords='offset points',xytext=off,fontsize=8,color='#52514e')
ax.set_xlabel(r'number of surface nodes $N$')
ax.set_ylabel(r'relative error of $|\mathbfit{E}^{\mathrm{sc}}|$ at $(0,0,1.2a)$')
ax.set_xlim(500,10000)
ax.legend(frameon=False,loc='lower left',fontsize=8)
ax.text(0.97,0.97,r'(b) $k_{\mathrm{out}}a=4.5$',transform=ax.transAxes,ha='right',va='top')
fig.tight_layout(w_pad=1.5)
fig.savefig(OUT/'Fig_validation_PEC.pdf'); fig.savefig(OUT/'Fig_validation_PEC.png',dpi=600)
print('saved output/Fig_validation_PEC.pdf/.png')
