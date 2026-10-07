#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Qiang Sun
# SPDX-License-Identifier: BSD-3-Clause
"""Exact 'LgRd' rod shape of the release generator (Geom_Mesh.f90, SearchLgRd_cvsa/_rdbscl).

The surface is a Rankine ovoid (axisymmetric source-sink pair at z = +-c) with tips at z = +-1 and
half-width b at z = 0, in the unzoomed local frame (axis z):
   F(z, rho) = (1-c^2)^2 [ (z+c)/sqrt((z+c)^2+rho^2) - (z-c)/sqrt((z-c)^2+rho^2) ] - 2 c rho^2 = 0,
   (1-c^2)^2 = b^2 sqrt(b^2+c^2).
Here the meridian is computed along rays from the centre, (rho, z) = R(alpha) (sin alpha, cos alpha),
which is robust at the tips (alpha = 0, pi).  Provides the meridian polyline, its arc length s and the
inverse map s -> (rho, z) by monotone cubic interpolation (accuracy ~1e-12 relative)."""
import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import PchipInterpolator, CubicSpline


def cvsa(b):
    return brentq(lambda c: (1 - c**2)**2 - b**2 * np.sqrt(b**2 + c**2), 1e-6, 1.01, xtol=1e-15)


class Rankine:
    def __init__(self, b, nalpha=40001):
        self.b = b; self.c = cvsa(b)
        a = np.linspace(0.0, np.pi, nalpha)
        # cluster rays toward the tips where the surface turns fastest
        a = 0.5 * np.pi * (1 - np.cos(a))
        R = np.array([self.R_of_alpha(x) for x in a])
        rho = R * np.sin(a); z = R * np.cos(a)
        rho[0] = 0.0; rho[-1] = 0.0
        ds = np.hypot(np.diff(rho), np.diff(z))
        s = np.concatenate([[0.0], np.cumsum(ds)])
        # refine arc length with a spline of (rho,z) in alpha (chord -> arc correction)
        sp_r = CubicSpline(a, rho); sp_z = CubicSpline(a, z)
        am = 0.5 * (a[1:] + a[:-1]); h = np.diff(a)
        # Simpson on each sub-interval
        def speed(x):
            return np.hypot(sp_r(x, 1), sp_z(x, 1))
        ds2 = h / 6 * (speed(a[:-1]) + 4 * speed(am) + speed(a[1:]))
        s = np.concatenate([[0.0], np.cumsum(ds2)])
        self.alpha = a; self.Ra = R; self.rho_a = rho; self.z_a = z; self.s_a = s
        self.S = s[-1]
        self._a_of_s = PchipInterpolator(s, a)
        self._sp_r = sp_r; self._sp_z = sp_z

    def F(self, z, r):
        c = self.c
        return (1 - c**2)**2 * ((z + c) / np.hypot(z + c, r) - (z - c) / np.hypot(z - c, r)) - 2 * c * r**2

    def G(self, z, r):
        """F / rho^2, evaluated without cancellation near the axis (|z| > c uses
        1 - u/sqrt(u^2+rho^2) = rho^2 / (sqrt(u^2+rho^2) (sqrt(u^2+rho^2) + u)))."""
        c = self.c; z = abs(z)
        if z > c:
            u1, u2 = z + c, z - c
            r1, r2 = np.hypot(u1, r), np.hypot(u2, r)
            return (1 - c**2)**2 * (1 / (r2 * (r2 + u2)) - 1 / (r1 * (r1 + u1))) - 2 * c
        return self.F(z, r) / r**2 if r > 0 else np.inf

    def R_of_alpha(self, a):
        sa, ca = np.sin(a), np.cos(a)
        if sa < 1e-300:
            return 1.0
        g = lambda R: self.G(R * ca, R * sa)
        return brentq(g, 0.5 * self.b, 1.2, xtol=1e-15, rtol=1e-15)

    def point_at_s(self, s):
        """(rho, z) on the meridian at arc length s from the z=+1 tip (s in [0, S])."""
        a = self._a_of_s(np.clip(s, 0, self.S))
        R = np.array([self.R_of_alpha(x) for x in np.atleast_1d(a)])
        return R * np.sin(np.atleast_1d(a)), R * np.cos(np.atleast_1d(a))

    def project(self, rho, z):
        """Exact surface point on the ray through (rho, z) from the origin."""
        a = np.arctan2(rho, z)
        R = self.R_of_alpha(a)
        return R * np.sin(a), R * np.cos(a)

    def nose_radius(self):
        z = 1 - np.array([1e-6, 1e-5])
        a = np.array([np.arctan2(1e-3, 1.0)])
        rr, zz = self.point_at_s(np.array([1e-4, 1e-3]))
        return rr**2 / (2 * (1 - zz))


if __name__ == '__main__':
    for b in (0.2, 0.3):
        R = Rankine(b, 4001)
        print('b=%.2f c=%.10f half-meridian length S/2=%.6f  nose radius(unzoomed) ~ %s  rho(0)=%.6f' %
              (b, R.c, R.S / 2, R.nose_radius(), R.point_at_s(R.S / 2)[0][0]))
