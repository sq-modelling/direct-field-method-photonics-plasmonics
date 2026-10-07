# Figure data and plotting scripts

This directory contains the numerical results plotted in Figs. 2, 3, 5 and 7
of the manuscript and the Python scripts that produce these figures from them.
The figure numbers in this directory follow the manuscript. (The `figNN`
prefixes of the directories in `inputs/` follow an earlier numbering; see the
top-level `README.md`.)

## Usage

The scripts need Python 3 with NumPy, SciPy, and Matplotlib. Run them from
this directory:

```sh
cd figure_data
python3 make_fig2_pec_validation.py
python3 make_fig3_three_rods.py
python3 make_fig5_au_ag_spectrum.py
python3 make_fig7_optical_force.py
```

Each script prints a few check values, and writes the figure as PDF and PNG
(600 dpi) and a table of plotted or derived values to `output/`, which is
created if necessary. The scripts locate their data relative to their own
location, so they can also be run from elsewhere. Nothing in this directory
other than `output/` is overwritten. The Fig. 5 script also reads the Au and
Ag optical constants from `../inputs/fig04_au_ag_spectrum/core_shell/`.

The figures use the serif font TeX Gyre Termes when Matplotlib can find it,
and otherwise Liberation Serif or DejaVu Serif; mathematical symbols use STIX.
The fallback fonts change only the text layout (for example, in Fig. 5(a) a
label may then touch the legend), not the data.

## Files

| Figure | Script | Data files | Table written to `output/` |
|---|---|---|---|
| Fig. 2 | `make_fig2_pec_validation.py` | `fig2a_pec_sphere_ka15_surface_field.csv`, `fig2b_pec_sphere_ka4p5_convergence.csv` | `fig2b_plotted_data.csv` |
| Fig. 3 | `make_fig3_three_rods.py` | `fig3_slice_y0_ka0001.csv.gz`, `fig3_slice_y0_ka5.csv.gz`, `fig3_slice_y0_ka50.csv.gz`, `fig3_rcs_zx.csv` | `fig3_summary.csv` |
| Fig. 5 | `make_fig5_au_ag_spectrum.py` | `fig5_Rslt_SPR_XS_archived.dat` (and `../inputs/fig04_au_ag_spectrum/core_shell/Input_nk_{Ag,Au}.dat`) | `fig5a_mie_check.csv` |
| Fig. 7 | `make_fig7_optical_force.py` | `fig7_Rslt_FrcTrq_archived.dat` | `fig7_plotted_data.csv` |

Copies of the four tables, as written when the manuscript figures were made,
are included here (`fig2b_plotted_data.csv`, `fig3_summary.csv`,
`fig5a_mie_check.csv`, and `fig7_plotted_data.csv`) for readers who want the
plotted values without running the scripts, and for comparison with a new run.

The archived results for Figs. 2, 5 and 7 were computed with earlier
development versions of the programs in `code/` (2021 and 2022), not with
release 0.1.0. The corresponding cases are supplied in `inputs/` and
`exercises/`, and the differences between the archived results and release
0.1.0 that have been checked are summarised in the top-level `README.md`
(Exercise and Validation status). The results for Fig. 3 were computed in
October 2026 with release 0.1.0 and the supplied inputs (see below).

## Figure 2: PEC sphere against Mie theory

- `fig2a_pec_sphere_ka15_surface_field.csv`: Fig. 2(a). Columns `x/a`, `y/a`,
  `z/a`, and `|E_tot|/E0` at all 2562 surface nodes of the PEC sphere,
  `k_out a = 15`, mesh level 8 (the case of Fig. 1(a); inputs
  `inputs/fig01_pec_sphere_cube/sphere`, run in 2021). The values are the node
  coordinates and the `|exE3|` column of the archived `Rslt_SurfCal3_EM.dat`
  of that run, copied as written; the remaining columns and the element
  connectivity are omitted. Rerunning the supplied case with release 0.1.0
  reproduces these values to `1e-12 E0` (top-level `README.md`, Validation
  status).
- `fig2b_pec_sphere_ka4p5_convergence.csv`: Fig. 2(b). For mesh levels 4–14
  (`N = 642` to `7842` nodes), the number of nodes and elements, `|E_sc|` at
  `(0, 0, 1.2a)` computed for `k_out a = 4.5`, and the time of the complete run
  in seconds. The values were collected from the archived summary files of the
  2022 convergence study (one file per mesh level for `|E_sc|`, and the
  run-time column, labelled `CPU (s)`, of the summary table). The 642-node run
  is recorded as 1 s and is labelled `<1 s` in the figure. The release case
  `exercises/pec_sphere_convergence_ka4p5` repeats this study; the differences
  from the archived values are given in the top-level `README.md` (Mesh
  convergence, Fig. 2(b)).

The script evaluates the exact Mie solution for a perfectly conducting sphere
(Bohren and Huffman 1983, ch. 4) at the same points. It prints the relative
`L2` error of `|E_tot|` over all 2562 nodes (1.15%), the largest nodal error
(0.0365 `E0`), the number of plotted E-plane nodes with `|y| < 0.03a` (64), and,
for panel (b), the Mie value of `|E_sc|` at `(0, 0, 1.2a)`, the relative errors
at 642 and 7842 nodes (`3.07e-3` and `2.15e-5`), and the fitted slope (−1.99).

## Figure 3: three closely spaced PEC rods

The results were computed in October 2026 with the Fortran sources of release
0.1.0, unchanged (built with Intel `ifx` 2024.0.2 and MKL 2024.0.0 and run on
48 cores of NCI Gadi), from the case directories `inputs/fig02_pec_rods/ka0001`,
`ka5` and `ka50` exactly as supplied. These cases import the gap-graded surface
meshes `Prtl_000{1,2,3}.msh` (see the top-level `README.md`, Figure 3).

- `fig3_slice_y0_ka5.csv.gz`, `fig3_slice_y0_ka50.csv.gz`: Figs. 3(b) and (c).
  The `plane - zx` zone of `Rslt_Dmn3_plot.dat` (total field of the
  electric-field solution) of the runs of `ka5` and `ka50`: 120 x 200 samples
  in the plane `y = 0`, `-3a <= x <= 3a`, `-2.5a <= z <= 7.5a`. Columns `x/a`,
  `z/a`, `pos` (from `RsltPos.dat`: 0 outside the rods, otherwise the number,
  in `Input_Geom.dat`, of the particle that contains the point: 1 and 2 are the
  outer rods at `x = -2.505a` and `+2.505a`, 3 the central rod, as labelled in
  the schematic of Fig. 3), and `Re(E_x^tot/E0)` (0 inside the rods),
  copied with 12 significant digits, i.e. as written by the solver.
- `fig3_slice_y0_ka0001.csv.gz`: Fig. 3(a), `k_out a = 1e-3`, on the same
  samples. At this frequency the electric field written by the solver is not
  used (top-level `README.md`, Figure 3). The field was computed from the
  surface current `J = n_out x H` of the magnetic-field solution with
  `../postprocessing/pec_current_route.py`, from the run of `ka0001` and a
  second run with the incident field multiplied by `i`. Same columns;
  `Re(E_x^tot/E0)` is the in-phase total field. At the samples more than
  `0.02a` from the rods it was checked against an electrostatic calculation
  for the uncharged rods on the same mesh (manuscript, Sec. 3.2.4): the
  differences are at most `0.32 E0` between `0.02a` and `0.1a` from the rods
  and below 1% of the largest field beyond `0.1a`. The largest value at these
  samples is 19.0, at `(x, z) = (+-0.983a, -0.138a)`, `0.027a` from the rod
  ends beside the gaps. The 129 samples within `0.02a` of the rods were not checked
  against this reference; they reach `Re(E_x^tot/E0) = 42.5` at
  `(+-0.983a, -0.088a)`, `0.0025a` from the ends of the central rod. All of
  these values lie far beyond the colour scale of the figure.
- `fig3_rcs_zx.csv`: Figs. 3(d)–(f) and two further curves quoted in the
  manuscript text. Angle `theta_deg` (0 to 360 degrees in 1-degree steps;
  `x = r cos(theta)`, `z = r sin(theta)`, so that 90 degrees is the forward
  direction) and `sigma/a^2`:
  - `sigma_over_a2_ka0001_far`: Fig. 3(d), far-field limit computed from the
    surface current with `pec_current_route.py` (file `rcs_far_zx.csv`);
  - `sigma_over_a2_ka5_r1000`, `sigma_over_a2_ka50_r1000`: Figs. 3(e) and (f),
    `Rslt_RCS_zx.dat` of the runs of `ka5` and `ka50` (electric field,
    `r = 1000a`);
  - `sigma_over_a2_ka50_r1e5` (not plotted): `ka50` with the RCS radius set to
    `1e5 a` in `Input_Phys_EM.dat` (far-field limit of the same evaluator),
    used in the text for the comparison of `r = 1000a` with the far field;
  - `sigma_over_a2_ka5_current_far` (not plotted): `ka5`, far field computed
    from the surface current with `pec_current_route.py`, used in the text for
    the comparison with the electric-field result (0.23 dB at most within
    20 dB of the maximum).

The script masks the samples inside the rods, plots `Re(E_x^tot/E0)` on
separate colour scales for (a) (0 to 2) and for (b), (c) (−2 to 2), and plots each radar cross section as
`10 log10(sigma/a^2)` over a 60 dB range given by the scale below each polar
panel. `fig3_summary.csv` lists, for each `k_out a`, the maximum of the radar
cross section and its direction, the forward and backward values (dB), and the
range of `Re(E_x^tot/E0)` over the samples outside the rods (for
`k_out a = 1e-3` its maximum, 42.46, is at the unchecked samples `0.0025a` from
the rods described above). The script prints the same values: maxima of
−95.62, 25.37 and 42.37 dB, at 270, 90 and 90 degrees, for `k_out a = 1e-3`, 5
and 50.

## Figure 5: Au-core/Ag-shell absorption spectra

- `fig5_Rslt_SPR_XS_archived.dat`: the output file `Rslt_SPR_XS.dat` of the
  archived 2021 run of the case `inputs/fig04_au_ag_spectrum/core_shell`
  (identical inputs and optical-constant tables), copied unchanged. One zone
  per Au-core radius, 5 to 45 nm in 5 nm steps; columns: free-space wavelength
  (nm, 300 to 900 nm in 10 nm steps), the closed-surface incident-flux residual
  `inc_xs` (ideally zero), and the scattering, extinction, and absorption cross
  sections in µm². Fig. 5(a) plots the absorption cross sections for core
  radii 10–45 nm between 300 and 800 nm.

The script computes the coated-sphere Mie solution (Aden and Kerker 1951;
Bohren and Huffman 1983, Sec. 8.1) for the concentric case `a_Au = 45 nm`
(open circles in Fig. 5(a)) and the homogeneous-sphere Mie solution for Ag and
Au spheres of radius 50 nm (Fig. 5(b)), with the same optical-constant tables,
linearly interpolated, and a refractive index of 1.33 for water.
`fig5a_mie_check.csv` lists the DFM and Mie absorption cross sections and their
relative difference at the 20 nm spaced wavelengths plotted as circles. The
script prints the largest relative difference at these wavelengths (1.72%
between 300 and 600 nm, 9.16% between 620 and 800 nm); the manuscript's value
of 2.8% for 300–600 nm refers to all 10 nm spaced wavelengths of the archived
run. It also prints the Mie absorption peaks of the homogeneous spheres
(Ag: 0.0193 µm² at 402 nm; Au: 0.0214 µm² at 552 nm).

## Figure 7: optical force and torque on two rods

- `fig7_Rslt_FrcTrq_archived.dat`: the output file `Rslt_FrcTrq.dat` of the
  archived 2021 run of the case `inputs/fig06_optical_force/pill`, copied
  unchanged. The run used the 25-point triangle rule of the development
  version; its comparison with release 0.1.0 (16-point rule) is described in
  the top-level `README.md` (Validation status). Columns: angle `theta` in
  degrees (0 to 360 in 5 degree steps), then `F_x`, `F_y`, `F_z`, `N_x`, `N_y`,
  `N_z` for rod 1 and for rod 2, normalised by `eps0 |E0|^2 a^2` (force) and
  `eps0 |E0|^2 a^3` (torque), `a = 0.5 µm`. The torques are about the origin.

The script plots the 0–180 degree subset and converts the torque to each
rod's own centre, `N_y(own centre) = N_y + x_c F_z` with `x_c = -1.01a` for
rod 1 and `+1.01a` for rod 2. `fig7_plotted_data.csv` lists the plotted forces,
the torques about each rod's own centre, and the torques about the origin. The
script prints the range and the zero crossings of the torque on rod 1 (from
−0.0246 to 0.0217, crossing zero at 87.0 and 170.0 degrees), the range of the
torque on rod 2, the extrema of the forces, and the largest difference between
the results at `theta` and `theta + 180°`, which describe the same geometry.

## Reproduction check

The scripts for Figs. 2, 5 and 7 were run in this directory with Python
3.12, NumPy 2.1, SciPy 1.17, and Matplotlib 3.9. They printed the values given
above, and the three tables written to `output/` were identical to the copies
supplied here. Rasterised at 150 dpi, the regenerated PDF files were
pixel-identical to the figures of the manuscript: Fig. 2 with the fallback font
DejaVu Serif, with which the manuscript figure was made, and Figs. 5 and 7 with
TeX Gyre Termes installed. The Fig. 3 script was run in the same way with
Matplotlib 3.9.2 and TeX Gyre Termes (which it also finds in a TeX
installation through `kpsewhich`); `output/fig3_summary.csv` was identical to
the copy supplied here, and the PDF, rasterised at 150 dpi, was pixel-identical
to the manuscript figure. Running `../postprocessing/pec_current_route.py` on
the two Gadi runs of `ka0001` reproduced `fig3_slice_y0_ka0001.csv.gz` to
`5e-10 E0` and `sigma_over_a2_ka0001_far` to `5e-11` (relative), i.e. to the
precision written.

## Licence

The scripts and data files in this directory are part of this repository and
are distributed under the BSD 3-Clause License; see `../LICENSE`.
