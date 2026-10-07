# Direct Field Method examples for photonics and plasmonics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23205796.svg)](https://doi.org/10.5281/zenodo.23205796)

This repository is the minimal source, input, and figure-data companion to the
manuscript *Computational photonics and plasmonics: A direct physical approach*
(manuscript submitted to *Journal of Physics: Photonics*, October 2026). It
contains the Fortran programs and case inputs needed for the numerical examples
in Figs. 1–7 (including the gap-graded surface meshes of Fig. 3 and the Python
script that generates them), a short Python exercise that compares a computed result with
Mie theory, the numerical results and Python plotting scripts for Figs. 2, 3,
5 and 7 (`figure_data/`), and a Python script that computes the fields of
perfect conductors from the surface current of the magnetic-field solution, as
used for the low-frequency case of Fig. 3 (`postprocessing/`), subject to the
explicit limitations below.

- Version: **0.1.0**, archived at Zenodo: [doi:10.5281/zenodo.23205796](https://doi.org/10.5281/zenodo.23205796)
  (all versions: [doi:10.5281/zenodo.23205795](https://doi.org/10.5281/zenodo.23205795))
- Primary software developer: **Qiang Sun**
- Recorded sphere-mesh and cube-mesh contributions: **Evert Klaseboer**
- Licence: **BSD 3-Clause License**
- Numerical discretisation: curved six-node quadratic triangular elements
  (Q6)

Apart from the results plotted in Figs. 2, 3, 5 and 7 (`figure_data/`), no
computed results are included, and no plotted figures, executables, object
files, job scripts, or machine-specific paths.

## Manuscript-to-code map

| Manuscript result | Program | Supplied case directories | Main numerical product |
|---|---|---|---|
| Fig. 1: PEC sphere and cube | `code/dfm_field_solver` | `inputs/fig01_pec_sphere_cube/{sphere,cube,cube_inset}` | Surface fields and planar field samples |
| Fig. 2: PEC sphere against Mie theory and mesh convergence | `code/dfm_field_solver` with `exercises/pec_sphere_mie.py` | (a) `inputs/fig01_pec_sphere_cube/sphere`; (b) `exercises/pec_sphere_convergence_ka4p5`, mesh levels 4–14 | Error of the surface field and of the scattered field at `(0, 0, 1.2a)` relative to Mie theory |
| Fig. 3: three closely spaced PEC rods | `code/dfm_field_solver`; for `k_out a = 1e-3` also `postprocessing/pec_current_route.py` | `inputs/fig02_pec_rods/{ka0001,ka5,ka50}` (gap-graded meshes) | Planar field samples and bistatic RCS cuts; at `k_out a = 1e-3` computed from the surface current |
| Fig. 4: two dielectric lenses | `code/dfm_field_solver` | `inputs/fig03_dielectric_lenses/{ka1,ka10,ka50}` | Planar total electric field |
| Fig. 5: eccentric Au-core/Ag-shell spectrum | `code/au_ag_spectrum` | `inputs/fig04_au_ag_spectrum/core_shell` | Scattering, extinction, and absorption cross sections |
| Fig. 6: Au-core/Ag-shell surface derivative and SERS proxy | `code/au_ag_surface_gradient` | `inputs/fig05_au_ag_surface_gradient/{level08,level10}` | Two-sided complex surface traces and Q6 connectivity |
| Fig. 7: optical force and torque | `code/optical_force` | `inputs/fig06_optical_force/pill` | Force and torque versus relative orientation |

The `figNN` prefixes of the input directories, and the figure numbers quoted
in comments and messages of the Fortran sources, follow an earlier numbering
of the manuscript, before the validation figure (now Fig. 2) was added: Figs.
2–6 of the sources are Figs. 3–7 of the manuscript, and Fig. 1 is unchanged.
The files in `figure_data/` use the manuscript numbering.

## Requirements

- a Fortran compiler with OpenMP support;
- LAPACK and BLAS, or Apple's Accelerate framework on macOS;
- GNU Make;
- enough memory for dense complex linear systems;
- for the exercise, Python 3 with NumPy and SciPy (Matplotlib is optional);
  for the plotting scripts in `figure_data/`, Python 3 with NumPy, SciPy, and
  Matplotlib; for `postprocessing/pec_current_route.py` and the mesh generator
  in `inputs/fig02_pec_rods/meshgen/`, Python 3 with NumPy and SciPy.

The Makefiles default to GNU Fortran. The retained programs were also built
with Intel `ifx` and Intel MKL and exercised on NCI Gadi, as described under
Validation status. Compiler and linker settings can be overridden with `FC`,
`FFLAGS`, `CHECKFLAGS`, `LDLIBS`, and `BUILD_DIR`.

## Build

From the repository root:

```sh
make -C code/dfm_field_solver
make -C code/au_ag_spectrum
make -C code/au_ag_surface_gradient
make -C code/optical_force
```

This creates one executable in each `.build` directory:

```text
code/dfm_field_solver/.build/field_solver
code/au_ag_spectrum/.build/au_ag_spectrum
code/au_ag_surface_gradient/.build/au_ag_surface_gradient
code/optical_force/.build/optical_force
```

Use `make -C code/<program> syntax` for a compiler syntax check and
`make -C code/<program> clean` to remove that program's `.build` directory.

## Run a case

Every solver reads fixed filenames from its current working directory. Run in
a fresh, writable copy of a case directory because several legacy `Rslt_*`
writers append to an existing file. For example:

```sh
repo="$(pwd)"
mkdir -p ../dfm-runs/rods-ka5
cp -R inputs/fig02_pec_rods/ka5/. ../dfm-runs/rods-ka5/
(
  cd ../dfm-runs/rods-ka5
  "$repo/code/dfm_field_solver/.build/field_solver"
)
```

Use the corresponding executable for the other case directories. The Au/Ag
surface-trace executable accepts either no argument or the optional lowercase
compatibility token `run`; both forms perform the same calculation:

```sh
mkdir -p ../dfm-runs/gradient-level08
cp -R inputs/fig05_au_ag_surface_gradient/level08/. ../dfm-runs/gradient-level08/
(
  cd ../dfm-runs/gradient-level08
  "$repo/code/au_ag_surface_gradient/.build/au_ag_surface_gradient" run
)
```

All supplied cases include the conventional files `Input_Geom.dat`,
`Input_Phys_EM.dat`, and `Input_Source_EM.dat`. These examples declare no
impressed sources, so only the geometry and physics files are read. A modified
case that declares impressed sources will also read `Input_Source_EM.dat`.

Cases containing `Prtl_*.msh` must keep those mesh files beside the inputs.
The Fig. 5 spectrum case additionally requires `Input_nk_Ag.dat` and
`Input_nk_Au.dat`.

## Principal outputs

| Program | Principal output files |
|---|---|
| `field_solver` | `Rslt_SurfCal1_EM.dat`, `Rslt_SurfCal2_EM.dat`, and `Rslt_SurfCal3_EM.dat` contain incident, scattered, and total surface fields. For each field, the `x`, `y`, and `z` columns are the real parts of the complex phasors (time dependence `exp(-i omega t)`) and the next column is the modulus of the complex vector. `Rslt_Dmn1_plot.dat`, `Rslt_Dmn2_plot.dat`, and `Rslt_Dmn3_plot.dat` contain the incident, scattered, and total fields at planar field samples when domain output is enabled. The Fig. 3 inputs also write `Rslt_RCS_xy.dat`, `Rslt_RCS_yz.dat`, and `Rslt_RCS_zx.dat`. |
| `au_ag_spectrum` | `Rslt_SPR_XS.dat` contains wavelength, the closed-surface incident-flux residual `inc_xs` (ideally zero), and the scattering, extinction, and absorption cross sections for each Au-core radius. |
| `au_ag_surface_gradient` | `surface_traces_both_sides.dat` contains the complex total `E`, `H`, `dE/dn`, and `dH/dn` traces on both interfaces. `surface_elements_q6.dat` contains interface-labelled Q6 connectivity; `cross_sections.dat` and `surface_trace_audit_summary.txt` provide cross sections and boundary diagnostics. |
| `optical_force` | `Rslt_FrcTrq.dat` contains angle followed by the three force and three torque components for each particle. For the supplied unit-amplitude field, the normalisations are `F/(epsilon_0 a^2)` and `N/(epsilon_0 a^3)`. Each torque is taken about the centre-of-mass point given for that particle in `Input_Geom.dat`, which is the origin for both rods in the supplied input. |

`Rslt_PrcdSmm.dat` is a short run record written by each driver. Generated
geometry scratch files named `Prtl_Orgnl.*` are internal and are not user
inputs.

## Exercise: the PEC sphere against Mie theory

`exercises/pec_sphere_mie.py` repeats the two checks of Fig. 2. It needs
only NumPy and SciPy, reads `k_out a`, the polarisation, and the direction of
incidence from `Input_Phys_EM.dat` in the run directory, and evaluates the
exact Mie solution for a perfectly conducting sphere at the same points.
`python3 exercises/pec_sphere_mie.py selftest` checks the Mie implementation
itself, and `python3 exercises/pec_sphere_mie.py --help` lists the options.

### Surface field, Fig. 2(a)

Run the sphere case of Fig. 1(a) (`k_out a = 15`, mesh level 8, 2562 nodes)
and compare its total surface field with Mie theory:

```sh
repo="$(pwd)"
make -C code/dfm_field_solver
mkdir -p ../dfm-runs/sphere
cp -R inputs/fig01_pec_sphere_cube/sphere/. ../dfm-runs/sphere/
(
  cd ../dfm-runs/sphere
  "$repo/code/dfm_field_solver/.build/field_solver"
)
python3 exercises/pec_sphere_mie.py surface ../dfm-runs/sphere
```

The script reports a relative L2 error of `|E_tot|` of 1.15% over all 2562
nodes, the value quoted in the manuscript for Fig. 2(a), and a largest nodal
error of 0.0365 `E_0`. On a PEC surface the total field is normal, so
`|E_tot| = |E_n|`, which is read directly from the `|exE3|` column of
`Rslt_SurfCal3_EM.dat`. Add `--plot eplane.png` to plot `|E_tot|` along the
great circle in the E-plane (requires Matplotlib).

To refine the mesh or change the frequency, edit the copied inputs before
running:

- `Input_Geom.dat`, record `'Icshdrl', 8, 1.0d0`: mesh level `L`, which gives
  `N = 40 L^2 + 2` nodes and `20 L^2` Q6 elements (642 nodes for `L = 4`,
  7842 nodes for `L = 14`);
- `Input_Phys_EM.dat`, record `"pwe", "p", 0, 0.0d0, "k", 15.0d0, ...`: the
  value after `"k"` is the wavenumber in inverse mesh units, which equals
  `k_out a` for this unit sphere in vacuum (`"l"` instead of `"k"` gives the
  free-space wavelength). The next six values are the polarisation and the
  direction of propagation.

### Mesh convergence, Fig. 2(b)

`exercises/pec_sphere_convergence_ka4p5` is the sphere case with
`k_out a = 4.5`, an unrotated mesh, and a 3 × 3 planar sample centred on the
observation point `(0, 0, 1.2a)` (first data record of `Input_Phys_EM.dat`),
from which the scattered field at that point is read. The study of Fig. 2(b)
is run by

```sh
make -C code/dfm_field_solver
sh exercises/run_sphere_convergence.sh ../dfm-runs/sphere-convergence
```

which runs mesh levels 4 to 14 (`N = 642` to `7842`) in separate directories
`level04`, ..., `level14` and then calls

```sh
python3 exercises/pec_sphere_mie.py point ../dfm-runs/sphere-convergence/level*
```

to tabulate `|E_sc|` at `(0, 0, 1.2a)`, its error relative to Mie theory, and
the fitted convergence rate. Selected levels can be run by listing them after
the run directory, for example
`sh exercises/run_sphere_convergence.sh ../dfm-runs/sphere-convergence 4 6 8`.
To set up a level by hand, copy the case directory and change the mesh level
in its `'Icshdrl'` record.

The cost grows quickly with the level. On the machine described under
Validation status, level 12 (5762 nodes) needed 11.8 GB of memory; the dense
system for level 14 has 23 526 complex unknowns and a system matrix of 8.9 GB,
and, scaling the level-12 value as `N^2`, needs about 22 GB.

The archived results plotted in Fig. 2(b), supplied in `figure_data/`, were
computed in 2022 with an earlier version of the solver, which differs from
release 0.1.0 in several respects, including its 25-point triangle rule
(release 0.1.0 uses the 16-point rule; see Numerical scope). Levels 4–12 were
rerun with this release. At levels 4–7 the errors are 0.7–6% larger than those
plotted (for example `3.25e-3` instead of `3.07e-3` for `N = 642`); from level
8 on, `|E_sc|` agrees with the archived value to within `1e-6` (relative), and
the errors differ from those plotted by less than 2%. The fitted rate over
levels 4–12 is `N^-2.00`.

## Figure data: Figs. 2, 3, 5 and 7

The directory `figure_data/` contains the numerical results plotted in Figs. 2,
3, 5 and 7 and, for each figure, a Python script that produces the figure from
them:

| Manuscript figure | Script | Results |
|---|---|---|
| Fig. 2: PEC sphere against Mie theory | `make_fig2_pec_validation.py` | `fig2a_pec_sphere_ka15_surface_field.csv`, `fig2b_pec_sphere_ka4p5_convergence.csv` |
| Fig. 3: three closely spaced PEC rods | `make_fig3_three_rods.py` | `fig3_slice_y0_ka{0001,5,50}.csv.gz`, `fig3_rcs_zx.csv` |
| Fig. 5: Au-core/Ag-shell spectra | `make_fig5_au_ag_spectrum.py` | `fig5_Rslt_SPR_XS_archived.dat` |
| Fig. 7: optical force and torque | `make_fig7_optical_force.py` | `fig7_Rslt_FrcTrq_archived.dat` |

The results for Figs. 2, 5 and 7 are archived results computed with earlier
development versions of the programs (2021 and 2022); those for Fig. 3 were
computed with release 0.1.0 and the supplied inputs (October 2026). Their
provenance, file formats, and the relation to release 0.1.0 are described in
`figure_data/README.md`. The scripts compute the Mie reference solutions shown
in Figs. 2 and 5 (perfectly conducting, homogeneous, and coated spheres) and
the torques of Fig. 7 about each rod's own centre. Run them from that
directory:

```sh
cd figure_data
python3 make_fig2_pec_validation.py
python3 make_fig3_three_rods.py
python3 make_fig5_au_ag_spectrum.py
python3 make_fig7_optical_force.py
```

Each script writes the figure (PDF and PNG) and a table of the plotted or
derived values to `figure_data/output/`; copies of these tables, as used for
the manuscript figures, are supplied beside the scripts. The Fig. 5 script
reads the optical constants from `inputs/fig04_au_ag_spectrum/core_shell`.

## Figure-specific notes

### Figure 1

The `cube` and `cube_inset` directories intentionally contain different cube
meshes and must both be retained. The sphere is generated internally; the cube
cases read `Prtl_0001.msh`.

### Figure 2

Figure 2(a) uses the `sphere` case of Fig. 1 and Fig. 2(b) the case
`exercises/pec_sphere_convergence_ka4p5`; both are described under Exercise.
The archived results and the plotting script are in `figure_data/`.

### Figure 3

The directories `ka0001`, `ka5` and `ka50` represent `k_out a = 0.001`, `5`
and `50`; all three are shown in the manuscript. Each contains the same
three-rod geometry, with two gaps `g = 0.005a`. The planar field is sampled in
the plane `y = 0` for `|x| <= 3a`, and the RCS cuts are evaluated at
`r = 1000a` in 1-degree steps.

**Gap-graded meshes.** The narrow gaps must be resolved by the mesh (Scope
and limitations section of the manuscript). Each case therefore imports the
curved six-node meshes `Prtl_0001.msh`, `Prtl_0002.msh` and `Prtl_0003.msh`
(`MeshRead = 1`, the third entry of the three `'LgRd'` records of
`Input_Geom.dat`). The element size equals the gap width `g` at the
gap-facing ends of the rods and grows to `0.08a` (`ka0001`, `ka5`) or `0.04a`
(`ka50`) away from them:

| Cases | Nodes (particles 1 / 2 / 3) | Elements | Mesh files |
|---|---|---|---|
| `ka0001`, `ka5` | 11 878 (4254 / 4254 / 3370) | 5936 | 1.7 MB per case |
| `ka50` | 41 798 (15 618 / 15 618 / 10 562), about 6.3 nodes per wavelength | 20 896 | 6.0 MB |

Particles 1 and 2 are the outer rods (length `3a`, centred at `x = -2.505a`
and `+2.505a`) and particle 3 the central rod (length `2a`); the schematic of
Fig. 3 uses the same numbers. The mesh files are written in the local frame of each rod;
the zoom, rotation and translation in `Input_Geom.dat` still apply, so that
the imported rods coincide with the generated `'LgRd'` shapes. The script
`inputs/fig02_pec_rods/meshgen/meshgen_rod.py` (with `rankine.py`)
regenerates the meshes byte for byte; the commands are given in its header.
Setting `MeshRead` back to `0` restores the internally generated uniform
meshes (`'Icshdrl'` level 10, 4002 nodes per rod), which do not resolve the
gaps and were used in earlier drafts of the manuscript; they are not
recommended. The complete runs of `ka5` and `ka0001` need about 30 GB of
memory (3 minutes on 48 cores of NCI Gadi), and that of `ka50` about 365 GB
(2 hours).

**`k_out a = 0.001`: use the magnetic field.** For perfect conductors at very
low frequencies, the net charge of each conductor is only weakly determined by
the electric-field formulation (manuscript, Secs. 2.3 and 7.3). The electric
field written by `field_solver` (the `exE3` columns of `Rslt_SurfCal3_EM.dat`,
`Rslt_Dmn*_plot.dat` and `Rslt_RCS_*.dat`) should therefore not be used for
the fields at very low `k`. For `ka0001`, for example, `Rslt_RCS_zx.dat` is
isotropic at 1.79 `a^2`, whereas the physical radar cross section is about
`3e-10 a^2`. (The RCS radius `r = 1000a` of the supplied `ka0001` input is,
moreover, only `k r = 1`, so that `Rslt_RCS_*.dat` would not be a far-field
quantity there even with a correct field; `pec_current_route.py` computes the
far-field limit.) For the rods of Fig. 3 the magnetic-field solution,
which the solver always computes for a plane wave, does not have this
problem, and the manuscript computes the fields and the radar cross section at
this frequency from the surface current `J = n_out x H` and the surface charge
`div_s J / (i k)`. Because the solver writes only the real parts of the
complex fields, two runs are needed: the case as supplied, and a copy in which
the incident field is multiplied by `i` (incident order 3 with phase `pi/2` in
`Input_Phys_EM.dat`), which writes `Re(i H) = -Im H`. From the repository
root:

```sh
repo="$(pwd)"
mkdir -p ../dfm-runs
cp -R inputs/fig02_pec_rods/ka0001 ../dfm-runs/rods-ka0001
python3 postprocessing/pec_current_route.py prepare inputs/fig02_pec_rods/ka0001 ../dfm-runs/rods-ka0001-ph90
for d in rods-ka0001 rods-ka0001-ph90; do
  ( cd ../dfm-runs/$d && "$repo/code/dfm_field_solver/.build/field_solver" )
done
python3 postprocessing/pec_current_route.py compute ../dfm-runs/rods-ka0001 ../dfm-runs/rods-ka0001-ph90 ../dfm-runs/rods-ka0001-J
```

`prepare` also switches off the planar samples in the second run. `compute`
writes the far-field radar cross section in the three cuts
(`rcs_far_{zx,xy,yz}.csv`), the real part of the total field at the planar
samples of the first run (`slice_zx.csv`), and a summary of the charge and the
dipole moments of each rod (`summary.txt`); the planar samples take about two
minutes. The method is described in the header of the script. The charge
obtained in this way integrates to zero on each rod by construction. The
division by `k` amplifies the discretisation error of the static current, but
this affects only the out-of-phase part of the field (about 4% of the dipole
moment for `ka0001`), which is not used: the slice of Fig. 3(a) is the
in-phase field, and the radar cross section changes by less than 0.01 dB
within 20 dB of its maximum. Applied to the
runs used for the manuscript, it reproduces `figure_data/fig3_slice_y0_ka0001.csv.gz`
and the `k_out a = 0.001` radar cross section in `figure_data/fig3_rcs_zx.csv`
to the precision written. The same two-run procedure applied to `ka5` gives a
radar cross section that agrees with the electric-field result of
`Rslt_RCS_zx.dat` to within 0.23 dB within 20 dB of its maximum.

### Figure 4

The three directories represent `k_out a = 1`, `10`, and `50` for two Q6
dielectric oblate spheroids with refractive-index ratio 1.5.

### Figure 5

The retained driver evaluates an internal superset: Au-core radii from 5 to
45 nm in 5 nm steps and wavelengths from 300 to 900 nm in 10 nm steps. The
manuscript plots the 10–45 nm and 300–800 nm subset. `Rslt_SPR_XS.dat` appends
to an existing file, so use a fresh case copy. The coated-sphere Mie results
in Fig. 5(a) and the homogeneous-sphere Mie curves in Fig. 5(b) are computed
by `figure_data/make_fig5_au_ag_spectrum.py` with the same Au and Ag optical
constants (the retained `Input_nk_*.dat` tables); the exercise script treats
only perfectly conducting spheres. The archived DFM spectra plotted in
Fig. 5(a) are in `figure_data/`.

### Figure 6

`level08` and `level10` contain the same eccentric Au-core/Ag-shell geometry
at 480 nm with two Q6 mesh resolutions, 2562 and 4002 nodes per interface.
Fig. 6 uses `level10`. The executable exports the boundary data needed for the
first-order finite-distance field reconstruction. Plotting and the
`delta = 2 nm` reconstruction are postprocessing steps and are not included
in this package.

### Figure 7

The driver generates 73 orientations from 0 to 360 degrees in 5-degree steps;
the manuscript plots the 0–180 degree subset. Both imported particle meshes
are required. Fig. 7(d) shows the torque on each rod about its own centre,
`x_c = -1.01a` for rod 1 and `+1.01a` for rod 2, whereas `Rslt_FrcTrq.dat`
gives the torque about the origin (see Principal outputs); in the normalised
units, `N_y(own centre) = N_y + x_c F_z`. Version 0.1.0 uses the same 16-point
symmetric triangle quadrature rule as the other examples. The figure-specific
quadrature check is summarised under Validation status. The archived 25-point
result plotted in Fig. 7 and the plotting script, which applies this change of
reference point, are in `figure_data/`.

## Numerical scope

The released cases use Q6 surface elements, a fixed six-point Gauss-Legendre
line rule, and a fixed 16-point degree-eight symmetric triangle rule. Unused
12-, 25-, and 49-point tables and selected unrelated solver branches have been
removed.
The triangle-rule provenance and other numerical-method attributions are listed
in `THIRD_PARTY_NOTICES.md`.

The present implementation stores dense matrices and solves the linear system
by LU factorisation, so memory grows as `N^2` and time as `N^3`. As a guide,
Table 1 of the manuscript lists, for the three rods of Fig. 3, a system matrix
of 20.3 GB and a peak memory of 29.9 GB at `k_out a = 0.001` and 5 (11 878
nodes), and 251.6 GB and 365 GB at `k_out a = 50` (41 798 nodes), and 76.5 GB
and 80.2 GB for the two lenses of Fig. 4.

Apart from the results and plotting scripts for Figs. 2, 3, 5 and 7 in
`figure_data/` and the post-processing script in `postprocessing/`, this
repository supplies numerical solvers and inputs, not a figure-production
pipeline. It intentionally excludes other historical
results, plot-composition scripts for the other figures, HPC launch files, and
machine-specific project paths.

## Validation status

The four retained executable targets were built with Intel `ifx` 2024.0.2 and
Intel MKL 2024.0.0 on NCI Gadi, where representative cases spanning Figs. 1
and 3–7 completed. The results of Fig. 3 were computed with this build of the
unchanged release sources (October 2026), from the supplied cases
`inputs/fig02_pec_rods/{ka0001,ka5,ka50}`; their provenance and checks are
described in `figure_data/README.md` and under Figure 3 above.

For the Fig. 7 quadrature check, the six plotted components (`Fx_1`, `Fz_1`,
`Ny_1`, `Fx_2`, `Fz_2`, and `Ny_2`) were compared over 0–180 degrees (37
samples) between the retained 16-point rule and an archived 25-point result.
The global relative L2 difference was `3.23e-4`. The largest component-wise,
peak-normalised maximum difference was `1.83e-3` for `Fx_2`; the other five
values were at most `5.32e-4`. Unplotted components that are near zero by
symmetry can show much larger normalised relative differences because their
reference peaks are also near zero, so they were not used for this
figure-specific assessment.

These checks establish build and run readiness and show that the quantities
plotted in Fig. 7 change only slightly between the two quadrature rules in this
test. They are not a formal quadrature-convergence study or a claim of bitwise
identity across compilers or platforms.

The four programs also build without compiler warnings with GNU Fortran 16.2
(Homebrew) and Apple Accelerate on macOS (Apple M3 Pro, 12 cores, 36 GB).
With this build, the Fig. 1(a) sphere case reproduced the archived total
surface field to the precision written (largest difference in `|E_tot|` of
`1e-12 E_0`), and hence the 1.15% error of Fig. 2(a), and the convergence case
of Fig. 2(b) was run for mesh levels 4–12, as described under Exercise. The
other cases were not rerun with this build.

## Developers, provenance, and licence

See `DEVELOPERS.md` for authorship and contribution information and
`THIRD_PARTY_NOTICES.md` for numerical-method and input-data provenance.

Copyright (c) 2026 Qiang Sun. The source code, the accompanying inputs, and
the data and scripts in `figure_data/` are available under the
`BSD 3-Clause License`; see `LICENSE` for the binding terms. The licence
permits use, modification, and redistribution in source or binary form,
including for commercial purposes, provided its conditions are met. In
particular, the copyright and licence notices must be retained, and the names
of the copyright holder and contributors may not be used to endorse or promote
derived products without specific prior written permission.

This package is open-source software under an OSI-approved licence. The
licence does not override the separate provenance and dependency notices in
`THIRD_PARTY_NOTICES.md`.

## Citation

E. Klaseboer, D. Y. C. Chan, A. J. Yuffa, A. C. Yucel, and Q. Sun,
“Computational photonics and plasmonics: A direct physical approach,”
manuscript submitted to *Journal of Physics: Photonics* (2026).

Citation metadata for the software are given in `CITATION.cff`, and the
metadata for the Zenodo archive in `.zenodo.json`. The release tag (`v0.1.0`)
and the Zenodo DOI will be added to these files and to this section after the
GitHub release has been archived; until then, please cite the repository URL
and version 0.1.0.
