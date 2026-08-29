# Direct Field Method examples for photonics and plasmonics

This repository is the minimal source-and-input companion to the manuscript
*Computational photonics and plasmonics: A direct physical approach* (manuscript
dated 17 August 2026). It contains the Fortran programs and case inputs needed
for the numerical examples in Figs. 1–6, subject to the explicit limitations
below.

- Version: **0.1.0**
- Primary software developer: **Qiang Sun**
- Recorded sphere-mesh and cube-mesh contributions: **Evert Klaseboer**
- Licence: **BSD 3-Clause License**
- Numerical discretisation: curved six-node quadratic triangular elements
  (Q6)

No computed results, plotted figures, executables, object files, job scripts,
or machine-specific paths are included.

## Manuscript-to-code map

| Manuscript result | Program | Supplied case directories | Main numerical product |
|---|---|---|---|
| Fig. 1: PEC sphere and cube | `code/dfm_field_solver` | `inputs/fig01_pec_sphere_cube/{sphere,cube,cube_inset}` | Surface fields and planar field samples |
| Fig. 2: three closely spaced PEC rods | `code/dfm_field_solver` | `inputs/fig02_pec_rods/{ka0001,ka5,ka50}` | Common DFM field/RCS route |
| Fig. 3: two dielectric lenses | `code/dfm_field_solver` | `inputs/fig03_dielectric_lenses/{ka1,ka10,ka50}` | Planar total electric field |
| Fig. 4: eccentric Au-core/Ag-shell spectrum | `code/au_ag_spectrum` | `inputs/fig04_au_ag_spectrum/core_shell` | Scattering, extinction, and absorption cross sections |
| Fig. 5: Au-core/Ag-shell surface derivative | `code/au_ag_surface_gradient` | `inputs/fig05_au_ag_surface_gradient/{level08,level10}` | Two-sided complex surface traces and Q6 connectivity |
| Fig. 6: optical force and torque | `code/optical_force` | `inputs/fig06_optical_force/pill` | Force and torque versus relative orientation |

## Requirements

- a Fortran compiler with OpenMP support;
- LAPACK and BLAS, or Apple's Accelerate framework on macOS;
- GNU Make;
- enough memory for dense complex linear systems.

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
mkdir -p ../dfm-runs/fig02-ka5
cp -R inputs/fig02_pec_rods/ka5/. ../dfm-runs/fig02-ka5/
(
  cd ../dfm-runs/fig02-ka5
  "$repo/code/dfm_field_solver/.build/field_solver"
)
```

Use the corresponding executable for the other case directories. The Au/Ag
surface-trace executable accepts either no argument or the optional lowercase
compatibility token `run`; both forms perform the same calculation:

```sh
mkdir -p ../dfm-runs/fig05-level08
cp -R inputs/fig05_au_ag_surface_gradient/level08/. ../dfm-runs/fig05-level08/
(
  cd ../dfm-runs/fig05-level08
  "$repo/code/au_ag_surface_gradient/.build/au_ag_surface_gradient" run
)
```

All supplied cases include the conventional files `Input_Geom.dat`,
`Input_Phys_EM.dat`, and `Input_Source_EM.dat`. These examples declare no
impressed sources, so only the geometry and physics files are read. A modified
case that declares impressed sources will also read `Input_Source_EM.dat`.

Cases containing `Prtl_*.msh` must keep those mesh files beside the inputs.
The Fig. 4 case additionally requires `Input_nk_Ag.dat` and `Input_nk_Au.dat`.

## Principal outputs

| Program | Principal output files |
|---|---|
| `field_solver` | `Rslt_SurfCal1_EM.dat`, `Rslt_SurfCal2_EM.dat`, and `Rslt_SurfCal3_EM.dat` contain incident, scattered, and total surface fields. `Rslt_Dmn1_plot.dat`, `Rslt_Dmn2_plot.dat`, and `Rslt_Dmn3_plot.dat` contain planar field samples when domain output is enabled. The Fig. 2 inputs also write `Rslt_RCS_xy.dat`, `Rslt_RCS_yz.dat`, and `Rslt_RCS_zx.dat`. |
| `au_ag_spectrum` | `Rslt_SPR_XS.dat` contains wavelength, the closed-surface incident-flux residual `inc_xs` (ideally zero), and the scattering, extinction, and absorption cross sections for each Au-core radius. |
| `au_ag_surface_gradient` | `surface_traces_both_sides.dat` contains the complex total `E`, `H`, `dE/dn`, and `dH/dn` traces on both interfaces. `surface_elements_q6.dat` contains interface-labelled Q6 connectivity; `cross_sections.dat` and `surface_trace_audit_summary.txt` provide cross sections and boundary diagnostics. |
| `optical_force` | `Rslt_FrcTrq.dat` contains angle followed by the three force and three torque components for each particle. For the supplied unit-amplitude field, the normalisations are `F/(epsilon_0 a^2)` and `N/(epsilon_0 a^3)`. |

`Rslt_PrcdSmm.dat` is a short run record written by each driver. Generated
geometry scratch files named `Prtl_Orgnl.*` are internal and are not user
inputs.

## Figure-specific notes

### Figure 1

The `cube` and `cube_inset` directories intentionally contain different cube
meshes and must both be retained. The sphere is generated internally; the cube
cases read `Prtl_0001.msh`.

### Figure 2

The directories represent `k_out a = 0.001`, `5`, and `50`. Each contains the
same three-rod dimensions with a gap `g = 0.005a`. 

### Figure 3

The three directories represent `k_out a = 1`, `10`, and `50` for two Q6
dielectric oblate spheroids with refractive-index ratio 1.5.

### Figure 4

The retained driver evaluates an internal superset: Au-core radii from 5 to
45 nm in 5 nm steps and wavelengths from 300 to 900 nm in 10 nm steps. The
manuscript plots the 10–45 nm and 300–800 nm subset. `Rslt_SPR_XS.dat` appends
to an existing file, so use a fresh case copy. The homogeneous-sphere Mie
curves in Fig. 4(b) require an independently licensed analytic Mie
implementation; the same Au and Ag optical constants are available in the
retained `Input_nk_*.dat` tables.

### Figure 5

`level08` and `level10` contain the same eccentric Au-core/Ag-shell geometry at
480 nm with two Q6 mesh resolutions. The executable exports the boundary data
needed for the first-order finite-distance field reconstruction. Plotting and
the `delta = 2 nm` reconstruction are postprocessing steps and are not included
in this source-and-input-only package.

### Figure 6

The driver generates 73 orientations from 0 to 360 degrees in 5-degree steps;
the manuscript plots the 0–180 degree subset. Both imported particle meshes
are required. Version 0.1.0 uses the same 16-point symmetric triangle
quadrature rule as the other examples. The figure-specific quadrature check is
summarised under Validation status.

## Numerical scope

The released cases use Q6 surface elements, a fixed six-point Gauss-Legendre
line rule, and a fixed 16-point degree-eight symmetric triangle rule. Unused
12-, 25-, and 49-point tables and selected unrelated solver branches have been
removed.
The triangle-rule provenance and other numerical-method attributions are listed
in `THIRD_PARTY_NOTICES.md`.

This repository supplies numerical solvers and inputs, not a pixel-identical
figure-production pipeline. It intentionally excludes historical results,
plot-composition scripts, HPC launch files, and machine-specific project paths.

## Validation status

The four retained executable targets were built with Intel `ifx` 2024.0.2 and
Intel MKL 2024.0.0 on NCI Gadi, where representative cases spanning Figs. 1–6
completed.

For the Fig. 6 quadrature check, the six plotted components (`Fx_1`, `Fz_1`,
`Ny_1`, `Fx_2`, `Fz_2`, and `Ny_2`) were compared over 0–180 degrees (37
samples) between the retained 16-point rule and an archived 25-point result.
The global relative L2 difference was `3.23e-4`. The largest component-wise,
peak-normalised maximum difference was `1.83e-3` for `Fx_2`; the other five
values were at most `5.32e-4`. Unplotted components that are near zero by
symmetry can show much larger normalised relative differences because their
reference peaks are also near zero, so they were not used for this
figure-specific assessment.

These checks establish build and run readiness and show that the quantities
plotted in Fig. 6 change only slightly between the two quadrature rules in this
test. They are not a formal quadrature-convergence study or a claim of bitwise
identity across compilers or platforms.

## Developers, provenance, and licence

See `DEVELOPERS.md` for authorship and contribution information and
`THIRD_PARTY_NOTICES.md` for numerical-method and input-data provenance.

Copyright (c) 2026 Qiang Sun. The source code and accompanying inputs are
available under the `BSD 3-Clause License`; see `LICENSE` for the binding
terms. The licence permits use, modification, and redistribution in source or
binary form, including for commercial purposes, provided its conditions are
met. In particular, the copyright and licence notices must be retained, and
the names of the copyright holder and contributors may not be used to endorse
or promote derived products without specific prior written permission.

This package is open-source software under an OSI-approved licence. The
licence does not override the separate provenance and dependency notices in
`THIRD_PARTY_NOTICES.md`.

## Manuscript citation

E. Klaseboer, D. Y. C. Chan, A. J. Yuffa, A. C. Yucel, and Q. Sun,
“Computational photonics and plasmonics: A direct physical approach,”
manuscript (2026).
