# Third-party notices and provenance

This directory is an open-source research-code package. Material for which
Qiang Sun is the licensor is made available under the BSD 3-Clause License in
`LICENSE`. That licence does not override third-party terms, ownership, or
attributions.

## External numerical libraries (not redistributed)

The Fortran programs call LAPACK routines, including `zgesv` and `dgesv`, and use OpenMP. The supplied Makefiles select Apple's Accelerate framework on macOS and compatible LAPACK/BLAS libraries on other systems. This repository does not redistribute those libraries, the OpenMP runtime, or compiler source or binaries. Users must obtain the dependencies under their respective terms.

The optional Python exercise in `exercises/` uses NumPy and SciPy, and Matplotlib for its optional plot. The plotting scripts in `figure_data/` use NumPy, SciPy, and Matplotlib. The post-processing script in `postprocessing/` and the mesh generator in `inputs/fig02_pec_rods/meshgen/` use NumPy and SciPy. These packages are not redistributed.

## Numerical-method and input-data attributions

- The retained 16-point, degree-eight symmetric triangle quadrature rule is tabulated by L. Zhang, T. Cui, and H. Liu, [“A set of symmetric quadrature rules on triangles and tetrahedra”](https://lsec.cc.ac.cn/~tcui/myinfo/paper/quad.pdf), *Journal of Computational Mathematics* **27** (2009), 89–96. The source stores weights scaled for a reference triangle of area one half. No Witherden–Vincent quadrature table is retained.
- `postprocessing/pec_current_route.py` uses the same 16-point triangle rule. The mesh generator `inputs/fig02_pec_rods/meshgen/meshgen_rod.py` uses the classical 7-point, degree-five triangle rule (Radon 1948) only to check the area and volume of the generated meshes; the rod shape is the Rankine ovoid of the solver's `'LgRd'` generator. These are method citations; no third-party code is included.
- The author-written finite-difference coefficient generator follows the recurrence described by B. Fornberg, “Generation of finite difference formulas on arbitrarily spaced grids,” *Mathematics of Computation* **51** (1988), 699–706.
- Surface-normal weighting follows methods cited in the source comments to N. Max, *Journal of Graphics Tools* **4**(2), and Chen and Wu, *Computer Aided Geometric Design* **21** (2004), 447–458. These are method citations, not bundled third-party libraries.
- The Au and Ag refractive-index inputs were generated from the Brendel–Bormann model of A. D. Rakić, A. B. Djurišić, J. M. Elazar, and M. L. Majewski, “Optical properties of metallic films for vertical-cavity optoelectronic devices,” *Applied Optics* **37** (1998), 5271–5283. Only the numerical case inputs are included; the legacy MATLAB generator is not included.
- The exercise script `exercises/pec_sphere_mie.py` evaluates the Mie series for a perfectly conducting sphere as given by C. F. Bohren and D. R. Huffman, *Absorption and Scattering of Light by Small Particles* (Wiley, 1983), chapter 4. This is a method citation; no third-party code is included.
- The plotting scripts in `figure_data/` evaluate the same Mie series for a perfectly conducting sphere (Fig. 2) and the Mie series for homogeneous and coated spheres (Fig. 5) as given by Bohren and Huffman (1983), chapter 4 and section 8.1, after A. L. Aden and M. Kerker, “Scattering of electromagnetic waves from two concentric spheres,” *Journal of Applied Physics* **22** (1951), 1242–1246. These are method citations; no third-party code is included.

## Legacy material excluded from this package

The working research tree contained optional or unused routines with separate or unresolved redistribution terms. They are deliberately excluded here: J.-P. Moreau/Numerical Recipes-derived LU and FFT routines; M. A. Botchev and D. R. Fokkema BiCGStab routines carrying a no-resale condition; Alan Miller/RANLIB random-number routines; SHTOOLS-derived Gaussian-random-sphere code; AMOS/Sandia special-function routines; and the legacy MATLAB optical-constant generator attributed to Collin Meierbachtol and modified by Qiang Sun.

These references identify mathematical methods and numerical input sources; they do not create or replace a licence grant.
