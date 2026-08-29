
! SPDX-FileCopyrightText: 2026 Qiang Sun
! SPDX-License-Identifier: BSD-3-Clause

MODULE EM_SurfCal_GlobalData

    IMPLICIT NONE

! Shared electromagnetic input, material, surface-field, and post-processing data.
!
! Naming and field conventions used throughout the EM modules:
!   ex / in : trace in the parent (exterior) domain / domain owned by this surface.
!             For a top-level surface the parent is domain 0; for a nested surface
!             it is corelnkshell(surface).  These labels do not specify the direction
!             of the geometric normal.
!   E1, H1  : prescribed incident/source field.
!   E2, H2  : solved scattered or correction field.
!   E3, H3  : total field, E3 = E1 + E2 and H3 = H1 + H2.
!   x,y,z   : Cartesian component of the field immediately preceding the suffix.
!   dnn     : derivative of that Cartesian component along n; dt1 and dt2 are
!             derivatives along the two surface tangents.
!
! The vectors n, t1, and t2 come from Geom_GlobalData.  In particular, n follows
! NrmlInOut in Input_Geom.dat, and the same stored n is used for ex and in traces.
! Complex fields are peak phasors with implicit exp(-i*omega*t); incident travelling
! waves therefore contain exp(+i*k.r).  No automatic unit conversion is performed.

! Optional pulse/FFT controls; ignored by workflows with FFTpulse_EM = 0.
    INTEGER :: FFTpulse_EM, FFTnmbrFrq_EM, FFTnmbrBIM_EM, FFTBIMCal_EM
    CHARACTER (LEN=4) :: FFTpulsetype_EM
    DOUBLE PRECISION :: FFTWidth_EM, FFTnmbrOsPulse_EM, FFTalfa_EM, FFTphsPulse_EM

! Excitation definition and derived free-space wavelength, wavenumber, and frequency.
    CHARACTER (LEN=3) :: excitetype_EM
    CHARACTER (LEN=1) :: poltype_EM
    CHARACTER (LEN=1) :: wlorwn_EM
    INTEGER :: incOrder_EM
    DOUBLE PRECISION :: incFeature_EM
    DOUBLE PRECISION :: incFieldx_EM,incFieldy_EM,incFieldz_EM,incFieldmdl_EM, &
    &                   inckx_EM,incky_EM,inckz_EM
    DOUBLE PRECISION :: phase_EM

! Unbounded exterior (domain 0) material; n,k pairs form relative eps and mu by squaring.
    DOUBLE PRECISION :: exepsn_EM, exepsk_EM, exmiun_EM, exmiuk_EM, prismn_EM, prismk_EM
    COMPLEX(KIND=KIND(1.0D0)) :: exeps_EM, exmiu_EM, exk_EM

    DOUBLE PRECISION :: vcmwl_EM, vcmwn_EM, AngFrqnc_EM

! Per-surface boundary data and material in the domain owned by each surface.
    CHARACTER (LEN = 3), ALLOCATABLE, DIMENSION (:) :: BCType_EM
    CHARACTER (LEN = 1), ALLOCATABLE, DIMENSION (:) :: BCRead_EM
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: BCValue_EM

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: inepsn_EM, inepsk_EM, &
    &                                               inmiun_EM, inmiuk_EM
    COMPLEX(KIND=KIND(1.0D0)), ALLOCATABLE, DIMENSION (:) :: ineps_EM, inmiu_EM, ink_EM

! Complex E/H surface traces and their directional derivatives at mesh nodes.
    COMPLEX(KIND=KIND(1.0D0)), ALLOCATABLE, DIMENSION (:) :: &
    &   exE1x_EM, exE1y_EM, exE1z_EM, inE1x_EM, inE1y_EM, inE1z_EM, &
    &   exE1xdnn_EM, exE1ydnn_EM, exE1zdnn_EM, inE1xdnn_EM, inE1ydnn_EM, inE1zdnn_EM, &
    &   exE1xdt1_EM, exE1ydt1_EM, exE1zdt1_EM, inE1xdt1_EM, inE1ydt1_EM, inE1zdt1_EM, &
    &   exE1xdt2_EM, exE1ydt2_EM, exE1zdt2_EM, inE1xdt2_EM, inE1ydt2_EM, inE1zdt2_EM, &
    &   exE2x_EM, exE2y_EM, exE2z_EM, inE2x_EM, inE2y_EM, inE2z_EM, &
    &   exE2xdnn_EM, exE2ydnn_EM, exE2zdnn_EM, inE2xdnn_EM, inE2ydnn_EM, inE2zdnn_EM, &
    &   exE3x_EM, exE3y_EM, exE3z_EM, inE3x_EM, inE3y_EM, inE3z_EM, &
    &   exE3xdnn_EM, exE3ydnn_EM, exE3zdnn_EM, inE3xdnn_EM, inE3ydnn_EM, inE3zdnn_EM, &
    &   &
    &   exH1x_EM, exH1y_EM, exH1z_EM, inH1x_EM, inH1y_EM, inH1z_EM, &
    &   exH1xdnn_EM, exH1ydnn_EM, exH1zdnn_EM, inH1xdnn_EM, inH1ydnn_EM, inH1zdnn_EM, &
    &   exH1xdt1_EM, exH1ydt1_EM, exH1zdt1_EM, inH1xdt1_EM, inH1ydt1_EM, inH1zdt1_EM, &
    &   exH1xdt2_EM, exH1ydt2_EM, exH1zdt2_EM, inH1xdt2_EM, inH1ydt2_EM, inH1zdt2_EM, &
    &   exH2x_EM, exH2y_EM, exH2z_EM, inH2x_EM, inH2y_EM, inH2z_EM, &
    &   exH2xdnn_EM, exH2ydnn_EM, exH2zdnn_EM, inH2xdnn_EM, inH2ydnn_EM, inH2zdnn_EM, &
    &   exH3x_EM, exH3y_EM, exH3z_EM, inH3x_EM, inH3y_EM, inH3z_EM, &
    &   exH3xdnn_EM, exH3ydnn_EM, exH3zdnn_EM, inH3xdnn_EM, inH3ydnn_EM, inH3zdnn_EM

! Source counts, contiguous source-index ranges, locations, strengths, and polarisations.
    INTEGER :: ttlnmbrsrc_EM, exnmbrsrc_EM
    INTEGER, ALLOCATABLE, DIMENSION (:) :: nmbrsrc_EM, srcstaID_EM, srcendID_EM
    INTEGER, ALLOCATABLE, DIMENSION (:) :: srcType_EM
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: &
    &   srcStrength_EM, xsrc_EM, ysrc_EM, zsrc_EM, polxsrc_EM, polysrc_EM, polzsrc_EM

! Per-surface force, torque, and surface-integral results.
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: &
    &   Frcx_EM, Frcy_EM, Frcz_EM, Trqx_EM, Trqy_EM, Trqz_EM, &
    &   FrcIncx_EM, FrcIncy_EM, FrcIncz_EM, TrqIncx_EM, TrqIncy_EM, TrqIncz_EM, &
    &   FrcScax_EM, FrcScay_EM, FrcScaz_EM, TrqScax_EM, TrqScay_EM, TrqScaz_EM, &
    &   FrcExtx_EM, FrcExty_EM, FrcExtz_EM, TrqExtx_EM, TrqExty_EM, TrqExtz_EM, &
    &   FrcTotx_EM, FrcToty_EM, FrcTotz_EM, TrqTotx_EM, TrqToty_EM, TrqTotz_EM, &
    &   FrcElcx_EM, FrcElcy_EM, FrcElcz_EM, TrqElcx_EM, TrqElcy_EM, TrqElcz_EM, &
    &   FrcMagx_EM, FrcMagy_EM, FrcMagz_EM, TrqMagx_EM, TrqMagy_EM, TrqMagz_EM, &
    &   FrcTtlx_EM, FrcTtly_EM, FrcTtlz_EM, TrqTtlx_EM, TrqTtly_EM, TrqTtlz_EM, &
    &   exsurfQ1_EM, exsurfQ2_EM, exsurfQ3_EM, insurfQ1_EM, insurfQ2_EM, insurfQ3_EM

! Cross-section sampling sphere: enable flag, radius, centre, and accumulated results.
    INTEGER :: on_xs
    DOUBLE PRECISION :: size_xs,xcen_xs,ycen_xs,zcen_xs,inc_xs,sct_xs,exc_xs,abs_xs

! Surface-enhancement controls and per-surface norms.
    INTEGER :: on_enhance
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: L1SurfIntg,L2SurfIntg,L4SurfIntg, &
    &                                               L1tip,L2tip,L4tip

END MODULE
