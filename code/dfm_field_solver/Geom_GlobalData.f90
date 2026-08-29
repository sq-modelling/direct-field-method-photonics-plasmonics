
! SPDX-FileCopyrightText: 2026 Qiang Sun
! SPDX-License-Identifier: BSD-3-Clause

MODULE Geom_GlobalData

    IMPLICIT NONE

! Shared case controls, quadrature tables, mesh topology, geometry, local frames,
! and differential-geometry arrays. Objects are concatenated into global node and
! element arrays; ndstaID/ndendID and elstaID/elendID delimit each object.
! corelnkshell=0 identifies a top-level surface; otherwise it stores the parent.

    INTEGER :: csStartID, csEndID, csStepNumber
    DOUBLE PRECISION :: csStartValue, csStepSize

! Observation-domain sampling controls.
    INTEGER :: on_fldcal
    CHARACTER (LEN = 2) :: Postxyz_RMpln
    DOUBLE PRECISION :: PostOfstx_RMpln,PostOfsty_RMpln,PostOfstz_RMpln
    DOUBLE PRECISION :: PostEgSZx_RMpln,PostEgSZy_RMpln,PostEgSZz_RMpln
    INTEGER :: PostEgnd_RMpln
! Plane orientation and type used by slice and cross-section outputs.
    CHARACTER (LEN = 2) :: Planexyz_Global, PlaneType_Global
! Quadrature tables used by surface integrals.
    INTEGER :: n_glqln1d,n_glqtr2d
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: wg_glqln1d, &
    &   xg_glqln1d
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: wg_glqtr2d, &
    &   xg_glqtr2d, yg_glqtr2d
! Geometry, connectivity, local frames, and surface derivatives.

    CHARACTER (LEN = 1) :: MeshType

    INTEGER :: nmbrprtl

    CHARACTER (LEN = 4), ALLOCATABLE, DIMENSION (:) :: PrtlType

    CHARACTER (LEN = 7), ALLOCATABLE, DIMENSION (:) :: MeshGnrtn

    INTEGER, ALLOCATABLE, DIMENSION (:) :: NrmlInOut

    INTEGER, ALLOCATABLE, DIMENSION (:) :: MeshRead

    INTEGER, ALLOCATABLE, DIMENSION (:) :: Meshnlvl

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: MeshRlvlStp

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: &
    &   bvsa, cvsa, dvsa, bowl_a, bowl_b, &
    &   dfsp_a, dfsp_b, dfsp_c, dfsp_l, dfsp_m, dfsp_n

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: volume, surfarea

    INTEGER, ALLOCATABLE, DIMENSION (:) ::  nmbrnd, nmbrelmnt

    INTEGER, ALLOCATABLE, DIMENSION (:) ::  ndstaID, ndendID
    INTEGER, ALLOCATABLE, DIMENSION (:) ::  elstaID, elendID

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: &
    &   sizezoom, xloctn, yloctn, zloctn, &
    &   xmssctr, ymssctr, zmssctr, &
    &   xnrmref, ynrmref, znrmref

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: &
    &   anglecal_x, anglecal_y, anglecal_z

    INTEGER, ALLOCATABLE, DIMENSION (:) :: corelnkshell, particle_group

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: &
    &   xnd, ynd, znd, nnx, nny, nnz, t1x, t1y, t1z, t2x, t2y, t2z, &
    &   curvt1, curvt2, curvmn, curvt1th, curvt2th, curvmnth, &
    &   elmntarea, nnxelmnt, nnyelmnt, nnzelmnt

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: d_dt1, d_dt2, d2_dt1, d2_dt2

    INTEGER :: ttlnmbrnd, ttlnmbrelmnt, ttlddtnd, ttlgroup

    INTEGER :: mxnmbrnd, mxnmbrelmnt

    INTEGER :: mxnmbrndlnkelmnt, mxnmbrndlnkelmntlnr

    INTEGER :: mxnmbrndlnknd1st, mxnmbrndlnknd1stslf, &
    &   mxnmbrndlnknd2nd, mxnmbrndlnknd2ndslf

    INTEGER, ALLOCATABLE, DIMENSION (:,:) ::  elmntlnknd, ndlnkelmnt

    INTEGER, ALLOCATABLE, DIMENSION (:,:) ::  elmntlnkndlnr, ndlnkelmntlnr

    INTEGER, ALLOCATABLE, DIMENSION (:,:) ::  &
    &   ndlnknd1st, ndlnknd1stslf, ndlnknd2nd, ndlnknd2ndslf

    INTEGER :: ttlsrcfmm
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: &
    &   srcfmm_vec,srcfmm_nrm,srcfmm_wtnd
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)   :: srcfmm_wght

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)   :: rho_den

END MODULE
