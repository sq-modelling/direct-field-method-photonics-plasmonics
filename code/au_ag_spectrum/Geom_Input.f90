
! SPDX-FileCopyrightText: 2026 Qiang Sun
! SPDX-License-Identifier: BSD-3-Clause

! Read the record-oriented geometry description and initialise quadrature.
! Input_Geom.dat uses the same length unit for coordinates, radii, offsets,
! and mesh parameters; the electromagnetic input must use that unit as well.
!
MODULE Geom_Input

    USE omp_lib

    USE Pre_Constants
    USE Pre_csvformat

    USE Geom_GlobalData

    IMPLICIT NONE

    CONTAINS

    ! First-pass reader: initialise the fixed 6-point line and 16-point
    ! triangle rules, then read the case range, object count, and mesh type.

    SUBROUTINE GetGeomInput_Int

        INTEGER :: IOS, ithprtl, GLQi, tpi1dGQ, tpi2dGQ

        DOUBLE PRECISION :: tp, tp1, tp2, tp3, tp4, tpcff, tpzoff
        DOUBLE PRECISION :: tp_phiincfar, tpvcmwlkw, tpwlkwdecay

        n_glqln1d = 6
        ALLOCATE (wg_glqln1d(n_glqln1d))
        ALLOCATE (xg_glqln1d(n_glqln1d))
        CALL GauLegCoeff1D(n_glqln1d,xg_glqln1d,wg_glqln1d)
        n_glqtr2d = 16
        ALLOCATE (wg_glqtr2d(n_glqtr2d))
        ALLOCATE (xg_glqtr2d(n_glqtr2d))
        ALLOCATE (yg_glqtr2d(n_glqtr2d))
        CALL GauLegUniTriCff2D(n_glqtr2d,xg_glqtr2d,yg_glqtr2d,wg_glqtr2d)
        OPEN (71, FILE = "Input_Geom.dat", STATUS = "OLD", IOSTAT = IOS)
        IF (IOS /= 0) THEN
            PRINT*, "'Input_Geom.dat' does not exist! Please check!"
            STOP
        END IF

        READ (71, *) !Case start value, start ID, end ID, step number, and step size
        READ (71, *) csStartValue, csStartID, csEndID, csStepNumber, csStepSize
        READ (71, *)
        READ (71, *) !Number of the objects, and type of mesh
        READ (71, *) nmbrprtl, MeshType
        READ (71, *)

        CLOSE (71)

    END SUBROUTINE

    ! Full reader: allocate per-object arrays and load geometry, mesh-source,
    ! orientation, scale, translation, and interface-topology records for the
    ! current case. Generated and imported meshes share the resulting arrays.

    SUBROUTINE GetGeomInput

        INTEGER :: i, j, k, IOS, ithprtl, id_tp

        DOUBLE PRECISION :: tp, tp1, tp2, tp3, tp4, tpcff, tpzoff
        DOUBLE PRECISION :: tp_phiincfar, tpvcmwlkw, tpwlkwdecay

        OPEN (71, FILE = "Input_Geom.dat", STATUS = "OLD", IOSTAT = IOS)
        IF (IOS /= 0) THEN
            PRINT*, "'Input_Geom.dat' does not exist! Please check!"
            STOP
        END IF

        DO i = 1, 6
            READ (71, *)
        END DO

        ALLOCATE (PrtlType(nmbrprtl))
        ALLOCATE (MeshRead(nmbrprtl))
        ALLOCATE (NrmlInOut(nmbrprtl))
        ALLOCATE (MeshGnrtn(nmbrprtl))
        ALLOCATE (Meshnlvl(nmbrprtl))
        ALLOCATE (MeshRlvlStp(nmbrprtl))
        ALLOCATE (xnrmref(nmbrprtl))
        ALLOCATE (ynrmref(nmbrprtl))
        ALLOCATE (znrmref(nmbrprtl))

        ALLOCATE (sizezoom(nmbrprtl))
        ALLOCATE (xloctn(nmbrprtl))
        ALLOCATE (yloctn(nmbrprtl))
        ALLOCATE (zloctn(nmbrprtl))

        ALLOCATE (corelnkshell(nmbrprtl))

        ALLOCATE (bvsa(nmbrprtl))
        ALLOCATE (cvsa(nmbrprtl))
        ALLOCATE (dvsa(nmbrprtl))
        ALLOCATE (bowl_a(nmbrprtl))
        ALLOCATE (bowl_b(nmbrprtl))
        ALLOCATE (dfsp_a(nmbrprtl))
        ALLOCATE (dfsp_b(nmbrprtl))
        ALLOCATE (dfsp_c(nmbrprtl))
        ALLOCATE (dfsp_l(nmbrprtl))
        ALLOCATE (dfsp_m(nmbrprtl))
        ALLOCATE (dfsp_n(nmbrprtl))
        ALLOCATE (anglecal_x(nmbrprtl))
        ALLOCATE (anglecal_y(nmbrprtl))
        ALLOCATE (anglecal_z(nmbrprtl))
        ALLOCATE (surfarea(nmbrprtl))
        ALLOCATE (volume(nmbrprtl))

        ALLOCATE (nmbrnd(nmbrprtl))
        ALLOCATE (nmbrelmnt(nmbrprtl))
        ALLOCATE (ndstaID(nmbrprtl))
        ALLOCATE (ndendID(nmbrprtl))
        ALLOCATE (elstaID(nmbrprtl))
        ALLOCATE (elendID(nmbrprtl))

        ALLOCATE (rho_den(nmbrprtl))
        ALLOCATE (xmssctr(nmbrprtl))
        ALLOCATE (ymssctr(nmbrprtl))
        ALLOCATE (zmssctr(nmbrprtl))


        DO ithprtl = 1, nmbrprtl

            READ (71, *)

            READ (71, *) !Type of particle:
            !   sphere Sphr; long Rod LgRd; Dumb-Bell DbBl; Prolong Spheroid PrSp
            !   link between core and shell particle id
            !   Type of mesh to read
            !   Normal direction
            READ (71, *) PrtlType(ithprtl), corelnkshell(ithprtl), &
            &   MeshRead(ithprtl), NrmlInOut(ithprtl), &
            &   xnrmref(ithprtl),ynrmref(ithprtl),znrmref(ithprtl)

            READ (71, *)
            READ (71, *) !Type of mesh: "Icshdrl" or "TwoCrcl";
            !   Level of mesh on hemisphere,
            !   uneven mesh step scale (0.01~1.99, 1.0->even)
            READ (71, *) MeshGnrtn(ithprtl), Meshnlvl(ithprtl), MeshRlvlStp(ithprtl)
            READ (71, *)

            READ (71, *) !if need to zoom and move the particle
            READ (71, *) sizezoom(ithprtl), &
            &            xloctn(ithprtl),yloctn(ithprtl),zloctn(ithprtl)
            READ (71, *)

            READ (71, *) !The shape function
            !   b for prolong spheroid, long rod;
            !   d for dumb-bell;
            !   bowl_a and bowl_b;
            !   other shape parameters.
            READ (71, *) tp, bowl_a(ithprtl), bowl_b(ithprtl), &
            &            dfsp_a(ithprtl), dfsp_b(ithprtl), dfsp_c(ithprtl), &
            &            dfsp_l(ithprtl), dfsp_m(ithprtl), dfsp_n(ithprtl)
            READ (71, *)

            bvsa(ithprtl) = 0.0d0
            cvsa(ithprtl) = 0.0d0
            dvsa(ithprtl) = 0.0d0
            IF (PrtlType(ithprtl) == "PrSp") bvsa(ithprtl) = tp     !prolate ellipsoid
            IF (PrtlType(ithprtl) == "LgRd") bvsa(ithprtl) = tp     !long rod
            IF (PrtlType(ithprtl) == "Disk") bvsa(ithprtl) = tp     !Disk
            IF (PrtlType(ithprtl) == "DbBl") dvsa(ithprtl) = tp     !Dumbbell-sphaped
            IF (PrtlType(ithprtl) == "ObSp") bvsa(ithprtl) = tp     !Oblate ellipsoid
            IF (PrtlType(ithprtl) == "Coin") dvsa(ithprtl) = tp     !Coin
            IF (PrtlType(ithprtl) == "JoJo") bvsa(ithprtl) = tp     !JoJo
            IF (PrtlType(ithprtl) == "Pnty") bvsa(ithprtl) = tp     !Pointy
            IF (PrtlType(ithprtl) == "ChPr") bvsa(ithprtl) = tp     !Ch prolate ellipsoid
            IF (PrtlType(ithprtl) == "ChOb") bvsa(ithprtl) = tp     !Ch oblate ellipsoid

            READ (71, *) !angle using the row (x), pitch (y) and yaw (z) idea
            READ (71, *) anglecal_x(ithprtl),anglecal_y(ithprtl),anglecal_z(ithprtl)
            READ (71, *)

            anglecal_x(ithprtl) = anglecal_x(ithprtl)/180.0d0*pai
            anglecal_y(ithprtl) = anglecal_y(ithprtl)/180.0d0*pai
            anglecal_z(ithprtl) = anglecal_z(ithprtl)/180.0d0*pai

            READ (71, *) !particle density
            READ (71, *) rho_den(ithprtl), &
            &            xmssctr(ithprtl), ymssctr(ithprtl), zmssctr(ithprtl)
            READ (71, *)

        END DO

        CLOSE (71)

    END SUBROUTINE

END MODULE
