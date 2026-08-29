
! SPDX-FileCopyrightText: 2026 Qiang Sun
! SPDX-License-Identifier: BSD-3-Clause

MODULE Geom_MeshSphereCircle

    ! Unit-sphere mesh generators used by Geom_Mesh.  Each generator writes
    ! Prtl_Orgnl.vrt (node_id,x,y,z), Prtl_Orgnl.cel (element connectivity),
    ! and Prtl_Orgnl.inp (node/element counts and the origin as reference point).
    !
    ! The icosahedral topology is derived from sphere-mesh code by Evert
    ! Klaseboer.  The Simple variant interpolates polar/azimuthal coordinates;
    ! the Paper variant uses minor great-circle interpolation.  TwoCrcl instead
    ! joins mirrored hemispheres built from concentric projected circles and is
    ! the only family here that applies the radial grading dmRlvlStpScl.
    !
    ! Simple, Paper, and TwoCrcl use an internal grid with 2*dmnlvlprtl
    ! intervals.  L/C elements use every interval; Q6 elements span two and
    ! therefore retain the intervening nodes as midside nodes.  Connectivity is
    ! [corner1,corner2,corner3] for L/C and [corner1,corner2,corner3,mid12,
    ! mid23,mid31] for Q6.  Geom_Mesh subsequently orients every element from
    ! (corner1-corner3)x(corner2-corner3): NrmlInOut=+1 points toward the
    ! reference origin and NrmlInOut=-1 points away from it.

    USE omp_lib

    USE Pre_Constants
    USE Pre_csvformat

    IMPLICIT NONE

    CONTAINS

    ! Generate the EK-derived icosahedral topology using angle interpolation.
    ! dmRlvlStpScl is unused; it is retained for the common generator interface.
    SUBROUTINE GetMeshNdElmntIcshdrlSimple(dmnlvlprtl, dmRlvlStpScl, dmMeshType)

        INTEGER, INTENT (IN) :: dmnlvlprtl

        DOUBLE PRECISION, INTENT (IN) :: dmRlvlStpScl

        CHARACTER (LEN=1), INTENT (IN) :: dmMeshType

        DOUBLE PRECISION :: SphereMesh_pi = DATAN(1.0d0)*4.0d0

        DOUBLE PRECISION, DIMENSION (12) :: xcshd, ycshd, zcshd

        DOUBLE PRECISION :: tp, tp1, tpAphi, tpBphi, tpphi, tpAtht, tpBtht, tptht, tprho, &
        &   cshdxy, cshdz, tpAx, tpAy, tpAz, tpBx, tpBy, tpBz, &
        &   tpi_n, tpxcdnt, tpycdnt, tpzcdnt

        INTEGER :: i, j, k, NdCnt, nlvl, ElmtCnt, lvlmk, lvlnd, lvli, pstni, pstnj

        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: xsphr, ysphr ,zsphr

        INTEGER, ALLOCATABLE, DIMENSION (:, :, :) :: NdIDRgbsCll

        INTEGER, ALLOCATABLE, DIMENSION (:, :) :: NdIDEg

        INTEGER, ALLOCATABLE, DIMENSION (:) :: ElmtLnk1, ElmtLnk2, ElmtLnk3

        ! Rotated unit-icosahedron vertices.

        tp = 2.0d0*DATAN( 2.0d0/( 1.0d0+DSQRT(5.0d0) ) )   !golden ratio with rotation
        cshdz  = DCOS(tp)
        cshdxy = DSIN(tp)

        tp = 0.4d0*SphereMesh_pi    !(2*pi/5)

        xcshd(1) = 0.0d0
        ycshd(1) = 0.0d0
        zcshd(1) = 1.0d0

        DO i = 2, 6

            xcshd(i) = cshdxy*DCOS( tp * DBLE(i-2)  )
            ycshd(i) = cshdxy*DSIN( tp * DBLE(i-2)  )
            zcshd(i) = cshdz

        END DO

        xcshd(7) =  0.0d0
        ycshd(7) =  0.0d0
        zcshd(7) = -1.0d0

        DO i = 8, 12

            xcshd(i) = cshdxy*DCOS( SphereMesh_pi + tp * DBLE(i-8)  )
            ycshd(i) = cshdxy*DSIN( SphereMesh_pi + tp * DBLE(i-8)  )
            zcshd(i) =-cshdz

        END DO

        ! Subdivide the icosahedron and place all nodes on the unit sphere.

        nlvl = 2*dmnlvlprtl

        ALLOCATE (xsphr(10*nlvl*nlvl+2))
        ALLOCATE (ysphr(10*nlvl*nlvl+2))
        ALLOCATE (zsphr(10*nlvl*nlvl+2))

        ALLOCATE (ElmtLnk1(10*2*nlvl*nlvl))
        ALLOCATE (ElmtLnk2(10*2*nlvl*nlvl))
        ALLOCATE (ElmtLnk3(10*2*nlvl*nlvl))

        ALLOCATE (NdIDRgbsCll(10, nlvl+1, nlvl+1))

        ALLOCATE (NdIDEg(20, nlvl+1))

        DO i = 1, 12
            xsphr(i) = xcshd(i)
            ysphr(i) = ycshd(i)
            zsphr(i) = zcshd(i)
        END DO

        NdCnt = 12

        ! Record the ordered nodes on the 20 icosahedron edges.
        DO k = 1, 20

            IF (k == 1) THEN ! 1_2
                NdIDEg(k,1) = 1
                NdIDEg(k,nlvl+1) = 2
            END IF
            IF (k == 2) THEN ! 1_3
                NdIDEg(k,1) = 1
                NdIDEg(k,nlvl+1) = 3
            END IF
            IF (k == 3) THEN ! 1_4
                NdIDEg(k,1) = 1
                NdIDEg(k,nlvl+1) = 4
            END IF
            IF (k == 4) THEN ! 1_5
                NdIDEg(k,1) = 1
                NdIDEg(k,nlvl+1) = 5
            END IF
            IF (k == 5) THEN ! 1_6
                NdIDEg(k,1) = 1
                NdIDEg(k,nlvl+1) = 6
            END IF

            IF (k == 6) THEN ! 2_11
                NdIDEg(k,1) = 2
                NdIDEg(k,nlvl+1) = 11
            END IF
            IF (k == 7) THEN ! 2_10
                NdIDEg(k,1) = 2
                NdIDEg(k,nlvl+1) = 10
            END IF
            IF (k == 8) THEN ! 3_11
                NdIDEg(k,1) = 3
                NdIDEg(k,nlvl+1) = 11
            END IF
            IF (k == 9) THEN ! 3_12
                NdIDEg(k,1) = 3
                NdIDEg(k,nlvl+1) = 12
            END IF
            IF (k == 10) THEN ! 4_12
                NdIDEg(k,1) = 4
                NdIDEg(k,nlvl+1) = 12
            END IF

            IF (k == 11) THEN ! 4_8
                NdIDEg(k,1) = 4
                NdIDEg(k,nlvl+1) = 8
            END IF
            IF (k == 12) THEN ! 5_8
                NdIDEg(k,1) = 5
                NdIDEg(k,nlvl+1) = 8
            END IF
            IF (k == 13) THEN ! 5_9
                NdIDEg(k,1) = 5
                NdIDEg(k,nlvl+1) = 9
            END IF
            IF (k == 14) THEN ! 6_9
                NdIDEg(k,1) = 6
                NdIDEg(k,nlvl+1) = 9
            END IF
            IF (k == 15) THEN ! 6_10
                NdIDEg(k,1) = 6
                NdIDEg(k,nlvl+1) = 10
            END IF

            IF (k == 16) THEN ! 8_7
                NdIDEg(k,1) = 8
                NdIDEg(k,nlvl+1) = 7
            END IF
            IF (k == 17) THEN ! 9_7
                NdIDEg(k,1) = 9
                NdIDEg(k,nlvl+1) = 7
            END IF
            IF (k == 18) THEN ! 10_7
                NdIDEg(k,1) = 10
                NdIDEg(k,nlvl+1) = 7
            END IF
            IF (k == 19) THEN ! 11_7
                NdIDEg(k,1) = 11
                NdIDEg(k,nlvl+1) = 7
            END IF
            IF (k == 20) THEN ! 12_7
                NdIDEg(k,1) = 12
                NdIDEg(k,nlvl+1) = 7
            END IF

            tpAx = xsphr(NdIDEg(k,1))
            tpAy = ysphr(NdIDEg(k,1))
            tpAz = zsphr(NdIDEg(k,1))
            tpBx = xsphr(NdIDEg(k,nlvl+1))
            tpBy = ysphr(NdIDEg(k,nlvl+1))
            tpBz = zsphr(NdIDEg(k,nlvl+1))

            IF (DABS(tpAy) < 1.0E-7) THEN
                IF (tpAx > 0.0d0) tpAphi = 0.0d0
                IF (tpAx < 0.0d0) tpAphi = SphereMesh_pi
            ELSE IF (tpAy > 0.0d0) THEN
                IF (DABS(tpAx) < 1.0E-7) THEN
                    tpAphi = 0.5d0*SphereMesh_pi
                ELSE IF (tpAx > 0.0d0) THEN
                    tpAphi = DATAN(DABS(tpAy/tpAx))
                ELSE IF (tpAx < 0.0d0) THEN
                    tpAphi = SphereMesh_pi-DATAN(DABS(tpAy/tpAx))
                END IF
            ELSE IF (tpAy < 0.0d0) THEN
                IF (DABS(tpAx) < 1.0E-7) THEN
                    tpAphi = 1.5d0*SphereMesh_pi
                ELSE IF (tpAx > 0.0d0) THEN
                    tpAphi = 2.0d0*SphereMesh_pi-DATAN(DABS(tpAy/tpAx))
                ELSE IF (tpAx < 0.0d0) THEN
                    tpAphi = SphereMesh_pi+DATAN(DABS(tpAy/tpAx))
                END IF
            END IF

            IF (DABS(tpBy) < 1.0E-7) THEN
                IF (tpBx > 0.0d0) tpBphi = 0.0d0
                IF (tpBx < 0.0d0) tpBphi = SphereMesh_pi
            ELSE IF (tpBy > 0.0d0) THEN
                IF (DABS(tpBx) < 1.0E-7) THEN
                    tpBphi = 0.5d0*SphereMesh_pi
                ELSE IF (tpBx > 0.0d0) THEN
                    tpBphi = DATAN(DABS(tpBy/tpBx))
                ELSE IF (tpBx < 0.0d0) THEN
                    tpBphi = SphereMesh_pi-DATAN(DABS(tpBy/tpBx))
                END IF
            ELSE IF (tpBy < 0.0d0) THEN
                IF (DABS(tpBx) < 1.0E-7) THEN
                    tpBphi = 1.5d0*SphereMesh_pi
                ELSE IF (tpBx > 0.0d0) THEN
                    tpBphi = 2.0d0*SphereMesh_pi-DATAN(DABS(tpBy/tpBx))
                ELSE IF (tpBx < 0.0d0) THEN
                    tpBphi = SphereMesh_pi+DATAN(DABS(tpBy/tpBx))
                END IF
            END IF

            IF (NdIDEg(k,1) == 1) THEN
                tpAphi = tpBphi
            END IF
            IF (NdIDEg(k,nlvl+1) == 7) THEN
                tpBphi = tpAphi
            END IF

            tpphi = (tpBphi-tpAphi)
            IF (tpphi > 0.40001*SphereMesh_pi) THEN
                tpphi = tpphi - 2.0d0*SphereMesh_pi
            END IF
            IF (tpphi <-0.40001*SphereMesh_pi) THEN
                tpphi = tpphi + 2.0d0*SphereMesh_pi
            END IF

            tpAtht = DACOS(tpAz)
            tpBtht = DACOS(tpBz)
            tptht = (tpBtht-tpAtht)

            DO i = 2, nlvl

                NdCnt = NdCnt + 1
                NdIDEg(k,i) = NdCnt

                tpi_n = DBLE(i-1)/DBLE(nlvl)

                tprho = DSIN( tpAtht + tptht*tpi_n )
                xsphr(NdIDEg(k,i)) = tprho * DCOS( tpAphi + tpphi*tpi_n )
                ysphr(NdIDEg(k,i)) = tprho * DSIN( tpAphi + tpphi*tpi_n )
                zsphr(NdIDEg(k,i)) = DCOS( tpAtht + tptht*tpi_n )

            END DO

        END DO

        DO k = 1, 10

            IF (k == 1) THEN !1_2, 1_3, 2_11, 3_11
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(1,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(2,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(6,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(8,i)
                END DO
            END IF
            IF (k == 2) THEN !1_3, 1_4, 3_12, 4_12
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(2,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(3,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(9,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(10,i)
                END DO
            END IF
            IF (k == 3) THEN !1_4, 1_5, 4_8, 5_8
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(3,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(4,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(11,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(12,i)
                END DO
            END IF
            IF (k == 4) THEN !1_5, 1_6, 5_9, 6_9
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(4,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(5,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(13,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(14,i)
                END DO
            END IF
            IF (k == 5) THEN !1_6, 1_2, 6_10, 2_10
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(5,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(1,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(15,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(7,i)
                END DO
            END IF
            IF (k == 6) THEN !2_11, 2_10, 11_7, 10_7
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(6,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(7,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(19,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(18,i)
                END DO
            END IF
            IF (k == 7) THEN !3_11, 3_12, 11_7, 12_7
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(8,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(9,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(19,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(20,i)
                END DO
            END IF
            IF (k == 8) THEN !4_12, 4_8, 12_7, 8_7
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(10,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(11,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(20,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(16,i)
                END DO
            END IF
            IF (k == 9) THEN !5_8, 5_9, 8_7, 9_7
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(12,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(13,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(16,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(17,i)
                END DO
            END IF
            IF (k == 10) THEN !6_9, 6_10, 9_7, 10_7
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(14,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(15,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(17,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(18,i)
                END DO
            END IF

            DO lvlmk = 2, nlvl

                tpAx = xsphr(NdIDRgbsCll(k,lvlmk+1,1))
                tpAy = ysphr(NdIDRgbsCll(k,lvlmk+1,1))
                tpAz = zsphr(NdIDRgbsCll(k,lvlmk+1,1))
                tpBx = xsphr(NdIDRgbsCll(k,1,lvlmk+1))
                tpBy = ysphr(NdIDRgbsCll(k,1,lvlmk+1))
                tpBz = zsphr(NdIDRgbsCll(k,1,lvlmk+1))

                IF (DABS(tpAy) < 1.0E-7) THEN
                    IF (tpAx > 0.0d0) tpAphi = 0.0d0
                    IF (tpAx < 0.0d0) tpAphi = SphereMesh_pi
                ELSE IF (tpAy > 0.0d0) THEN
                    IF (DABS(tpAx) < 1.0E-7) THEN
                        tpAphi = 0.5d0*SphereMesh_pi
                    ELSE IF (tpAx > 0.0d0) THEN
                        tpAphi = DATAN(DABS(tpAy/tpAx))
                    ELSE IF (tpAx < 0.0d0) THEN
                        tpAphi = SphereMesh_pi-DATAN(DABS(tpAy/tpAx))
                    END IF
                ELSE IF (tpAy < 0.0d0) THEN
                    IF (DABS(tpAx) < 1.0E-7) THEN
                        tpAphi = 1.5d0*SphereMesh_pi
                    ELSE IF (tpAx > 0.0d0) THEN
                        tpAphi = 2.0d0*SphereMesh_pi-DATAN(DABS(tpAy/tpAx))
                    ELSE IF (tpAx < 0.0d0) THEN
                        tpAphi = SphereMesh_pi+DATAN(DABS(tpAy/tpAx))
                    END IF
                END IF

                IF (DABS(tpBy) < 1.0E-7) THEN
                    IF (tpBx > 0.0d0) tpBphi = 0.0d0
                    IF (tpBx < 0.0d0) tpBphi = SphereMesh_pi
                ELSE IF (tpBy > 0.0d0) THEN
                    IF (DABS(tpBx) < 1.0E-7) THEN
                        tpBphi = 0.5d0*SphereMesh_pi
                    ELSE IF (tpBx > 0.0d0) THEN
                        tpBphi = DATAN(DABS(tpBy/tpBx))
                    ELSE IF (tpBx < 0.0d0) THEN
                        tpBphi = SphereMesh_pi-DATAN(DABS(tpBy/tpBx))
                    END IF
                ELSE IF (tpBy < 0.0d0) THEN
                    IF (DABS(tpBx) < 1.0E-7) THEN
                        tpBphi = 1.5d0*SphereMesh_pi
                    ELSE IF (tpBx > 0.0d0) THEN
                        tpBphi = 2.0d0*SphereMesh_pi-DATAN(DABS(tpBy/tpBx))
                    ELSE IF (tpBx < 0.0d0) THEN
                        tpBphi = SphereMesh_pi+DATAN(DABS(tpBy/tpBx))
                    END IF
                END IF

                tpphi = (tpBphi-tpAphi)
                IF (tpphi > 0.40001*SphereMesh_pi) THEN
                    tpphi = tpphi - 2.0d0*SphereMesh_pi
                END IF
                IF (tpphi <-0.40001*SphereMesh_pi) THEN
                    tpphi = tpphi + 2.0d0*SphereMesh_pi
                END IF
                tprho = DSQRT(tpAx**2+tpAy**2)

                lvli = 0

                DO lvlnd = 1, lvlmk-1

                    pstni = lvlmk-(lvlnd-1)
                    pstnj = (lvlmk+2)-pstni

                    NdCnt = NdCnt + 1
                    NdIDRgbsCll(k,pstni,pstnj) = NdCnt

                    lvli = lvli + 1

                    tpi_n = DBLE(lvli)/DBLE(lvlmk)

                    tpxcdnt = tprho*DCOS(tpAphi + tpphi*tpi_n)
                    tpycdnt = tprho*DSIN(tpAphi + tpphi*tpi_n)
                    xsphr(NdIDRgbsCll(k,pstni,pstnj)) = tpxcdnt
                    ysphr(NdIDRgbsCll(k,pstni,pstnj)) = tpycdnt
                    zsphr(NdIDRgbsCll(k,pstni,pstnj)) = tpAz


                END DO

            END DO
            DO lvlmk = nlvl-1, 2, -1

                tpAx = xsphr(NdIDRgbsCll(k,nlvl+1,nlvl+1-lvlmk))
                tpAy = ysphr(NdIDRgbsCll(k,nlvl+1,nlvl+1-lvlmk))
                tpAz = zsphr(NdIDRgbsCll(k,nlvl+1,nlvl+1-lvlmk))
                tpBx = xsphr(NdIDRgbsCll(k,nlvl+1-lvlmk,nlvl+1))
                tpBy = ysphr(NdIDRgbsCll(k,nlvl+1-lvlmk,nlvl+1))
                tpBz = zsphr(NdIDRgbsCll(k,nlvl+1-lvlmk,nlvl+1))

                IF (DABS(tpAy) < 1.0E-7) THEN
                    IF (tpAx > 0.0d0) tpAphi = 0.0d0
                    IF (tpAx < 0.0d0) tpAphi = SphereMesh_pi
                ELSE IF (tpAy > 0.0d0) THEN
                    IF (DABS(tpAx) < 1.0E-7) THEN
                        tpAphi = 0.5d0*SphereMesh_pi
                    ELSE IF (tpAx > 0.0d0) THEN
                        tpAphi = DATAN(DABS(tpAy/tpAx))
                    ELSE IF (tpAx < 0.0d0) THEN
                        tpAphi = SphereMesh_pi-DATAN(DABS(tpAy/tpAx))
                    END IF
                ELSE IF (tpAy < 0.0d0) THEN
                    IF (DABS(tpAx) < 1.0E-7) THEN
                        tpAphi = 1.5d0*SphereMesh_pi
                    ELSE IF (tpAx > 0.0d0) THEN
                        tpAphi = 2.0d0*SphereMesh_pi-DATAN(DABS(tpAy/tpAx))
                    ELSE IF (tpAx < 0.0d0) THEN
                        tpAphi = SphereMesh_pi+DATAN(DABS(tpAy/tpAx))
                    END IF
                END IF

                IF (DABS(tpBy) < 1.0E-7) THEN
                    IF (tpBx > 0.0d0) tpBphi = 0.0d0
                    IF (tpBx < 0.0d0) tpBphi = SphereMesh_pi
                ELSE IF (tpBy > 0.0d0) THEN
                    IF (DABS(tpBx) < 1.0E-7) THEN
                        tpBphi = 0.5d0*SphereMesh_pi
                    ELSE IF (tpBx > 0.0d0) THEN
                        tpBphi = DATAN(DABS(tpBy/tpBx))
                    ELSE IF (tpBx < 0.0d0) THEN
                        tpBphi = SphereMesh_pi-DATAN(DABS(tpBy/tpBx))
                    END IF
                ELSE IF (tpBy < 0.0d0) THEN
                    IF (DABS(tpBx) < 1.0E-7) THEN
                        tpBphi = 1.5d0*SphereMesh_pi
                    ELSE IF (tpBx > 0.0d0) THEN
                        tpBphi = 2.0d0*SphereMesh_pi-DATAN(DABS(tpBy/tpBx))
                    ELSE IF (tpBx < 0.0d0) THEN
                        tpBphi = SphereMesh_pi+DATAN(DABS(tpBy/tpBx))
                    END IF
                END IF

                tpphi = (tpBphi-tpAphi)
                IF (tpphi > 0.40001*SphereMesh_pi) THEN
                    tpphi = tpphi - 2.0d0*SphereMesh_pi
                END IF
                IF (tpphi <-0.40001*SphereMesh_pi) THEN
                    tpphi = tpphi + 2.0d0*SphereMesh_pi
                END IF
                tprho = DSQRT(tpAx**2+tpAy**2)

                lvli = 0

                DO lvlnd = 1, lvlmk-1

                    pstnj = (nlvl+1)-(lvlmk-lvlnd)
                    pstni = 2*(nlvl+1)-lvlmk-pstnj

                    NdCnt = NdCnt + 1
                    NdIDRgbsCll(k,pstni,pstnj) = NdCnt

                    lvli = lvli + 1

                    tpi_n = DBLE(lvli)/DBLE(lvlmk)

                    tpxcdnt = tprho*DCOS(tpAphi + tpphi*tpi_n)
                    tpycdnt = tprho*DSIN(tpAphi + tpphi*tpi_n)
                    xsphr(NdIDRgbsCll(k,pstni,pstnj)) = tpxcdnt
                    ysphr(NdIDRgbsCll(k,pstni,pstnj)) = tpycdnt
                    zsphr(NdIDRgbsCll(k,pstni,pstnj)) = tpAz

                END DO

            END DO

        END DO

        OPEN (111, FILE = "Prtl_Orgnl.vrt", STATUS = "REPLACE")

        DO i = 1, 10*nlvl*nlvl+2
            WRITE (111, *) i, xsphr(i), ysphr(i), zsphr(i)
        END DO

        CLOSE (111)

        OPEN (121, FILE = "Prtl_Orgnl.cel", STATUS = "REPLACE")

        IF (dmMeshType == "L" .OR. dmMeshType == "C") THEN

            ElmtCnt = 0
            DO k =1, 10
                DO i = 1, nlvl
                    DO j = 1, nlvl
                        ElmtCnt = ElmtCnt + 1
                        WRITE (121, *) ElmtCnt, NdIDRgbsCll(k,i,j), &
                                     & NdIDRgbsCll(k,i,j+1), NdIDRgbsCll(k,i+1,j)
                        ElmtCnt = ElmtCnt + 1
                        WRITE (121, *) ElmtCnt, NdIDRgbsCll(k,i,j+1), &
                                    &  NdIDRgbsCll(k,i+1,j+1), NdIDRgbsCll(k,i+1,j)
                    END DO
                END DO
            END DO

        END IF

        IF (dmMeshType == "Q") THEN

            ElmtCnt = 0
            DO k =1, 10
                DO i = 1, nlvl+1-2,2
                    DO j = 1, nlvl+1-2,2
                        ElmtCnt = ElmtCnt + 1
                        WRITE (121, *) ElmtCnt, &
                                    &  NdIDRgbsCll(k,i,j), &
                                    &  NdIDRgbsCll(k,i,j+2), &
                                    &  NdIDRgbsCll(k,i+2,j), &
                                    &  NdIDRgbsCll(k,i,j+1), &
                                    &  NdIDRgbsCll(k,i+1,j+1), &
                                    &  NdIDRgbsCll(k,i+1,j)
                        ElmtCnt = ElmtCnt + 1
                        WRITE (121, *) ElmtCnt, &
                                    &  NdIDRgbsCll(k,i,j+2), &
                                    &  NdIDRgbsCll(k,i+2,j+2), &
                                    &  NdIDRgbsCll(k,i+2,j), &
                                    &  NdIDRgbsCll(k,i+1,j+2), &
                                    &  NdIDRgbsCll(k,i+2,j+1), &
                                    &  NdIDRgbsCll(k,i+1,j+1)
                    END DO
                END DO
            END DO

        END IF

        CLOSE (121)

        OPEN (101, FILE = "Prtl_Orgnl.inp", STATUS = "REPLACE")

        WRITE (101, *)
        WRITE (101, *)
        WRITE (101, *) 10*nlvl*nlvl+2
        WRITE (101, *) ElmtCnt
        WRITE (101, *) 0.0d0, 0.0d0, 0.0d0

        CLOSE (101)

        DEALLOCATE (xsphr)
        DEALLOCATE (ysphr)
        DEALLOCATE (zsphr)
        DEALLOCATE (NdIDRgbsCll)
        DEALLOCATE (NdIDEg)
        DEALLOCATE (ElmtLnk1)
        DEALLOCATE (ElmtLnk2)
        DEALLOCATE (ElmtLnk3)

    END SUBROUTINE


    ! Generate the same icosahedral topology as the Simple variant, placing
    ! edge and patch nodes along minor great-circle arcs via NdCodntssphr.
    ! dmRlvlStpScl is unused; it is retained for the common generator interface.
    SUBROUTINE GetMeshNdElmntIcshdrlPaper(dmnlvlprtl, dmRlvlStpScl, dmMeshType)

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: dmnlvlprtl

        DOUBLE PRECISION, INTENT (IN) :: dmRlvlStpScl

        CHARACTER (LEN=1), INTENT (IN) :: dmMeshType

        DOUBLE PRECISION :: SphereMesh_pi = DATAN(1.0d0)*4.0d0

        DOUBLE PRECISION, DIMENSION (12) :: xcshd, ycshd, zcshd

        DOUBLE PRECISION :: tp, tp1, cshdxy, cshdz, &
        &                   tpAx, tpAy, tpAz, tpBx, tpBy, tpBz, &
        &                   tpi_n, tpxcdnt, tpycdnt, tpzcdnt

        INTEGER :: i, j, k, NdCnt, nlvl, ElmtCnt, lvlmk, lvlnd, lvli, pstni, pstnj

        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: xsphr, ysphr ,zsphr

        INTEGER, ALLOCATABLE, DIMENSION (:, :, :) :: NdIDRgbsCll

        INTEGER, ALLOCATABLE, DIMENSION (:, :) :: NdIDEg

        INTEGER, ALLOCATABLE, DIMENSION (:) :: ElmtLnk1, ElmtLnk2, ElmtLnk3

        ! Rotated unit-icosahedron vertices.

        tp = 2.0d0*DATAN( 2.0d0/( 1.0d0+DSQRT(5.0d0) ) )   !golden ratio with rotation
        cshdz  = DCOS(tp)
        cshdxy = DSIN(tp)

        tp = 0.4d0*SphereMesh_pi    !(2*pi/5)

        xcshd(1) = 0.0d0
        ycshd(1) = 0.0d0
        zcshd(1) = 1.0d0

        DO i = 2, 6

            xcshd(i) = cshdxy*DCOS( tp * DBLE(i-2)  )
            ycshd(i) = cshdxy*DSIN( tp * DBLE(i-2)  )
            zcshd(i) = cshdz

        END DO

        xcshd(7) =  0.0d0
        ycshd(7) =  0.0d0
        zcshd(7) = -1.0d0

        DO i = 8, 12

            xcshd(i) = cshdxy*DCOS( SphereMesh_pi + tp * DBLE(i-8)  )
            ycshd(i) = cshdxy*DSIN( SphereMesh_pi + tp * DBLE(i-8)  )
            zcshd(i) =-cshdz

        END DO

        ! Subdivide the icosahedron and place all nodes on great-circle arcs.

        nlvl = 2*dmnlvlprtl

        ALLOCATE (xsphr(10*nlvl*nlvl+2))
        ALLOCATE (ysphr(10*nlvl*nlvl+2))
        ALLOCATE (zsphr(10*nlvl*nlvl+2))

        ALLOCATE (ElmtLnk1(10*2*nlvl*nlvl))
        ALLOCATE (ElmtLnk2(10*2*nlvl*nlvl))
        ALLOCATE (ElmtLnk3(10*2*nlvl*nlvl))

        ALLOCATE (NdIDRgbsCll(10, nlvl+1, nlvl+1))

        ALLOCATE (NdIDEg(20, nlvl+1))

        DO i = 1, 12
            xsphr(i) = xcshd(i)
            ysphr(i) = ycshd(i)
            zsphr(i) = zcshd(i)
        END DO

        NdCnt = 12

        ! Record the ordered nodes on the 20 icosahedron edges.
        DO k = 1, 20

            IF (k == 1) THEN ! 1_2
                NdIDEg(k,1) = 1
                NdIDEg(k,nlvl+1) = 2
            END IF
            IF (k == 2) THEN ! 1_3
                NdIDEg(k,1) = 1
                NdIDEg(k,nlvl+1) = 3
            END IF
            IF (k == 3) THEN ! 1_4
                NdIDEg(k,1) = 1
                NdIDEg(k,nlvl+1) = 4
            END IF
            IF (k == 4) THEN ! 1_5
                NdIDEg(k,1) = 1
                NdIDEg(k,nlvl+1) = 5
            END IF
            IF (k == 5) THEN ! 1_6
                NdIDEg(k,1) = 1
                NdIDEg(k,nlvl+1) = 6
            END IF

            IF (k == 6) THEN ! 2_11
                NdIDEg(k,1) = 2
                NdIDEg(k,nlvl+1) = 11
            END IF
            IF (k == 7) THEN ! 2_10
                NdIDEg(k,1) = 2
                NdIDEg(k,nlvl+1) = 10
            END IF
            IF (k == 8) THEN ! 3_11
                NdIDEg(k,1) = 3
                NdIDEg(k,nlvl+1) = 11
            END IF
            IF (k == 9) THEN ! 3_12
                NdIDEg(k,1) = 3
                NdIDEg(k,nlvl+1) = 12
            END IF
            IF (k == 10) THEN ! 4_12
                NdIDEg(k,1) = 4
                NdIDEg(k,nlvl+1) = 12
            END IF

            IF (k == 11) THEN ! 4_8
                NdIDEg(k,1) = 4
                NdIDEg(k,nlvl+1) = 8
            END IF
            IF (k == 12) THEN ! 5_8
                NdIDEg(k,1) = 5
                NdIDEg(k,nlvl+1) = 8
            END IF
            IF (k == 13) THEN ! 5_9
                NdIDEg(k,1) = 5
                NdIDEg(k,nlvl+1) = 9
            END IF
            IF (k == 14) THEN ! 6_9
                NdIDEg(k,1) = 6
                NdIDEg(k,nlvl+1) = 9
            END IF
            IF (k == 15) THEN ! 6_10
                NdIDEg(k,1) = 6
                NdIDEg(k,nlvl+1) = 10
            END IF

            IF (k == 16) THEN ! 8_7
                NdIDEg(k,1) = 8
                NdIDEg(k,nlvl+1) = 7
            END IF
            IF (k == 17) THEN ! 9_7
                NdIDEg(k,1) = 9
                NdIDEg(k,nlvl+1) = 7
            END IF
            IF (k == 18) THEN ! 10_7
                NdIDEg(k,1) = 10
                NdIDEg(k,nlvl+1) = 7
            END IF
            IF (k == 19) THEN ! 11_7
                NdIDEg(k,1) = 11
                NdIDEg(k,nlvl+1) = 7
            END IF
            IF (k == 20) THEN ! 12_7
                NdIDEg(k,1) = 12
                NdIDEg(k,nlvl+1) = 7
            END IF

            tpAx = xsphr(NdIDEg(k,1))
            tpAy = ysphr(NdIDEg(k,1))
            tpAz = zsphr(NdIDEg(k,1))
            tpBx = xsphr(NdIDEg(k,nlvl+1))
            tpBy = ysphr(NdIDEg(k,nlvl+1))
            tpBz = zsphr(NdIDEg(k,nlvl+1))

            DO i = 2, nlvl

                NdCnt = NdCnt + 1
                NdIDEg(k,i) = NdCnt


                tpi_n = DBLE(i-1)/DBLE(nlvl)
                CALL NdCodntssphr(tpAx, tpAy, tpAz, tpBx, tpBy, tpBz, &
                &                 tpi_n, tpxcdnt, tpycdnt, tpzcdnt)
                xsphr(NdIDEg(k,i)) = tpxcdnt
                ysphr(NdIDEg(k,i)) = tpycdnt
                zsphr(NdIDEg(k,i)) = tpzcdnt

            END DO

        END DO

        DO k = 1, 10

            IF (k == 1) THEN !1_2, 1_3, 2_11, 3_11
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(1,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(2,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(6,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(8,i)
                END DO
            END IF
            IF (k == 2) THEN !1_3, 1_4, 3_12, 4_12
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(2,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(3,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(9,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(10,i)
                END DO
            END IF
            IF (k == 3) THEN !1_4, 1_5, 4_8, 5_8
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(3,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(4,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(11,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(12,i)
                END DO
            END IF
            IF (k == 4) THEN !1_5, 1_6, 5_9, 6_9
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(4,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(5,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(13,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(14,i)
                END DO
            END IF
            IF (k == 5) THEN !1_6, 1_2, 6_10, 2_10
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(5,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(1,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(15,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(7,i)
                END DO
            END IF
            IF (k == 6) THEN !2_11, 2_10, 11_7, 10_7
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(6,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(7,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(19,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(18,i)
                END DO
            END IF
            IF (k == 7) THEN !3_11, 3_12, 11_7, 12_7
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(8,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(9,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(19,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(20,i)
                END DO
            END IF
            IF (k == 8) THEN !4_12, 4_8, 12_7, 8_7
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(10,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(11,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(20,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(16,i)
                END DO
            END IF
            IF (k == 9) THEN !5_8, 5_9, 8_7, 9_7
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(12,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(13,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(16,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(17,i)
                END DO
            END IF
            IF (k == 10) THEN !6_9, 6_10, 9_7, 10_7
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(14,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(15,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(17,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(18,i)
                END DO
            END IF

            DO lvlmk = 2, nlvl

                tpAx = xsphr(NdIDRgbsCll(k,lvlmk+1,1))
                tpAy = ysphr(NdIDRgbsCll(k,lvlmk+1,1))
                tpAz = zsphr(NdIDRgbsCll(k,lvlmk+1,1))
                tpBx = xsphr(NdIDRgbsCll(k,1,lvlmk+1))
                tpBy = ysphr(NdIDRgbsCll(k,1,lvlmk+1))
                tpBz = zsphr(NdIDRgbsCll(k,1,lvlmk+1))

                lvli = 0

                DO lvlnd = 1, lvlmk-1

                    pstni = lvlmk-(lvlnd-1)
                    pstnj = (lvlmk+2)-pstni

                    NdCnt = NdCnt + 1
                    NdIDRgbsCll(k,pstni,pstnj) = NdCnt

                    lvli = lvli + 1

                    tpi_n = DBLE(lvli)/DBLE(lvlmk)
                    CALL NdCodntssphr(tpAx, tpAy, tpAz, tpBx, tpBy, tpBz, &
                    &                 tpi_n, tpxcdnt, tpycdnt, tpzcdnt)
                    xsphr(NdIDRgbsCll(k,pstni,pstnj)) = tpxcdnt
                    ysphr(NdIDRgbsCll(k,pstni,pstnj)) = tpycdnt
                    zsphr(NdIDRgbsCll(k,pstni,pstnj)) = tpzcdnt

                END DO

            END DO
            DO lvlmk = nlvl-1, 2, -1

                tpAx = xsphr(NdIDRgbsCll(k,nlvl+1,nlvl+1-lvlmk))
                tpAy = ysphr(NdIDRgbsCll(k,nlvl+1,nlvl+1-lvlmk))
                tpAz = zsphr(NdIDRgbsCll(k,nlvl+1,nlvl+1-lvlmk))
                tpBx = xsphr(NdIDRgbsCll(k,nlvl+1-lvlmk,nlvl+1))
                tpBy = ysphr(NdIDRgbsCll(k,nlvl+1-lvlmk,nlvl+1))
                tpBz = zsphr(NdIDRgbsCll(k,nlvl+1-lvlmk,nlvl+1))

                lvli = 0

                DO lvlnd = 1, lvlmk-1

                    pstnj = (nlvl+1)-(lvlmk-lvlnd)
                    pstni = 2*(nlvl+1)-lvlmk-pstnj

                    NdCnt = NdCnt + 1
                    NdIDRgbsCll(k,pstni,pstnj) = NdCnt

                    lvli = lvli + 1

                    tpi_n = DBLE(lvli)/DBLE(lvlmk)
                    CALL NdCodntssphr(tpAx, tpAy, tpAz, tpBx, tpBy, tpBz, &
                    &                 tpi_n, tpxcdnt, tpycdnt, tpzcdnt)
                    xsphr(NdIDRgbsCll(k,pstni,pstnj)) = tpxcdnt
                    ysphr(NdIDRgbsCll(k,pstni,pstnj)) = tpycdnt
                    zsphr(NdIDRgbsCll(k,pstni,pstnj)) = tpzcdnt

                END DO

            END DO

        END DO

        OPEN (111, FILE = "Prtl_Orgnl.vrt", STATUS = "REPLACE")

        DO i = 1, 10*nlvl*nlvl+2
            WRITE (111, *) i, xsphr(i), ysphr(i), zsphr(i)
        END DO

        CLOSE (111)

        OPEN (121, FILE = "Prtl_Orgnl.cel", STATUS = "REPLACE")

        IF (dmMeshType == "L" .OR. dmMeshType == "C") THEN

            ElmtCnt = 0
            DO k =1, 10
                DO i = 1, nlvl
                    DO j = 1, nlvl
                        ElmtCnt = ElmtCnt + 1
                        WRITE (121, *) ElmtCnt, NdIDRgbsCll(k,i,j), &
                        &              NdIDRgbsCll(k,i,j+1), NdIDRgbsCll(k,i+1,j)
                        ElmtCnt = ElmtCnt + 1
                        WRITE (121, *) ElmtCnt, NdIDRgbsCll(k,i,j+1), &
                        &              NdIDRgbsCll(k,i+1,j+1), NdIDRgbsCll(k,i+1,j)
                    END DO
                END DO
            END DO

        END IF

        IF (dmMeshType == "Q") THEN

            ElmtCnt = 0
            DO k =1, 10
                DO i = 1, nlvl+1-2,2
                    DO j = 1, nlvl+1-2,2
                        ElmtCnt = ElmtCnt + 1
                        WRITE (121, *) ElmtCnt, &
                        &   NdIDRgbsCll(k,i,j), &
                        &   NdIDRgbsCll(k,i,j+2), &
                        &   NdIDRgbsCll(k,i+2,j), &
                        &   NdIDRgbsCll(k,i,j+1), &
                        &   NdIDRgbsCll(k,i+1,j+1), &
                        &   NdIDRgbsCll(k,i+1,j)
                        ElmtCnt = ElmtCnt + 1
                        WRITE (121, *) ElmtCnt, &
                        &   NdIDRgbsCll(k,i,j+2), &
                        &   NdIDRgbsCll(k,i+2,j+2), &
                        &   NdIDRgbsCll(k,i+2,j), &
                        &   NdIDRgbsCll(k,i+1,j+2), &
                        &   NdIDRgbsCll(k,i+2,j+1),&
                        &   NdIDRgbsCll(k,i+1,j+1)
                    END DO
                END DO
            END DO

        END IF

        CLOSE (121)

        OPEN (101, FILE = "Prtl_Orgnl.inp", STATUS = "REPLACE")

        WRITE (101, *)
        WRITE (101, *)
        WRITE (101, *) 10*nlvl*nlvl+2
        WRITE (101, *) ElmtCnt
        WRITE (101, *) 0.0d0, 0.0d0, 0.0d0

        CLOSE (101)

        DEALLOCATE (xsphr)
        DEALLOCATE (ysphr)
        DEALLOCATE (zsphr)
        DEALLOCATE (NdIDRgbsCll)
        DEALLOCATE (NdIDEg)
        DEALLOCATE (ElmtLnk1)
        DEALLOCATE (ElmtLnk2)
        DEALLOCATE (ElmtLnk3)

    END SUBROUTINE

    ! Return the unit-sphere point at fraction dmtpi_nlvl of the minor
    ! great-circle arc from unit vector A to unit vector B.
    SUBROUTINE NdCodntssphr(dmtpAx, dmtpAy, dmtpAz, dmtpBx, dmtpBy, dmtpBz, &
                            &dmtpi_nlvl, dmtpxcdnt, dmtpycdnt, dmtpzcdnt)

        DOUBLE PRECISION, INTENT(IN) :: dmtpAx,dmtpAy,dmtpAz,dmtpBx,dmtpBy,dmtpBz, &
        &   dmtpi_nlvl
        DOUBLE PRECISION, INTENT(OUT) :: dmtpxcdnt, dmtpycdnt, dmtpzcdnt

        DOUBLE PRECISION, DIMENSION (3,4):: GEclM
        DOUBLE PRECISION, DIMENSION (3) :: GErslt
        DOUBLE PRECISION, DIMENSION (4) :: GETp

        DOUBLE PRECISION :: ArcAB, tp

        INTEGER :: i, j, k, p, remp

        ArcAB = DACOS(DABS(dmtpAx*dmtpBx+dmtpAy*dmtpBy+dmtpAz*dmtpBz))

        GEclM(1,4) = DCOS(dmtpi_nlvl*ArcAB)
        GEclM(2,4) = DCOS((1.0d0-dmtpi_nlvl)*ArcAB)
        GEclM(3,4) = 0.0d0

        GEclM(1,1) = dmtpAx
        GEclM(1,2) = dmtpAy
        GEclM(1,3) = dmtpAz

        GEclM(2,1) = dmtpBx
        GEclM(2,2) = dmtpBy
        GEclM(2,3) = dmtpBz

        GEclM(3,1) = (dmtpAy*dmtpBz - dmtpAz*dmtpBy)
        GEclM(3,2) = (dmtpAz*dmtpBx - dmtpAx*dmtpBz)
        GEclM(3,3) = (dmtpAx*dmtpBy - dmtpAy*dmtpBx)

        DO j=1,3-1
            tp=DABS(GEclM(j,j))
            remp = j
            DO k=j+1,3
                IF (tp<DABS(GEclM(k,j))) THEN
                    tp=DABS(GEclM(k,j))
                    remp=k
                END IF
            END DO
            IF (remp>j) THEN
                GETp(:)=GEclM(remp,:)
                GEclM(remp,:)=GEclM(j,:)
                GEclM(j,:)=GETp(:)
            END IF
            DO i=j+1,3
                GEclM(i,:)=GEclM(i,:)-GEclM(j,:)*GEclM(i,j)/GEclM(j,j)
            END DO
        END DO

        GErslt(3)=GEclM(3,3+1)/GEclM(3,3)
        DO i=3-1, 1, -1
            tp=0.0d0
            DO j=i+1,3
                tp=tp+GEclM(i,j)*GErslt(j)
            END DO
            GErslt(i)=(GEclM(i,3+1)-tp)/GEclM(i,i)
        END DO

        tp = DSQRT(GErslt(1)**2+GErslt(2)**2+GErslt(3)**2)
        dmtpxcdnt = GErslt(1)/tp
        dmtpycdnt = GErslt(2)/tp
        dmtpzcdnt = GErslt(3)/tp

    END SUBROUTINE

    ! Generate a sphere by joining two mirrored hemispheres.  Each hemisphere
    ! comprises three rotated patches whose levels project concentric circles
    ! onto the sphere.  dmRlvlStpScl is the geometric ratio of successive
    ! polar-angle steps (1 gives uniform angular spacing).
    SUBROUTINE GetMeshNdElmntTwoCrcl(dmnlvlprtl, dmRlvlStpScl, dmMeshType)

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: dmnlvlprtl

        DOUBLE PRECISION, INTENT (IN) :: dmRlvlStpScl

        CHARACTER (LEN=1), INTENT (IN) :: dmMeshType

        DOUBLE PRECISION :: SphereMesh_pi = DATAN(1.0d0)*4.0d0

        INTEGER :: i,j, k, ndrd, ndlist, ElmntCnt, NmbrMeridianNd, MeridianNdNo, &
        &   ndrdbtm, nlvlsphr

        DOUBLE PRECISION, PARAMETER :: rsdl = 1.0E-5

        DOUBLE PRECISION :: RlvlStpInt

        DOUBLE PRECISION :: Rz, p23, p13, p16, theta, tp1, tp2, Rsphr

        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: Rlvl

        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: x1, y1, z1, x2, y2, z2, &
        &   x3, y3, z3

        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: xtop, ytop, ztop
        INTEGER, ALLOCATABLE, DIMENSION (:,:) :: ndtopindex1, ndtopindex2, ndtopindex3

        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: xbtm, ybtm, zbtm
        INTEGER, ALLOCATABLE, DIMENSION (:,:) :: ndbtmindex1, ndbtmindex2, ndbtmindex3

        INTEGER :: cntctRingIndexNo, znlv

        INTEGER, ALLOCATABLE, DIMENSION(:) :: MeridianNdindex

        DOUBLE PRECISION :: fR00, fR01, fR02, xR00, xR01, xR02

        p23 = 2.0d0*SphereMesh_pi/3.0d0
        p13 = 1.0d0*SphereMesh_pi/3.0d0
        p16 = 1.0d0*SphereMesh_pi/6.0d0

        ! Radii of the projected circles from the pole to the equator.

        Rsphr = 1.0d0
        nlvlsphr = 2*dmnlvlprtl

        xR00 = 0.0d0
        xR01 = 0.5d0*SphereMesh_pi

        fR00 = 0.0d0
        DO k = 1, nlvlsphr
            fR00 = fR00 + xR00*(dmRlvlStpScl**(k-1))
        END DO
        fR00 = 0.5d0*SphereMesh_pi - fR00

        fR01 = 0.0d0
        DO k = 1, nlvlsphr
            fR01 = fR01 + xR01*(dmRlvlStpScl**(k-1))
        END DO
        fR01 = 0.5d0*SphereMesh_pi - fR01

        i = 1
        DO
            xR02 = (xR00 + xR01)*0.5d0
            fR02 = 0.0d0
            DO k = 1, nlvlsphr
                fR02 = fR02 + xR02*(dmRlvlStpScl**(k-1))
            END DO
            fR02 = 0.5d0*SphereMesh_pi - fR02

            IF (fR00*fR02<0.0) THEN
                xR01 = xR02
                fR01 = fR02
            ELSE
                xR00 = xR02
                fR00 = fR02
            END IF

            IF (DABS(fR02) < 1.e-12) THEN
                RlvlStpInt = xR02
                EXIT
            END IF
            i = i + 1
            IF (i > 100000) THEN
                PRINT*, "No R0_bbl obtained! Please check!!"
                STOP
            END IF
        END DO

        ALLOCATE (Rlvl(nlvlsphr))
        fR00 = 0.0d0
        DO  i = 1, nlvlsphr
            fR00 = fR00 + RlvlStpInt*(dmRlvlStpScl**(i-1))
            Rlvl(i) = Rsphr*DSIN(fR00)
        END DO

        IF (dmMeshType == "Q") THEN
            DO  i = 1, nlvlsphr-1, 2
                tp1 = DASIN(Rlvl(i+1))
                IF (i == 1) THEN
                    tp2 = 0.0d0
                ELSE
                    tp2 = DASIN(Rlvl(i-1))
                END IF
                Rlvl(i) = DSIN(0.5d0*(tp1+tp2))
            END DO
        END IF

        ! Coordinates and shared node indices on the upper hemisphere.

        ALLOCATE (x1(nlvlsphr+1,nlvlsphr+1))
        ALLOCATE (y1(nlvlsphr+1,nlvlsphr+1))
        ALLOCATE (z1(nlvlsphr+1,nlvlsphr+1))

        ALLOCATE (ndtopindex1(nlvlsphr+1,nlvlsphr+1))

        ALLOCATE (x2(nlvlsphr+1,nlvlsphr+1))
        ALLOCATE (y2(nlvlsphr+1,nlvlsphr+1))
        ALLOCATE (z2(nlvlsphr+1,nlvlsphr+1))

        ALLOCATE (ndtopindex2(nlvlsphr+1,nlvlsphr+1))

        ALLOCATE (x3(nlvlsphr+1,nlvlsphr+1))
        ALLOCATE (y3(nlvlsphr+1,nlvlsphr+1))
        ALLOCATE (z3(nlvlsphr+1,nlvlsphr+1))

        ALLOCATE (ndtopindex3(nlvlsphr+1,nlvlsphr+1))

        x1(1,1) = 0.0d0
        y1(1,1) = 0.0d0
        z1(1,1) = Rsphr

        DO i = 2, nlvlsphr+1
            DO j = 1, i
                theta = 0.0d0 + (DBLE(j-1))*p13/(DBLE(i-1))
                x1(i,j) = RLvl(i-1)*DCOS(theta)
                y1(i,j) = RLvl(i-1)*DSIN(theta)
                IF (DABS(Rsphr-RLvl(i-1)) < rsdl) THEN
                    z1(i,j) = 0.0d0
                ELSE
                    z1(i,j) = DSQRT(Rsphr**2 - RLvl(i-1)**2)
                END IF
            END DO
        END DO

        DO j = 2, nlvlsphr+1
            DO i = 1, j-1
                theta = p23 - (DBLE(i-1))*p13/(DBLE(j-1))
                x1(i,j) = RLvl(j-1)*DCOS(theta)
                y1(i,j) = RLvl(j-1)*DSIN(theta)
                IF (DABS(Rsphr - RLvl(j-1)) < rsdl) THEN
                    z1(i,j) = 0.0d0
                ELSE
                    z1(i,j) = DSQRT(Rsphr**2 - RLvl(j-1)**2)
                END IF
            END DO
        END DO

        DO i = 1, nlvlsphr+1
            DO j = 1, nlvlsphr+1
                IF (j == 1) THEN
                    x2(i,1) = x1(1,i)
                    y2(i,1) = y1(1,i)
                ELSE
                    x2(i,j) = -DCOS(p13)*x1(i,j) - DSIN(p13)*y1(i,j)
                    y2(i,j) = DSIN(p13)*x1(i,j) - DCOS(p13)*y1(i,j)
                END IF
                z2(i,j) = z1(i,j)
            END DO
        END DO

        DO i = 1, nlvlsphr+1
            IF (i == 1) THEN
                DO j = 1, nlvlsphr+1
                    x3(1,j) = x1(j,1)
                    y3(1,j) = y1(j,1)
                    z3(i,j) = z1(i,j)
                END DO
            ELSE
                DO j = 1, nlvlsphr+1
                    IF (j == 1) THEN
                        x3(i,1) = x2(1,i)
                        y3(i,1) = y2(1,i)
                    ELSE
                        x3(i,j) = -DCOS(p13)*x2(i,j) - DSIN(p13)*y2(i,j)
                        y3(i,j) = DSIN(p13)*x2(i,j) - DCOS(p13)*y2(i,j)
                    END IF
                    z3(i,j) = z1(i,j)
                END DO
            END IF
        END DO

        ndrd = 0
        DO i = 1, nlvlsphr+1
            DO j = 1, nlvlsphr+1
                ndrd = ndrd + 1
                ndtopindex1(i,j) = ndrd
            END DO
        END DO

        DO i = 1, nlvlsphr+1
            DO j = 1, nlvlsphr+1
                IF (j == 1) THEN
                    ndtopindex2(i,1) = ndtopindex1(1,i)
                ELSE
                    ndrd = ndrd + 1
                    ndtopindex2(i,j) = ndrd
                END IF
            END DO
        END DO

        DO i = 1, nlvlsphr+1
            IF (i == 1) THEN
                DO j = 1, nlvlsphr+1
                    ndtopindex3(1,j) = ndtopindex1(j,1)
                END DO
            ELSE
                DO j = 1, nlvlsphr+1
                    IF (j == 1) THEN
                        ndtopindex3(i,1) = ndtopindex2(1,i)
                    ELSE
                        ndrd = ndrd + 1
                        ndtopindex3(i,j) = ndrd
                    END IF
                END DO
            END IF
        END DO

        ! Flatten the three upper-hemisphere patch arrays into node-ID order.

        ALLOCATE (xtop(ndrd))
        ALLOCATE (ytop(ndrd))
        ALLOCATE (ztop(ndrd))

        ND: DO ndlist  = 1, ndrd

            DO i = 1, nlvlsphr+1
                DO j = 1, nlvlsphr+1
                    IF (ndtopindex1(i,j) == ndlist) THEN
                        xtop(ndtopindex1(i,j)) = x1(i,j)
                        ytop(ndtopindex1(i,j)) = y1(i,j)
                        ztop(ndtopindex1(i,j)) = z1(i,j)
                        CYCLE ND
                    END IF
                END DO
            END DO

            DO i = 1, nlvlsphr+1
                DO j = 1, nlvlsphr+1
                    IF (j/=1) THEN
                        IF (ndtopindex2(i,j) == ndlist) THEN
                            xtop(ndtopindex2(i,j)) = x2(i,j)
                            ytop(ndtopindex2(i,j)) = y2(i,j)
                            ztop(ndtopindex2(i,j)) = z2(i,j)
                            CYCLE ND
                        END IF
                    END IF
                END DO
            END DO

            DO i = 1, nlvlsphr+1
                IF (i/=1) THEN
                    DO j = 1, nlvlsphr+1
                        IF (j/=1) THEN
                            IF (ndtopindex3(i,j) == ndlist) THEN
                                xtop(ndtopindex3(i,j)) = x3(i,j)
                                ytop(ndtopindex3(i,j)) = y3(i,j)
                                ztop(ndtopindex3(i,j)) = z3(i,j)
                                CYCLE ND
                            END IF
                        END IF
                    END DO
                END IF
            END DO

        END DO ND

        ! Reuse the equatorial ring when reflecting nodes into the lower hemisphere.

        MeridianNdNo = 0

        NmbrMeridianNd = 6*nlvlsphr

        ALLOCATE (MeridianNdindex(6*nlvlsphr))

        i = nlvlsphr+1
        DO j = 1, nlvlsphr+1
            MeridianNdNo = MeridianNdNo + 1
            MeridianNdindex(MeridianNdNo) = ndtopindex1(i,j)
        END DO

        j = nlvlsphr+1
        DO i = nlvlsphr, 1, -1
            MeridianNdNo = MeridianNdNo + 1
            MeridianNdindex(MeridianNdNo) = ndtopindex1(i,j)
        END DO

        i = nlvlsphr+1
        DO j = 2, nlvlsphr+1
            MeridianNdNo = MeridianNdNo + 1
            MeridianNdindex(MeridianNdNo) = ndtopindex2(i,j)
        END DO

        j = nlvlsphr+1
        DO i = nlvlsphr, 1, -1
            MeridianNdNo = MeridianNdNo + 1
            MeridianNdindex(MeridianNdNo) = ndtopindex2(i,j)
        END DO

        i = nlvlsphr+1
        DO j = 2, nlvlsphr+1
            MeridianNdNo = MeridianNdNo + 1
            MeridianNdindex(MeridianNdNo) = ndtopindex3(i,j)
        END DO

        j = nlvlsphr+1
        DO i = nlvlsphr, 2, -1
            MeridianNdNo = MeridianNdNo + 1
            MeridianNdindex(MeridianNdNo) = ndtopindex3(i,j)
        END DO

        ALLOCATE (xbtm(ndrd - MeridianNdNo))
        ALLOCATE (ybtm(ndrd - MeridianNdNo))
        ALLOCATE (zbtm(ndrd - MeridianNdNo))

        ALLOCATE (ndbtmindex1(nlvlsphr+1,nlvlsphr+1))
        ALLOCATE (ndbtmindex2(nlvlsphr+1,nlvlsphr+1))
        ALLOCATE (ndbtmindex3(nlvlsphr+1,nlvlsphr+1))

        ndrdbtm = 0
        DO i = 1, nlvlsphr+1
            DO j = 1, nlvlsphr+1
                IF (i == nlvlsphr+1 .or. j == nlvlsphr+1) THEN
                    ndbtmindex1(i,j) = ndtopindex1(i,j)
                ELSE
                    ndrdbtm = ndrdbtm + 1
                    xbtm(ndrdbtm) = xtop(ndtopindex1(i,j))
                    ybtm(ndrdbtm) = ytop(ndtopindex1(i,j))
                    zbtm(ndrdbtm) = -ztop(ndtopindex1(i,j))
                    ndbtmindex1(i,j) = ndrdbtm + ndrd
                END IF
            END DO
        END DO

        DO i = 1, nlvlsphr+1
            DO j = 1, nlvlsphr+1
                IF (j == 1) THEN
                    ndbtmindex2(i,1) = ndbtmindex1(1,i)
                ELSE IF (j/=1) THEN
                    IF (i == nlvlsphr+1 .or. j == nlvlsphr+1) THEN
                        ndbtmindex2(i,j) = ndtopindex2(i,j)
                    ELSE
                        ndrdbtm = ndrdbtm + 1
                        xbtm(ndrdbtm) = xtop(ndtopindex2(i,j))
                        ybtm(ndrdbtm) = ytop(ndtopindex2(i,j))
                        zbtm(ndrdbtm) = -ztop(ndtopindex2(i,j))
                        ndbtmindex2(i,j) = ndrdbtm + ndrd
                    END IF
                END IF
            END DO
        END DO

        DO i = 1, nlvlsphr+1
            IF (i == 1) THEN
                DO j = 1, nlvlsphr+1
                    ndbtmindex3(1,j) = ndbtmindex1(j,1)
                END DO
            ELSE
                DO j = 1, nlvlsphr+1
                    IF (j == 1) THEN
                        ndbtmindex3(i,1) = ndbtmindex2(1,i)
                    ELSE IF (j/=1) THEN
                        IF (i == nlvlsphr+1 .or. j == nlvlsphr+1) THEN
                            ndbtmindex3(i,j) = ndtopindex3(i,j)
                        ELSE
                            ndrdbtm = ndrdbtm + 1
                            xbtm(ndrdbtm) = xtop(ndtopindex3(i,j))
                            ybtm(ndrdbtm) = ytop(ndtopindex3(i,j))
                            zbtm(ndrdbtm) = -ztop(ndtopindex3(i,j))
                            ndbtmindex3(i,j) = ndrdbtm + ndrd
                        END IF
                    END IF
                END DO
            END IF
        END DO

        OPEN (111, FILE = "Prtl_Orgnl.vrt", STATUS = "REPLACE")

        DO i = 1, ndrd
            WRITE (111, *) i, xtop(i), ytop(i), ztop(i)
        END DO
        DO i = 1, ndrd - MeridianNdNo
            WRITE (111, *) ndrd + i, xbtm(i), ybtm(i), zbtm(i)
        END DO
        CLOSE (111)

        ! Write triangle connectivity in the module-level L/C or Q6 node order.

        OPEN (121, FILE = "Prtl_Orgnl.cel", STATUS = "REPLACE")

        IF (dmMeshType == "L" .OR. dmMeshType == "C") THEN

            ElmntCnt = 0
            DO i = 1, nlvlsphr+1-1
                DO j = 1, nlvlsphr+1-1
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt, ndtopindex1(i,j), ndtopindex1(i,j+1), ndtopindex1(i+1,j+1)
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt, ndtopindex1(i,j), ndtopindex1(i+1,j+1), ndtopindex1(i+1,j)
                END DO
            END DO

            DO i = 1, nlvlsphr+1-1
                DO j = 1, nlvlsphr+1-1
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt, ndtopindex2(i,j), ndtopindex2(i,j+1), ndtopindex2(i+1,j+1)
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt, ndtopindex2(i,j), ndtopindex2(i+1,j+1), ndtopindex2(i+1,j)
                END DO
            END DO

            DO i = 1, nlvlsphr+1-1
                DO j = 1, nlvlsphr+1-1
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt, ndtopindex3(i,j), ndtopindex3(i,j+1), ndtopindex3(i+1,j+1)
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt, ndtopindex3(i,j), ndtopindex3(i+1,j+1), ndtopindex3(i+1,j)
                END DO
            END DO

            DO i = 1, nlvlsphr+1-1
                DO j = 1, nlvlsphr+1-1
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt, ndbtmindex1(i,j), ndbtmindex1(i,j+1), ndbtmindex1(i+1,j+1)
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt, ndbtmindex1(i,j), ndbtmindex1(i+1,j+1), ndbtmindex1(i+1,j)
                END DO
            END DO

            DO i = 1, nlvlsphr+1-1
                DO j = 1, nlvlsphr+1-1
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt, ndbtmindex2(i,j), ndbtmindex2(i,j+1), ndbtmindex2(i+1,j+1)
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt, ndbtmindex2(i,j), ndbtmindex2(i+1,j+1), ndbtmindex2(i+1,j)
                END DO
            END DO

            DO i = 1, nlvlsphr+1-1
                DO j = 1, nlvlsphr+1-1
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt, ndbtmindex3(i,j), ndbtmindex3(i,j+1), ndbtmindex3(i+1,j+1)
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt, ndbtmindex3(i,j), ndbtmindex3(i+1,j+1), ndbtmindex3(i+1,j)
                END DO
            END DO

        END IF

        IF (dmMeshType == "Q") THEN
            ElmntCnt = 0
            DO i = 1, nlvlsphr+1-2, 2
                DO j = 1, nlvlsphr+1-2, 2
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt, ndtopindex1(i,j),ndtopindex1(i,j+2),ndtopindex1(i+2,j+2), &
                    &ndtopindex1(i,j+1), ndtopindex1(i+1,j+2), ndtopindex1(i+1,j+1)
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt, ndtopindex1(i,j),ndtopindex1(i+2,j+2),ndtopindex1(i+2,j),&
                    &ndtopindex1(i+1,j+1), ndtopindex1(i+2,j+1), ndtopindex1(i+1,j)
                END DO
            END DO

            DO i = 1, nlvlsphr+1-2, 2
                DO j = 1, nlvlsphr+1-2, 2
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt,ndtopindex2(i,j),ndtopindex2(i,j+2),ndtopindex2(i+2,j+2), &
                    &ndtopindex2(i,j+1), ndtopindex2(i+1,j+2), ndtopindex2(i+1,j+1)
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt,ndtopindex2(i,j),ndtopindex2(i+2,j+2),ndtopindex2(i+2,j), &
                    &ndtopindex2(i+1,j+1),ndtopindex2(i+2,j+1),ndtopindex2(i+1,j)
                END DO
            END DO

            DO i = 1, nlvlsphr+1-2, 2
                DO j = 1, nlvlsphr+1-2, 2
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt,ndtopindex3(i,j),ndtopindex3(i,j+2),ndtopindex3(i+2,j+2), &
                    &ndtopindex3(i,j+1), ndtopindex3(i+1,j+2), ndtopindex3(i+1,j+1)
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt,ndtopindex3(i,j),ndtopindex3(i+2,j+2),ndtopindex3(i+2,j), &
                    &ndtopindex3(i+1,j+1), ndtopindex3(i+2,j+1), ndtopindex3(i+1,j)
                END DO
            END DO

            DO i = 1, nlvlsphr+1-2, 2
                DO j = 1, nlvlsphr+1-2, 2
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt,ndbtmindex1(i,j),ndbtmindex1(i,j+2),ndbtmindex1(i+2,j+2), &
                    &ndbtmindex1(i,j+1), ndbtmindex1(i+1,j+2), ndbtmindex1(i+1,j+1)
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt,ndbtmindex1(i,j),ndbtmindex1(i+2,j+2),ndbtmindex1(i+2,j), &
                    &ndbtmindex1(i+1,j+1), ndbtmindex1(i+2,j+1), ndbtmindex1(i+1,j)
                END DO
            END DO

            DO i = 1, nlvlsphr+1-2, 2
                DO j = 1, nlvlsphr+1-2, 2
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt,ndbtmindex2(i,j),ndbtmindex2(i,j+2),ndbtmindex2(i+2,j+2), &
                    &ndbtmindex2(i,j+1), ndbtmindex2(i+1,j+2), ndbtmindex2(i+1,j+1)
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt,ndbtmindex2(i,j),ndbtmindex2(i+2,j+2),ndbtmindex2(i+2,j), &
                    &ndbtmindex2(i+1,j+1), ndbtmindex2(i+2,j+1), ndbtmindex2(i+1,j)
                END DO
            END DO

            DO i = 1, nlvlsphr+1-2, 2
                DO j = 1, nlvlsphr+1-2, 2
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt,ndbtmindex3(i,j),ndbtmindex3(i,j+2),ndbtmindex3(i+2,j+2), &
                    &ndbtmindex3(i,j+1),ndbtmindex3(i+1,j+2),ndbtmindex3(i+1,j+1)
                    ElmntCnt = ElmntCnt + 1
                    WRITE (121, *) &
                    &ElmntCnt,ndbtmindex3(i,j),ndbtmindex3(i+2,j+2),ndbtmindex3(i+2,j), &
                    &ndbtmindex3(i+1,j+1),ndbtmindex3(i+2,j+1),ndbtmindex3(i+1,j)
                END DO
            END DO
        END IF


        CLOSE (121)

        OPEN (101, FILE = "Prtl_Orgnl.inp", STATUS = "REPLACE")

        WRITE (101, *)
        WRITE (101, *)
        WRITE (101, *) 2*ndrd - MeridianNdNo
        WRITE (101, *) ElmntCnt
        WRITE (101, *) 0.0d0, 0.0d0, 0.0d0

        CLOSE (101)

        DEALLOCATE(Rlvl)
        DEALLOCATE(x1)
        DEALLOCATE(y1)
        DEALLOCATE(z1)
        DEALLOCATE(x2)
        DEALLOCATE(y2)
        DEALLOCATE(z2)
        DEALLOCATE(x3)
        DEALLOCATE(y3)
        DEALLOCATE(z3)
        DEALLOCATE(xtop)
        DEALLOCATE(ytop)
        DEALLOCATE(ztop)
        DEALLOCATE(xbtm)
        DEALLOCATE(ybtm)
        DEALLOCATE(zbtm)
        DEALLOCATE(ndtopindex1)
        DEALLOCATE(ndtopindex2)
        DEALLOCATE(ndtopindex3)
        DEALLOCATE(ndbtmindex1)
        DEALLOCATE(ndbtmindex2)
        DEALLOCATE(ndbtmindex3)
        DEALLOCATE(MeridianNdindex)

    END SUBROUTINE




    ! Solve the scalar long-rod constraint for dmca by bracketed bisection.
    SUBROUTINE SearchLgRd_cvsa(dmi, dmba,dmca)

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: dmi
        DOUBLE PRECISION, INTENT(IN) :: dmba
        DOUBLE PRECISION, INTENT(OUT) :: dmca
        DOUBLE PRECISION :: fR00, fR01, fR02, xR00, xR01, xR02
        INTEGER :: i

        xR00 = 0.000001d0
        xR01 = 1.01d0

        fR00 = (1.0d0 - xR00**2)**2 - 1.0d0*dmba**2 * (dmba**2 + xR00**2)**0.5d0
        fR01 = (1.0d0 - xR01**2)**2 - 1.0d0*dmba**2 * (dmba**2 + xR01**2)**0.5d0

        i = 1
        DO
            xR02 = (xR00 + xR01)*0.5d0
            fR02 =  (1.0d0 - xR02**2)**2 - 1.0d0*dmba**2 * (dmba**2 + xR02**2)**0.5d0

            IF (fR00*fR02<0.0d0) THEN
                xR01 = xR02
                fR01 = fR02
            ELSE
                xR00 = xR02
                fR00 = fR02
            END IF

            IF (DABS(fR02) < 1.0E-12) THEN
                dmca = xR02
                EXIT
            END IF
            i = i + 1
            IF (i > 1E7) THEN
                PRINT*, dmi, "No 'cvsa' obtained! Please check!!"
                STOP
            END IF
        END DO

    END SUBROUTINE

    ! Solve the scalar long-rod radial-scale constraint by bracketed bisection.
    SUBROUTINE SearchLgRd_rdbscl(dmi, dmx, dmca, dmrdbscl)

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: dmi
        DOUBLE PRECISION, INTENT (IN) :: dmx, dmca
        DOUBLE PRECISION, INTENT (OUT) :: dmrdbscl
        DOUBLE PRECISION :: fR00, fR01, fR02, xR00, xR01, xR02
        INTEGER :: i

        xR00 =  0.001d0
        xR01 =  1.001d0

        fR00 = ((1.0d0-dmca**2)**2)*( (dmx+dmca)/DSQRT((dmx+dmca)**2 + xR00**2)&
            &                        -(dmx-dmca)/DSQRT((dmx-dmca)**2 + xR00**2) ) &
            & -2.0d0*1.0d0*dmca*(xR00**2)
        fR01 = ((1.0d0-dmca**2)**2)*( (dmx+dmca)/DSQRT((dmx+dmca)**2 + xR01**2)&
            &                        -(dmx-dmca)/DSQRT((dmx-dmca)**2 + xR01**2) ) &
            & -2.0d0*1.0d0*dmca*(xR01**2)

        i = 1
        DO
            xR02 = (xR00 + xR01)*0.5d0
            fR02 = ((1.0d0-dmca**2)**2)*( (dmx+dmca)/DSQRT((dmx+dmca)**2 + xR02**2)&
                &                        -(dmx-dmca)/DSQRT((dmx-dmca)**2 + xR02**2) ) &
                & -2.0d0*1.0d0*dmca*(xR02**2)

            IF (fR00*fR02<0.0d0) THEN
                xR01 = xR02
                fR01 = fR02
            ELSE
                xR00 = xR02
                fR00 = fR02
            END IF

            IF (DABS(fR02) < 1.0E-12) THEN
                dmrdbscl = xR02
                EXIT
            END IF
            i = i + 1
            IF (i > 1E7) THEN
                PRINT*, dmi, "No 'dmrdbscl' obtained! Please check!!"
                STOP
            END IF
        END DO

    END SUBROUTINE

    ! Solve the scalar dumbbell constraint for dmca by bracketed bisection.
    SUBROUTINE SearchDbBl_cvsa(dmi, dmda, dmca)

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: dmi
        DOUBLE PRECISION, INTENT(IN) :: dmda
        DOUBLE PRECISION, INTENT(OUT) :: dmca
        DOUBLE PRECISION :: fR00, fR01, fR02, xR00, xR01, xR02
        INTEGER :: i

        xR00 = 0.0001d0
        xR01 = 1.0001d0

        fR00 = (1.0d0 - xR00**2)**3 &
        &     -(1.0d0 + 3.0d0*xR00**2)*(xR00**2 + dmda**2)**(3.0d0/2.0d0)
        fR01 = (1.0d0 - xR01**2)**3 &
        &     -(1.0d0 + 3.0d0*xR01**2)*(xR01**2 + dmda**2)**(3.0d0/2.0d0)

        i = 1
        DO
            xR02 = (xR00 + xR01)*0.5d0
            fR02 = (1.0d0 - xR02**2)**3 &
            &     -(1.0d0 + 3.0d0*xR02**2)*(xR02**2 + dmda**2)**(3.0d0/2.0d0)

            IF (fR00*fR02<0.0d0) THEN
                xR01 = xR02
                fR01 = fR02
            ELSE
                xR00 = xR02
                fR00 = fR02
            END IF

            IF (DABS(fR02) < 1.0E-12) THEN
                dmca = xR02
                EXIT
            END IF
            i = i + 1
            IF (i > 1E7) THEN
                PRINT*, dmi, "No 'cvsa' obtained! Please check!!"
                STOP
            END IF
        END DO

    END SUBROUTINE

    ! Solve the scalar dumbbell constraint for dmba by bracketed bisection.
    SUBROUTINE SearchDbBl_bvsa(dmi, dmda, dmca, dmba)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: dmi
        DOUBLE PRECISION, INTENT(IN) :: dmda, dmca
        DOUBLE PRECISION, INTENT(OUT) :: dmba
        DOUBLE PRECISION :: fR00, fR01, fR02, xR00, xR01, xR02
        INTEGER :: i

        xR00 =  0.0001d0
        xR01 =  1.0001d0

        fR00 = (DSQRT((dmca+dmca)**2 + xR00**2))**(-3) + xR00**(-3) &
            & - 2.0d0*(DSQRT(dmca**2 + dmda**2))**(-3)
        fR01 = (DSQRT((dmca+dmca)**2 + xR01**2))**(-3) + xR01**(-3) &
            & - 2.0d0*(DSQRT(dmca**2 + dmda**2))**(-3)

        i = 1
        DO
            xR02 = (xR00 + xR01)*0.5d0
            fR02 = (DSQRT((dmca+dmca)**2 + xR02**2))**(-3) + xR02**(-3) &
                & - 2.0d0*(DSQRT(dmca**2 + dmda**2))**(-3)

            IF (fR00*fR02<0.0d0) THEN
                xR01 = xR02
                fR01 = fR02
            ELSE
                xR00 = xR02
                fR00 = fR02
            END IF

            IF (DABS(fR02) < 1.0E-12) THEN
                dmba = xR02
                EXIT
            END IF
            i = i + 1
            IF (i > 1E7) THEN
                PRINT*, dmi, "No 'bvsa' obtained! Please check!!"
                STOP
            END IF
        END DO

    END SUBROUTINE


    ! Solve the scalar dumbbell radial-scale constraint by bracketed bisection.
    SUBROUTINE SearchDbBl_rdbscl(dmi, dmx, dmda, dmca, dmrdbscl)

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: dmi
        DOUBLE PRECISION, INTENT (IN) :: dmx, dmda, dmca
        DOUBLE PRECISION, INTENT (OUT) :: dmrdbscl
        DOUBLE PRECISION :: fR00, fR01, fR02, xR00, xR01, xR02
        INTEGER :: i

        xR00 =  0.0001d0
        xR01 =  1.0001d0

        fR00 = (DSQRT((dmx+dmca)**2 + xR00**2))**(-3) &
        &     +(DSQRT((dmx-dmca)**2 + xR00**2))**(-3) &
        &     -2.0d0*(DSQRT(dmca**2 + dmda**2))**(-3)
        fR01 = (DSQRT((dmx+dmca)**2 + xR01**2))**(-3) &
        &     +(DSQRT((dmx-dmca)**2 + xR01**2))**(-3) &
        &     -2.0d0*(DSQRT(dmca**2 + dmda**2))**(-3)

        i = 1
        DO
            xR02 = (xR00 + xR01)*0.5d0
            fR02 = (DSQRT((dmx+dmca)**2 + xR02**2))**(-3) &
            &     +(DSQRT((dmx-dmca)**2 + xR02**2))**(-3) &
            &     -2.0d0*(DSQRT(dmca**2 + dmda**2))**(-3)

            IF (fR00*fR02<0.0d0) THEN
                xR01 = xR02
                fR01 = fR02
            ELSE
                xR00 = xR02
                fR00 = fR02
            END IF

            IF (DABS(fR02) < 1.0E-12) THEN
                dmrdbscl = xR02
                EXIT
            END IF
            i = i + 1
            IF (i > 1E7) THEN
                PRINT*,dmi, "No 'dmrdbscl' obtained! Please check!!"
                STOP
            END IF
        END DO

    END SUBROUTINE



END MODULE
