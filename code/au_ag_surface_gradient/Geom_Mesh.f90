
! SPDX-FileCopyrightText: 2026 Qiang Sun
! SPDX-License-Identifier: BSD-3-Clause

! Build the global surface mesh used by the boundary-integral solver.
! The module either generates supported analytic surfaces or imports Prtl_*.msh,
! applies object transforms, and constructs node-element adjacency for L3 or Q6
! triangular elements. Global node and element offsets follow object order.
!
MODULE Geom_Mesh

    USE omp_lib

    USE Pre_Constants

    USE Geom_GlobalData
    USE Geom_MeshSphereCircle

    IMPLICIT NONE

    DOUBLE PRECISION :: dmMeshRlvlStp

    CONTAINS

    ! Generate or import every object mesh, apply scale/rotation/translation,
    ! orient its triangles according to NrmlInOut, and populate global arrays.

    SUBROUTINE Meshediting

        CHARACTER (LEN=50) :: filename

        LOGICAL :: filexist

        DOUBLE PRECISION :: xrefnn, yrefnn, zrefnn

        INTEGER :: i, j, k, ithprtl, ndindex, elmntindex
        INTEGER :: ndA, ndB, ndC, ndD, ndE, ndF, IOS

        INTEGER :: nd_i, nd_j, nd_k, elmnt_i, elmnt_j, elmnt_k

        INTEGER :: ndoffset1, elmntoffset1

        INTEGER :: thti, tpswpndid

        INTEGER :: nmbrndsznd, nmbrelmntsznd

        INTEGER :: ii,jj,iflag,imove,icount,ttlhostlist

        INTEGER, ALLOCATABLE, DIMENSION (:) :: reorder,hostlist,hostlistrcd
        INTEGER, ALLOCATABLE, DIMENSION (:,:) :: rcrelmntlnknd
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: rcdxnd,rcdynd,rcdznd

        DOUBLE PRECISION :: tpvctACx, tpvctACy, tpvctACz, &
        &                   tpvctBCx, tpvctBCy, tpvctBCz, &
        &                   tpvctrefCx, tpvctrefCy, tpvctrefCz, &
        &                   tpnrmlchck, rdbscl, tpxyscl, tpyzscl, &
        &                   tpx, tpy, tpz, tp1, tp2, tp3

        DOUBLE PRECISION :: tpswp, tpthtx, tpthty, tpthtz

        ttlnmbrnd = 0
        ttlnmbrelmnt = 0

        DO ithprtl = 1, nmbrprtl

            IF (MeshRead(ithprtl) == 0) THEN

                IF (DABS(MeshRlvlStp(ithprtl)-1.0d0)<0.01d0) THEN
                    IF (MeshGnrtn(ithprtl) == "Icshdrl") &
                        & CALL GetMeshNdElmntIcshdrlPaper &
                            &(Meshnlvl(ithprtl),MeshRlvlStp(ithprtl),MeshType)
                    IF (MeshGnrtn(ithprtl) == "Icshsim") &
                        & CALL GetMeshNdElmntIcshdrlSimple &
                            &(Meshnlvl(ithprtl),MeshRlvlStp(ithprtl),MeshType)
                    IF (MeshGnrtn(ithprtl) == "TwoCrcl") &
                        & CALL GetMeshNdElmntTwoCrcl &
                            &(Meshnlvl(ithprtl),MeshRlvlStp(ithprtl),MeshType)
                ELSE
                    CALL GetMeshNdElmntTwoCrcl &
                        &(Meshnlvl(ithprtl),MeshRlvlStp(ithprtl),MeshType)
                END IF

                OPEN (101, FILE = "Prtl_Orgnl.inp", STATUS = "OLD", IOSTAT = IOS)
                IF (IOS /= 0) THEN
                    PRINT*, '"Prtl_Orgnl.inp" does not exist! Please check!'
                    STOP
                END IF

                READ (101, *)
                READ (101, *)
                READ (101, *) nmbrndsznd
                READ (101, *) nmbrelmntsznd
                READ (101, *) xrefnn, yrefnn, zrefnn

                CLOSE (101)

            END IF

            IF (MeshRead(ithprtl) == 1 .or. MeshRead(ithprtl) ==-1) THEN

                IF (     ithprtl < 10  ) THEN
                    WRITE(filename,'(I1)') ithprtl
                    filename = "000"//TRIM(filename)
                ELSE IF (ithprtl < 100 ) THEN
                    WRITE(filename,'(I2)') ithprtl
                    filename =  "00"//TRIM(filename)
                ELSE IF (ithprtl < 1000) THEN
                    WRITE(filename,'(I3)') ithprtl
                    filename =   "0"//TRIM(filename)
                ELSE
                    WRITE(filename,'(I4)') ithprtl
                END IF

                filename = "Prtl_"//TRIM(filename)//".msh"

                OPEN (101, FILE = filename, STATUS = "OLD", IOSTAT = IOS)
                IF (IOS /= 0) THEN
                    PRINT*, filename, ' does not exist! Please check!'
                    STOP
                END IF

                DO i = 1, 8
                    READ (101, *)
                END DO
                READ (101, *) nmbrndsznd
                DO i = 1, nmbrndsznd
                    READ (101, *)
                END DO
                READ (101, *)
                READ (101, *)
                READ (101, *) nmbrelmntsznd

                CLOSE (101)

            END IF

            IF (MeshRead(ithprtl) == 2 .or. MeshRead(ithprtl) ==-2) THEN

                IF (     ithprtl < 10  ) THEN
                    WRITE(filename,'(I1)') ithprtl
                    filename = "000"//TRIM(filename)
                ELSE IF (ithprtl < 100 ) THEN
                    WRITE(filename,'(I2)') ithprtl
                    filename =  "00"//TRIM(filename)
                ELSE IF (ithprtl < 1000) THEN
                    WRITE(filename,'(I3)') ithprtl
                    filename =   "0"//TRIM(filename)
                ELSE
                    WRITE(filename,'(I4)') ithprtl
                END IF

                filename = "Prtl_"//TRIM(filename)//".msh"

                OPEN (101, FILE = filename, STATUS = "OLD", IOSTAT = IOS)
                IF (IOS /= 0) THEN
                    PRINT*, filename, ' does not exist! Please check!'
                    STOP
                END IF

                DO i = 1, 8
                    READ (101, *)
                END DO
                READ (101, *) nmbrndsznd
                DO i = 1, nmbrndsznd
                    READ (101, *)
                END DO
                READ (101, *)
                READ (101, *)
                READ (101, *) nmbrelmntsznd

                nmbrelmntsznd = nmbrelmntsznd*2

                CLOSE (101)

            END IF

            ttlnmbrnd = ttlnmbrnd + nmbrndsznd
            ttlnmbrelmnt = ttlnmbrelmnt + nmbrelmntsznd

        END DO

        ALLOCATE(xnd(ttlnmbrnd))
        ALLOCATE(ynd(ttlnmbrnd))
        ALLOCATE(znd(ttlnmbrnd))
        ALLOCATE(nnx(ttlnmbrnd))
        ALLOCATE(nny(ttlnmbrnd))
        ALLOCATE(nnz(ttlnmbrnd))
        ALLOCATE(t1x(ttlnmbrnd))
        ALLOCATE(t1y(ttlnmbrnd))
        ALLOCATE(t1z(ttlnmbrnd))
        ALLOCATE(t2x(ttlnmbrnd))
        ALLOCATE(t2y(ttlnmbrnd))
        ALLOCATE(t2z(ttlnmbrnd))
        ALLOCATE(curvt1(ttlnmbrnd))
        ALLOCATE(curvt2(ttlnmbrnd))
        ALLOCATE(curvmn(ttlnmbrnd))
        ALLOCATE(curvt1th(ttlnmbrnd))
        ALLOCATE(curvt2th(ttlnmbrnd))
        ALLOCATE(curvmnth(ttlnmbrnd))

!$OMP PARALLEL PRIVATE (i)
!$OMP DO
        DO i = 1, ttlnmbrnd
            xnd(i) = 0.0d0
            ynd(i) = 0.0d0
            znd(i) = 0.0d0
            nnx(i) = 0.0d0
            nny(i) = 0.0d0
            nnz(i) = 0.0d0
            t1x(i) = 0.0d0
            t1y(i) = 0.0d0
            t1z(i) = 0.0d0
            t2x(i) = 0.0d0
            t2y(i) = 0.0d0
            t2z(i) = 0.0d0
            curvt1(i) = 0.0d0
            curvt2(i) = 0.0d0
            curvmn(i) = 0.0d0
            curvt1th(i) = 0.0d0
            curvt2th(i) = 0.0d0
            curvmnth(i) = 0.0d0
        END DO
!$OMP END DO
!$OMP END PARALLEL

        ALLOCATE(elmntarea(ttlnmbrelmnt))
        ALLOCATE(nnxelmnt (ttlnmbrelmnt))
        ALLOCATE(nnyelmnt (ttlnmbrelmnt))
        ALLOCATE(nnzelmnt (ttlnmbrelmnt))

!$OMP PARALLEL PRIVATE (i)
!$OMP DO
        DO i = 1, ttlnmbrelmnt
            elmntarea(i) = 0.0d0
            nnxelmnt(i)  = 0.0d0
            nnyelmnt(i)  = 0.0d0
            nnzelmnt(i)  = 0.0d0
        END DO
!$OMP END DO
!$OMP END PARALLEL

        ttlsrcfmm = ttlnmbrelmnt*n_glqtr2d
        ALLOCATE(srcfmm_vec(3,ttlsrcfmm))
        ALLOCATE(srcfmm_nrm(3,ttlsrcfmm))
        ALLOCATE(srcfmm_wght(ttlsrcfmm))
        IF (MeshType == "L") THEN
            ALLOCATE(elmntlnknd(ttlnmbrelmnt,3))
            ALLOCATE(srcfmm_wtnd(3,ttlsrcfmm))
        END IF
        IF (MeshType == "Q") THEN
            ALLOCATE(elmntlnknd(ttlnmbrelmnt,6))
            ALLOCATE(elmntlnkndlnr(4*ttlnmbrelmnt,3))
            elmntlnkndlnr(:,:) = 0
            ALLOCATE(srcfmm_wtnd(6,ttlsrcfmm))
        END IF
        elmntlnknd(:,:) = 0



        ndoffset1 = 0
        elmntoffset1 = 0
        DO ithprtl = 1, nmbrprtl

            IF (MeshRead(ithprtl) == 0) THEN

                IF (DABS(MeshRlvlStp(ithprtl)-1.0d0)<0.01d0) THEN
                    IF (MeshGnrtn(ithprtl) == "Icshdrl") &
                        & CALL GetMeshNdElmntIcshdrlPaper &
                            &(Meshnlvl(ithprtl),MeshRlvlStp(ithprtl),MeshType)
                    IF (MeshGnrtn(ithprtl) == "Icshsim") &
                        & CALL GetMeshNdElmntIcshdrlSimple &
                            &(Meshnlvl(ithprtl),MeshRlvlStp(ithprtl),MeshType)
                    IF (MeshGnrtn(ithprtl) == "TwoCrcl") &
                        & CALL GetMeshNdElmntTwoCrcl &
                            & (Meshnlvl(ithprtl),MeshRlvlStp(ithprtl),MeshType)
                ELSE
                    CALL GetMeshNdElmntTwoCrcl &
                        & (Meshnlvl(ithprtl),MeshRlvlStp(ithprtl),MeshType)
                END IF

                OPEN (101, FILE = "Prtl_Orgnl.inp", STATUS = "OLD", IOSTAT = IOS)
                IF (IOS /= 0) THEN
                    PRINT*, '"Prtl_Orgnl.inp" does not exist! Please check!'
                    STOP
                END IF

                READ (101, *)
                READ (101, *)
                READ (101, *) nmbrnd(ithprtl)
                READ (101, *) nmbrelmnt(ithprtl)
                READ (101, *) xrefnn, yrefnn, zrefnn

                CLOSE (101, STATUS = "DELETE")

                OPEN (111, FILE = "Prtl_Orgnl.vrt", STATUS = "OLD", IOSTAT = IOS)
                IF (IOS /= 0) THEN
                    PRINT*, "'Prtl_Orgnl.vrt' does not exist! Please check!"
                    STOP
                END IF

                DO i = 1, nmbrnd(ithprtl)
                    nd_i = i + ndoffset1
                    READ (111, *) ndindex, xnd(nd_i), ynd(nd_i), znd(nd_i)
                END DO

                CLOSE (111, STATUS = "DELETE")

                OPEN (121, FILE = "Prtl_Orgnl.cel", STATUS = "OLD", IOSTAT = IOS)
                IF (IOS /= 0) THEN
                    PRINT*, '"Prtl_Orgnl.cel" does not exist! Please check!'
                    STOP
                END IF

                DO k = 1, nmbrelmnt(ithprtl)

                    elmnt_k = k + elmntoffset1

                    IF (MeshType == "L") THEN

                        READ (121, *) elmntindex, ndA, ndB, ndC

                        ndA = ndA + ndoffset1
                        ndB = ndB + ndoffset1
                        ndC = ndC + ndoffset1

                        tpvctACx = xnd(ndA) - xnd(ndC)
                        tpvctACy = ynd(ndA) - ynd(ndC)
                        tpvctACz = znd(ndA) - znd(ndC)

                        tpvctBCx = xnd(ndB) - xnd(ndC)
                        tpvctBCy = ynd(ndB) - ynd(ndC)
                        tpvctBCz = znd(ndB) - znd(ndC)

                        tpvctRefCx = xrefnn - xnd(ndC)
                        tpvctRefCy = yrefnn - ynd(ndC)
                        tpvctRefCz = zrefnn - znd(ndC)

                        tpnrmlchck =tpvctRefCx * (tpvctACy*tpvctBCz - tpvctACz*tpvctBCy)&
                                & + tpvctRefCy * (tpvctACz*tpvctBCx - tpvctACx*tpvctBCz)&
                                & + tpvctRefCz * (tpvctACx*tpvctBCy - tpvctACy*tpvctBCx)

                        IF (NrmlInOut(ithprtl) == 1) THEN
                            IF (tpnrmlchck > 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndA
                                elmntlnknd(elmnt_k,2) = ndB
                                elmntlnknd(elmnt_k,3) = ndC
                            END IF

                            IF (tpnrmlchck < 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndB
                                elmntlnknd(elmnt_k,2) = ndA
                                elmntlnknd(elmnt_k,3) = ndC
                            END IF
                        END IF

                        IF (NrmlInOut(ithprtl) ==-1) THEN
                            IF (tpnrmlchck > 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndB
                                elmntlnknd(elmnt_k,2) = ndA
                                elmntlnknd(elmnt_k,3) = ndC
                            END IF

                            IF (tpnrmlchck < 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndA
                                elmntlnknd(elmnt_k,2) = ndB
                                elmntlnknd(elmnt_k,3) = ndC
                            END IF
                        END IF

                    END IF

                    IF (MeshType == "Q") THEN

                        READ (121, *) elmntindex, ndA, ndB, ndC, ndD, ndE, ndF

                        ndA = ndA + ndoffset1
                        ndB = ndB + ndoffset1
                        ndC = ndC + ndoffset1
                        ndD = ndD + ndoffset1
                        ndE = ndE + ndoffset1
                        ndF = ndF + ndoffset1

                        tpvctACx = xnd(ndA) - xnd(ndC)
                        tpvctACy = ynd(ndA) - ynd(ndC)
                        tpvctACz = znd(ndA) - znd(ndC)

                        tpvctBCx = xnd(ndB) - xnd(ndC)
                        tpvctBCy = ynd(ndB) - ynd(ndC)
                        tpvctBCz = znd(ndB) - znd(ndC)

                        tpvctRefCx = xrefnn - xnd(ndC)
                        tpvctRefCy = yrefnn - ynd(ndC)
                        tpvctRefCz = zrefnn - znd(ndC)

                        tpnrmlchck =tpvctRefCx * (tpvctACy*tpvctBCz - tpvctACz*tpvctBCy)&
                                & + tpvctRefCy * (tpvctACz*tpvctBCx - tpvctACx*tpvctBCz)&
                                & + tpvctRefCz * (tpvctACx*tpvctBCy - tpvctACy*tpvctBCx)

                        IF (NrmlInOut(ithprtl) == 1) THEN
                            IF (tpnrmlchck > 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndA
                                elmntlnknd(elmnt_k,2) = ndB
                                elmntlnknd(elmnt_k,3) = ndC
                                elmntlnknd(elmnt_k,4) = ndD
                                elmntlnknd(elmnt_k,5) = ndE
                                elmntlnknd(elmnt_k,6) = ndF
                            END IF

                            IF (tpnrmlchck < 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndB
                                elmntlnknd(elmnt_k,2) = ndA
                                elmntlnknd(elmnt_k,3) = ndC
                                elmntlnknd(elmnt_k,4) = ndD
                                elmntlnknd(elmnt_k,5) = ndF
                                elmntlnknd(elmnt_k,6) = ndE
                            END IF
                        END IF

                        IF (NrmlInOut(ithprtl) ==-1) THEN
                            IF (tpnrmlchck < 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndA
                                elmntlnknd(elmnt_k,2) = ndB
                                elmntlnknd(elmnt_k,3) = ndC
                                elmntlnknd(elmnt_k,4) = ndD
                                elmntlnknd(elmnt_k,5) = ndE
                                elmntlnknd(elmnt_k,6) = ndF
                            END IF

                            IF (tpnrmlchck > 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndB
                                elmntlnknd(elmnt_k,2) = ndA
                                elmntlnknd(elmnt_k,3) = ndC
                                elmntlnknd(elmnt_k,4) = ndD
                                elmntlnknd(elmnt_k,5) = ndF
                                elmntlnknd(elmnt_k,6) = ndE
                            END IF
                        END IF

                        elmntlnkndlnr(4*(elmnt_k-1)+1,1)=elmntlnknd(elmnt_k,1)
                        elmntlnkndlnr(4*(elmnt_k-1)+1,2)=elmntlnknd(elmnt_k,4)
                        elmntlnkndlnr(4*(elmnt_k-1)+1,3)=elmntlnknd(elmnt_k,6)

                        elmntlnkndlnr(4*(elmnt_k-1)+2,1)=elmntlnknd(elmnt_k,4)
                        elmntlnkndlnr(4*(elmnt_k-1)+2,2)=elmntlnknd(elmnt_k,2)
                        elmntlnkndlnr(4*(elmnt_k-1)+2,3)=elmntlnknd(elmnt_k,5)

                        elmntlnkndlnr(4*(elmnt_k-1)+3,1)=elmntlnknd(elmnt_k,5)
                        elmntlnkndlnr(4*(elmnt_k-1)+3,2)=elmntlnknd(elmnt_k,3)
                        elmntlnkndlnr(4*(elmnt_k-1)+3,3)=elmntlnknd(elmnt_k,6)

                        elmntlnkndlnr(4*(elmnt_k-1)+4,1)=elmntlnknd(elmnt_k,4)
                        elmntlnkndlnr(4*(elmnt_k-1)+4,2)=elmntlnknd(elmnt_k,5)
                        elmntlnkndlnr(4*(elmnt_k-1)+4,3)=elmntlnknd(elmnt_k,6)

                    END IF

                END DO

                CLOSE (121, STATUS = "DELETE")

            END IF


            IF (MeshRead(ithprtl) == 1 .or. MeshRead(ithprtl) ==-1) THEN

                IF (     ithprtl < 10  ) THEN
                    WRITE(filename,'(I1)') ithprtl
                    filename = "000"//TRIM(filename)
                ELSE IF (ithprtl < 100 ) THEN
                    WRITE(filename,'(I2)') ithprtl
                    filename =  "00"//TRIM(filename)
                ELSE IF (ithprtl < 1000) THEN
                    WRITE(filename,'(I3)') ithprtl
                    filename =   "0"//TRIM(filename)
                ELSE
                    WRITE(filename,'(I4)') ithprtl
                END IF

                filename = "Prtl_"//TRIM(filename)//".msh"

                OPEN (101, FILE = filename, STATUS = "OLD", IOSTAT = IOS)
                IF (IOS /= 0) THEN
                    PRINT*, filename, ' does not exist! Please check!'
                    STOP
                END IF

                DO i = 1, 8
                    READ (101, *)
                END DO
                READ (101, *) nmbrnd(ithprtl)
                DO i = 1, nmbrnd(ithprtl)
                    nd_i = i + ndoffset1
                    READ (101, *) ndindex, xnd(nd_i), ynd(nd_i), znd(nd_i)
                END DO
                READ (101, *)
                READ (101, *)
                READ (101, *) nmbrelmnt(ithprtl)

                xrefnn = xnrmref(ithprtl)
                yrefnn = ynrmref(ithprtl)
                zrefnn = znrmref(ithprtl)

                DO k = 1, nmbrelmnt(ithprtl)

                    elmnt_k = k + elmntoffset1

                    IF (MeshType == "L") THEN

                        READ (101, *) elmntindex, elmntindex, elmntindex, elmntindex, &
                                    & elmntindex, ndA, ndB, ndC

                        ndA = ndA + ndoffset1
                        ndB = ndB + ndoffset1
                        ndC = ndC + ndoffset1

                        tpvctACx = xnd(ndA) - xnd(ndC)
                        tpvctACy = ynd(ndA) - ynd(ndC)
                        tpvctACz = znd(ndA) - znd(ndC)

                        tpvctBCx = xnd(ndB) - xnd(ndC)
                        tpvctBCy = ynd(ndB) - ynd(ndC)
                        tpvctBCz = znd(ndB) - znd(ndC)

                        tpvctRefCx = xrefnn - xnd(ndC)
                        tpvctRefCy = yrefnn - ynd(ndC)
                        tpvctRefCz = zrefnn - znd(ndC)

                        tpnrmlchck =tpvctRefCx * (tpvctACy*tpvctBCz - tpvctACz*tpvctBCy)&
                                & + tpvctRefCy * (tpvctACz*tpvctBCx - tpvctACx*tpvctBCz)&
                                & + tpvctRefCz * (tpvctACx*tpvctBCy - tpvctACy*tpvctBCx)

                        IF (NrmlInOut(ithprtl) == 1) THEN
                            IF (tpnrmlchck > 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndA
                                elmntlnknd(elmnt_k,2) = ndB
                                elmntlnknd(elmnt_k,3) = ndC
                            END IF

                            IF (tpnrmlchck < 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndB
                                elmntlnknd(elmnt_k,2) = ndA
                                elmntlnknd(elmnt_k,3) = ndC
                            END IF
                        END IF

                        IF (NrmlInOut(ithprtl) ==-1) THEN
                            IF (tpnrmlchck > 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndB
                                elmntlnknd(elmnt_k,2) = ndA
                                elmntlnknd(elmnt_k,3) = ndC
                            END IF

                            IF (tpnrmlchck < 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndA
                                elmntlnknd(elmnt_k,2) = ndB
                                elmntlnknd(elmnt_k,3) = ndC
                            END IF
                        END IF

                    END IF

                    IF (MeshType == "Q") THEN

                        READ (101, *) elmntindex, elmntindex, elmntindex, elmntindex, &
                                    & elmntindex, ndA, ndB, ndC, ndD, ndE, ndF

                        ndA = ndA + ndoffset1
                        ndB = ndB + ndoffset1
                        ndC = ndC + ndoffset1
                        ndD = ndD + ndoffset1
                        ndE = ndE + ndoffset1
                        ndF = ndF + ndoffset1

                        tpvctACx = xnd(ndA) - xnd(ndC)
                        tpvctACy = ynd(ndA) - ynd(ndC)
                        tpvctACz = znd(ndA) - znd(ndC)

                        tpvctBCx = xnd(ndB) - xnd(ndC)
                        tpvctBCy = ynd(ndB) - ynd(ndC)
                        tpvctBCz = znd(ndB) - znd(ndC)

                        tpvctRefCx = xrefnn - xnd(ndC)
                        tpvctRefCy = yrefnn - ynd(ndC)
                        tpvctRefCz = zrefnn - znd(ndC)

                        tpnrmlchck =tpvctRefCx * (tpvctACy*tpvctBCz - tpvctACz*tpvctBCy)&
                                & + tpvctRefCy * (tpvctACz*tpvctBCx - tpvctACx*tpvctBCz)&
                                & + tpvctRefCz * (tpvctACx*tpvctBCy - tpvctACy*tpvctBCx)

                        IF (NrmlInOut(ithprtl) == 1) THEN
                            IF (tpnrmlchck > 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndA
                                elmntlnknd(elmnt_k,2) = ndB
                                elmntlnknd(elmnt_k,3) = ndC
                                elmntlnknd(elmnt_k,4) = ndD
                                elmntlnknd(elmnt_k,5) = ndE
                                elmntlnknd(elmnt_k,6) = ndF
                            END IF

                            IF (tpnrmlchck < 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndB
                                elmntlnknd(elmnt_k,2) = ndA
                                elmntlnknd(elmnt_k,3) = ndC
                                elmntlnknd(elmnt_k,4) = ndD
                                elmntlnknd(elmnt_k,5) = ndF
                                elmntlnknd(elmnt_k,6) = ndE
                            END IF
                        END IF

                        IF (NrmlInOut(ithprtl) ==-1) THEN
                            IF (tpnrmlchck < 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndA
                                elmntlnknd(elmnt_k,2) = ndB
                                elmntlnknd(elmnt_k,3) = ndC
                                elmntlnknd(elmnt_k,4) = ndD
                                elmntlnknd(elmnt_k,5) = ndE
                                elmntlnknd(elmnt_k,6) = ndF
                            END IF

                            IF (tpnrmlchck > 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndB
                                elmntlnknd(elmnt_k,2) = ndA
                                elmntlnknd(elmnt_k,3) = ndC
                                elmntlnknd(elmnt_k,4) = ndD
                                elmntlnknd(elmnt_k,5) = ndF
                                elmntlnknd(elmnt_k,6) = ndE
                            END IF
                        END IF

                        elmntlnkndlnr(4*(elmnt_k-1)+1,1)=elmntlnknd(elmnt_k,1)
                        elmntlnkndlnr(4*(elmnt_k-1)+1,2)=elmntlnknd(elmnt_k,4)
                        elmntlnkndlnr(4*(elmnt_k-1)+1,3)=elmntlnknd(elmnt_k,6)

                        elmntlnkndlnr(4*(elmnt_k-1)+2,1)=elmntlnknd(elmnt_k,4)
                        elmntlnkndlnr(4*(elmnt_k-1)+2,2)=elmntlnknd(elmnt_k,2)
                        elmntlnkndlnr(4*(elmnt_k-1)+2,3)=elmntlnknd(elmnt_k,5)

                        elmntlnkndlnr(4*(elmnt_k-1)+3,1)=elmntlnknd(elmnt_k,5)
                        elmntlnkndlnr(4*(elmnt_k-1)+3,2)=elmntlnknd(elmnt_k,3)
                        elmntlnkndlnr(4*(elmnt_k-1)+3,3)=elmntlnknd(elmnt_k,6)

                        elmntlnkndlnr(4*(elmnt_k-1)+4,1)=elmntlnknd(elmnt_k,4)
                        elmntlnkndlnr(4*(elmnt_k-1)+4,2)=elmntlnknd(elmnt_k,5)
                        elmntlnkndlnr(4*(elmnt_k-1)+4,3)=elmntlnknd(elmnt_k,6)

                    END IF

                END DO

                CLOSE (101)

            END IF

            IF (MeshRead(ithprtl) == 2 .or. MeshRead(ithprtl) ==-2) THEN

                IF (     ithprtl < 10  ) THEN
                    WRITE(filename,'(I1)') ithprtl
                    filename = "000"//TRIM(filename)
                ELSE IF (ithprtl < 100 ) THEN
                    WRITE(filename,'(I2)') ithprtl
                    filename =  "00"//TRIM(filename)
                ELSE IF (ithprtl < 1000) THEN
                    WRITE(filename,'(I3)') ithprtl
                    filename =   "0"//TRIM(filename)
                ELSE
                    WRITE(filename,'(I4)') ithprtl
                END IF

                filename = "Prtl_"//TRIM(filename)//".msh"

                OPEN (101, FILE = filename, STATUS = "OLD", IOSTAT = IOS)
                IF (IOS /= 0) THEN
                    PRINT*, filename, ' does not exist! Please check!'
                    STOP
                END IF

                DO i = 1, 8
                    READ (101, *)
                END DO
                READ (101, *) nmbrnd(ithprtl)
                DO i = 1, nmbrnd(ithprtl)
                    nd_i = i + ndoffset1
                    READ (101, *) ndindex, xnd(nd_i), ynd(nd_i), znd(nd_i)
                END DO
                READ (101, *)
                READ (101, *)
                READ (101, *) nmbrelmnt(ithprtl)

                nmbrelmnt(ithprtl) = 2*nmbrelmnt(ithprtl)

                xrefnn = xnrmref(ithprtl)
                yrefnn = ynrmref(ithprtl)
                zrefnn = znrmref(ithprtl)

                DO k = 1, nmbrelmnt(ithprtl)/2

                    READ (101, *) elmntindex, elmntindex, elmntindex, elmntindex, &
                    &             elmntindex, ndA, ndB, ndC, ndD

                    ndA = ndA + ndoffset1
                    ndB = ndB + ndoffset1
                    ndC = ndC + ndoffset1
                    ndD = ndD + ndoffset1

                    elmnt_k = 2*(k-1)+1 + elmntoffset1

                    IF (MeshType == "L") THEN

                        tpvctACx = xnd(ndA) - xnd(ndC)
                        tpvctACy = ynd(ndA) - ynd(ndC)
                        tpvctACz = znd(ndA) - znd(ndC)

                        tpvctBCx = xnd(ndB) - xnd(ndC)
                        tpvctBCy = ynd(ndB) - ynd(ndC)
                        tpvctBCz = znd(ndB) - znd(ndC)

                        tpvctRefCx = xrefnn - xnd(ndC)
                        tpvctRefCy = yrefnn - ynd(ndC)
                        tpvctRefCz = zrefnn - znd(ndC)

                        tpnrmlchck =tpvctRefCx * (tpvctACy*tpvctBCz - tpvctACz*tpvctBCy)&
                                & + tpvctRefCy * (tpvctACz*tpvctBCx - tpvctACx*tpvctBCz)&
                                & + tpvctRefCz * (tpvctACx*tpvctBCy - tpvctACy*tpvctBCx)

                        IF (NrmlInOut(ithprtl) == 1) THEN
                            IF (tpnrmlchck > 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndA
                                elmntlnknd(elmnt_k,2) = ndB
                                elmntlnknd(elmnt_k,3) = ndC
                            END IF

                            IF (tpnrmlchck < 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndB
                                elmntlnknd(elmnt_k,2) = ndA
                                elmntlnknd(elmnt_k,3) = ndC
                            END IF
                        END IF

                        IF (NrmlInOut(ithprtl) ==-1) THEN
                            IF (tpnrmlchck > 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndB
                                elmntlnknd(elmnt_k,2) = ndA
                                elmntlnknd(elmnt_k,3) = ndC
                            END IF

                            IF (tpnrmlchck < 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndA
                                elmntlnknd(elmnt_k,2) = ndB
                                elmntlnknd(elmnt_k,3) = ndC
                            END IF
                        END IF

                    END IF

                    ndE = ndA
                    ndA = ndC
                    ndB = ndD
                    ndC = ndE

                    elmnt_k = 2*(k-1)+2 + elmntoffset1

                    IF (MeshType == "L") THEN

                        tpvctACx = xnd(ndA) - xnd(ndC)
                        tpvctACy = ynd(ndA) - ynd(ndC)
                        tpvctACz = znd(ndA) - znd(ndC)

                        tpvctBCx = xnd(ndB) - xnd(ndC)
                        tpvctBCy = ynd(ndB) - ynd(ndC)
                        tpvctBCz = znd(ndB) - znd(ndC)

                        tpvctRefCx = xrefnn - xnd(ndC)
                        tpvctRefCy = yrefnn - ynd(ndC)
                        tpvctRefCz = zrefnn - znd(ndC)

                        tpnrmlchck =tpvctRefCx * (tpvctACy*tpvctBCz - tpvctACz*tpvctBCy)&
                                & + tpvctRefCy * (tpvctACz*tpvctBCx - tpvctACx*tpvctBCz)&
                                & + tpvctRefCz * (tpvctACx*tpvctBCy - tpvctACy*tpvctBCx)

                        IF (NrmlInOut(ithprtl) == 1) THEN
                            IF (tpnrmlchck > 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndA
                                elmntlnknd(elmnt_k,2) = ndB
                                elmntlnknd(elmnt_k,3) = ndC
                            END IF

                            IF (tpnrmlchck < 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndB
                                elmntlnknd(elmnt_k,2) = ndA
                                elmntlnknd(elmnt_k,3) = ndC
                            END IF
                        END IF

                        IF (NrmlInOut(ithprtl) ==-1) THEN
                            IF (tpnrmlchck > 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndB
                                elmntlnknd(elmnt_k,2) = ndA
                                elmntlnknd(elmnt_k,3) = ndC
                            END IF

                            IF (tpnrmlchck < 0.0d0) THEN
                                elmntlnknd(elmnt_k,1) = ndA
                                elmntlnknd(elmnt_k,2) = ndB
                                elmntlnknd(elmnt_k,3) = ndC
                            END IF
                        END IF

                    END IF

                END DO

                CLOSE (101)

            END IF


            ndoffset1 = ndoffset1 + nmbrnd(ithprtl)
            elmntoffset1 = elmntoffset1 + nmbrelmnt(ithprtl)

        END DO

        ndoffset1 = 0
        elmntoffset1 = 0
        DO ithprtl = 1, nmbrprtl

            ndstaID(ithprtl) = ndoffset1 + 1
            elstaID(ithprtl) = elmntoffset1 + 1

            ndoffset1 = ndoffset1 + nmbrnd(ithprtl)
            elmntoffset1 = elmntoffset1 + nmbrelmnt(ithprtl)

            ndendID(ithprtl) = ndoffset1
            elendID(ithprtl) = elmntoffset1

        END DO

! Generate or import the base surface geometry.

        DO ithprtl = 1, nmbrprtl

            IF (MeshRead(ithprtl) == 0) THEN


                IF (PrtlType(ithprtl) == "PrSp") THEN
!$OMP PARALLEL PRIVATE (i)
!$OMP DO
                    DO i = ndstaID(ithprtl), ndendID(ithprtl)

                        xnd(i) = xnd(i)*bvsa(ithprtl)
                        ynd(i) = ynd(i)*bvsa(ithprtl)
                    END DO
!$OMP END DO
!$OMP END PARALLEL
                END IF

                IF (PrtlType(ithprtl) == "ObSp") THEN
!$OMP PARALLEL PRIVATE (i)
!$OMP DO
                    DO i = ndstaID(ithprtl), ndendID(ithprtl)

                        znd(i) = znd(i)*bvsa(ithprtl)
                    END DO
!$OMP END DO
!$OMP END PARALLEL
                END IF

                IF (PrtlType(ithprtl) == "LgRd") THEN

                    IF (DABS(bvsa(ithprtl) - 1.0d0) < 1.E-3) THEN
                        cvsa(ithprtl) = 0.0d0
                    ELSE
                        CALL SearchLgRd_cvsa(ithprtl,bvsa(ithprtl),cvsa(ithprtl))
                    END IF

!$OMP PARALLEL PRIVATE (i,rdbscl,tpxyscl)
!$OMP DO
                    DO i = ndstaID(ithprtl), ndendID(ithprtl)

                        IF (DABS(bvsa(ithprtl) - 1.0d0) > 1.E-3) THEN
                            IF (DABS(znd(i)) < 1.E-7) THEN
                                rdbscl = bvsa(ithprtl)
                                tpxyscl = DSQRT(xnd(i)**2+ynd(i)**2)
                                xnd(i) = xnd(i) * rdbscl/tpxyscl
                                ynd(i) = ynd(i) * rdbscl/tpxyscl
                            ELSE IF (     DABS(znd(i) - 1.0d0) < 1.E-7 &
                                    &.OR. DABS(znd(i) + 1.0d0) < 1.E-7) THEN
                                rdbscl = 0.0d0
                            ELSE
                                CALL SearchLgRd_rdbscl(ithprtl, znd(i), &
                                                    &  cvsa(ithprtl), rdbscl)
                                tpxyscl = DSQRT(xnd(i)**2+ynd(i)**2)
                                xnd(i) = xnd(i) * rdbscl/tpxyscl
                                ynd(i) = ynd(i) * rdbscl/tpxyscl
                            END IF
                        END IF
                    END DO
!$OMP END DO
!$OMP END PARALLEL

                END IF

                IF (PrtlType(ithprtl) == "Disk") THEN

                    IF (DABS(bvsa(ithprtl) - 1.0d0) < 1.E-3) THEN
                        cvsa(ithprtl) = 0.0d0
                    ELSE
                        CALL SearchLgRd_cvsa(ithprtl,bvsa(ithprtl),cvsa(ithprtl))
                    END IF

!$OMP PARALLEL PRIVATE (i,rdbscl,tpxyscl)
!$OMP DO
                    DO i = ndstaID(ithprtl), ndendID(ithprtl)

                        tpxyscl = DSQRT(xnd(i)**2+ynd(i)**2)
                        IF (DABS(bvsa(ithprtl) - 1.0d0) > 1.E-3) THEN
                            IF (tpxyscl < 1.E-7) THEN
                                rdbscl = bvsa(ithprtl)
                                znd(i) = znd(i) * rdbscl/DABS(znd(i))
                            ELSE IF (     DABS(tpxyscl - 1.0d0) < 1.E-7 &
                                    &.OR. DABS(tpxyscl + 1.0d0) < 1.E-7) THEN
                                rdbscl = 0.0d0
                            ELSE
                                CALL SearchLgRd_rdbscl(ithprtl, tpxyscl, &
                                                    &  cvsa(ithprtl), rdbscl)
                                znd(i) = znd(i) * rdbscl/DABS(znd(i))
                            END IF
                        END IF
                    END DO
!$OMP END DO
!$OMP END PARALLEL

                END IF

                IF (PrtlType(ithprtl) == "DbBl") THEN

                    IF (DABS(dvsa(ithprtl) - 1.0d0) < 1.E-3) THEN
                        cvsa(ithprtl) = 0.0d0
                        bvsa(ithprtl) = 1.0d0
                    ELSE
                        CALL SearchDbBl_cvsa(ithprtl,dvsa(ithprtl),cvsa(ithprtl))
                        CALL SearchDbBl_bvsa(ithprtl,dvsa(ithprtl),cvsa(ithprtl),&
                                            &bvsa(ithprtl))
                    END IF

!$OMP PARALLEL PRIVATE (i,rdbscl,tpxyscl)
!$OMP DO
                    DO i = ndstaID(ithprtl), ndendID(ithprtl)

                        IF (DABS(dvsa(ithprtl) - 1.0d0) > 1.E-3) THEN
                            IF (DABS(znd(i)) < 1.E-7) THEN
                                rdbscl = dvsa(ithprtl)
                                tpxyscl = DSQRT(xnd(i)**2+ynd(i)**2)
                                xnd(i) = xnd(i) * rdbscl/tpxyscl
                                ynd(i) = ynd(i) * rdbscl/tpxyscl
                            ELSE IF (     DABS(znd(i) - 1.0d0) < 1.E-7 &
                                    &.OR. DABS(znd(i) + 1.0d0) < 1.E-7) THEN
                                rdbscl = 0.0d0
                            ELSE
                                CALL SearchDbBl_rdbscl(ithprtl,znd(i),&
                                                    &  dvsa(ithprtl),cvsa(ithprtl),rdbscl)
                                tpxyscl = DSQRT(xnd(i)**2+ynd(i)**2)
                                xnd(i) = xnd(i) * rdbscl/tpxyscl
                                ynd(i) = ynd(i) * rdbscl/tpxyscl
                            END IF
                        END IF
                    END DO
!$OMP END DO
!$OMP END PARALLEL

                END IF

                IF (PrtlType(ithprtl) == "Coin") THEN

                    IF (DABS(dvsa(ithprtl) - 1.0d0) < 1.E-3) THEN
                        cvsa(ithprtl) = 0.0d0
                        bvsa(ithprtl) = 1.0d0
                    ELSE
                        CALL SearchDbBl_cvsa(ithprtl,dvsa(ithprtl),cvsa(ithprtl))
                        CALL SearchDbBl_bvsa(ithprtl,dvsa(ithprtl),cvsa(ithprtl),&
                                            &bvsa(ithprtl))
                    END IF

!$OMP PARALLEL PRIVATE (i,rdbscl,tpxyscl)
!$OMP DO
                    DO i = ndstaID(ithprtl), ndendID(ithprtl)

                        IF (DABS(dvsa(ithprtl) - 1.0d0) > 1.E-3) THEN
                            tpxyscl = DSQRT(xnd(i)**2+ynd(i)**2)
                            IF (DABS(tpxyscl) < 1.E-7) THEN
                                rdbscl = dvsa(ithprtl)
                                znd(i) = znd(i)/DABS(znd(i))&
                                                &  *rdbscl
                            ELSE IF (     DABS(tpxyscl - 1.0d0) < 1.E-7 ) THEN
                                rdbscl = 0.0d0
                            ELSE
                                CALL SearchDbBl_rdbscl(ithprtl,tpxyscl,dvsa(ithprtl),&
                                                    &  cvsa(ithprtl),rdbscl)
                                znd(i) = znd(i)/DABS(znd(i))&
                                                &  *rdbscl
                            END IF
                        END IF
                    END DO
!$OMP END DO
!$OMP END PARALLEL

                END IF


                IF (PrtlType(ithprtl) == "Bowl") THEN

!$OMP PARALLEL PRIVATE (i,rdbscl)
!$OMP DO
                    DO i = ndstaID(ithprtl), ndendID(ithprtl)
                        znd(i) = -znd(i)
                    END DO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE (i,rdbscl)
!$OMP DO
                    DO i = ndstaID(ithprtl), ndendID(ithprtl)

                        rdbscl = DSQRT(xnd(i)**2+ynd(i)**2)
                        xnd(i) = 2.0d0*xnd(i)   !2.0*sin[theta]
                        ynd(i) = 2.0d0*ynd(i)   !2.0*sin[theta]
                        znd(i) = bowl_a(ithprtl)*znd(i) &
                        &   -bowl_b(ithprtl)*rdbscl*rdbscl   !a*cos[theta]-b*sin[theta]**2

                    END DO
!$OMP END DO
!$OMP END PARALLEL

                END IF

                IF (PrtlType(ithprtl) == "DfSp") THEN

!$OMP PARALLEL PRIVATE (i,rdbscl)
!$OMP DO
                    DO i = ndstaID(ithprtl), ndendID(ithprtl)

                        rdbscl = xnd(i)
                        xnd(i) =  SIGN(1.0d0,rdbscl)*dfsp_a(ithprtl)&
                        &             *(rdbscl**2)**(1.0d0/dfsp_l(ithprtl))
                        rdbscl = ynd(i)
                        ynd(i) =  SIGN(1.0d0,rdbscl)*dfsp_b(ithprtl)&
                        &             *(rdbscl**2)**(1.0d0/dfsp_m(ithprtl))
                        rdbscl = znd(i)
                        znd(i) =  SIGN(1.0d0,rdbscl)*dfsp_c(ithprtl)&
                        &             *(rdbscl**2)**(1.0d0/dfsp_n(ithprtl))

                    END DO
!$OMP END DO
!$OMP END PARALLEL

                END IF

                IF (PrtlType(ithprtl) == "HSph") THEN

!$OMP PARALLEL PRIVATE (i,rdbscl)
!$OMP DO
                    DO i = ndstaID(ithprtl), ndendID(ithprtl)

                        IF (znd(i) < 0.0d0) znd(i) = 0.0d0
                    END DO
!$OMP END DO
!$OMP END PARALLEL

                END IF

                IF (PrtlType(ithprtl) == "JoJo") THEN

!$OMP PARALLEL PRIVATE (i,rdbscl)
!$OMP DO
                    DO i = ndstaID(ithprtl), ndendID(ithprtl)

                        rdbscl = DACOS(znd(i))
                        xnd(i) = xnd(i)&
                        &*(1.0d0+bvsa(ithprtl)*DCOS(rdbscl)*DCOS(rdbscl))
                        ynd(i) = ynd(i)&
                        &*(1.0d0+bvsa(ithprtl)*DCOS(rdbscl)*DCOS(rdbscl))
                    END DO
!$OMP END DO
!$OMP END PARALLEL

                END IF

                IF (PrtlType(ithprtl) == "Pnty") THEN

!$OMP PARALLEL PRIVATE (i,rdbscl)
!$OMP DO
                    DO i = ndstaID(ithprtl), ndendID(ithprtl)

                        rdbscl = 0.5d0*DACOS(znd(i))
                        xnd(i) = xnd(i)&
                        &*(1.0d0+bvsa(ithprtl)*DCOS(rdbscl)*DCOS(rdbscl))
                        xnd(i) = 0.1d0*xnd(i)
                        ynd(i) = ynd(i)&
                        &*(1.0d0+bvsa(ithprtl)*DCOS(rdbscl)*DCOS(rdbscl))
                        ynd(i) = 0.1d0*ynd(i)
                    END DO
!$OMP END DO
!$OMP END PARALLEL

                END IF

                IF (PrtlType(ithprtl) == "Cheb" .OR. PrtlType(ithprtl) == "ChPr" .OR. PrtlType(ithprtl) == "ChOb") THEN

!$OMP PARALLEL PRIVATE (i,j,tp1,rdbscl)
!$OMP DO
                    DO i = ndstaID(ithprtl), ndendID(ithprtl)

                        rdbscl = DACOS(znd(i))  !theta
                        tp1 = 0.0d0
                        DO j = NINT(dfsp_a(ithprtl)), NINT(dfsp_b(ithprtl))
                            tp1 = tp1 + dfsp_c(ithprtl)*DCOS(rdbscl*DBLE(j))  !sum{eps*cos[n*theta]}
                        END DO
                        rdbscl = 1.0d0 + tp1 !R(theta) = 1 + sum{[(-1)^j]*eps*cos[n*theta]}
                        xnd(i) = xnd(i)*rdbscl  !new x
                        ynd(i) = ynd(i)*rdbscl  !new y
                        znd(i) = znd(i)*rdbscl  !new z

                        IF (PrtlType(ithprtl) == "ChOb") THEN
                            znd(i) = znd(i)*bvsa(ithprtl)
                        END IF

                        IF (PrtlType(ithprtl) == "ChPr") THEN
                            xnd(i) = xnd(i)*bvsa(ithprtl)
                            ynd(i) = ynd(i)*bvsa(ithprtl)
                        END IF

                    END DO
!$OMP END DO
!$OMP END PARALLEL

                END IF

            END IF

        END DO

! Apply object scale factors.

        DO ithprtl = 1, nmbrprtl

!$OMP PARALLEL PRIVATE (i)
!$OMP DO
            DO i = ndstaID(ithprtl), ndendID(ithprtl)

                xnd(i) = xnd(i)*sizezoom(ithprtl)
                ynd(i) = ynd(i)*sizezoom(ithprtl)
                znd(i) = znd(i)*sizezoom(ithprtl)

            END DO
!$OMP END DO
!$OMP END PARALLEL

! Apply object rotations.
!$OMP PARALLEL PRIVATE (i,tpx,tpy,tpz,tp1,tp2,tp3)
!$OMP DO
            DO i = ndstaID(ithprtl), ndendID(ithprtl)

                tpx = xnd(i)
                tpy = ynd(i)
                tpz = znd(i)

                tp1 = DCOS(anglecal_y(ithprtl))*DCOS(anglecal_z(ithprtl))
                tp2 =-DCOS(anglecal_y(ithprtl))*DSIN(anglecal_z(ithprtl))
                tp3 = DSIN(anglecal_y(ithprtl))
                xnd(i) = tpx*tp1+tpy*tp2+tpz*tp3

                tp1 = DSIN(anglecal_x(ithprtl))*DSIN(anglecal_y(ithprtl))&
                    &*DCOS(anglecal_z(ithprtl)) &
                    &+DCOS(anglecal_x(ithprtl))*DSIN(anglecal_z(ithprtl))
                tp2 =-DSIN(anglecal_x(ithprtl))*DSIN(anglecal_y(ithprtl))&
                    &*DSIN(anglecal_z(ithprtl)) &
                    &+DCOS(anglecal_x(ithprtl))*DCOS(anglecal_z(ithprtl))
                tp3 =-DSIN(anglecal_x(ithprtl))*DCOS(anglecal_y(ithprtl))
                ynd(i) = tpx*tp1+tpy*tp2+tpz*tp3

                tp1 =-DCOS(anglecal_x(ithprtl))*DSIN(anglecal_y(ithprtl))&
                    &*DCOS(anglecal_z(ithprtl)) &
                    &+DSIN(anglecal_x(ithprtl))*DSIN(anglecal_z(ithprtl))
                tp2 = DCOS(anglecal_x(ithprtl))*DSIN(anglecal_y(ithprtl))&
                    &*DSIN(anglecal_z(ithprtl)) &
                    &+DSIN(anglecal_x(ithprtl))*DCOS(anglecal_z(ithprtl))
                tp3 = DCOS(anglecal_x(ithprtl))*DCOS(anglecal_y(ithprtl))
                znd(i) = tpx*tp1+tpy*tp2+tpz*tp3

            END DO
!$OMP END DO
!$OMP END PARALLEL

! Shift bowl geometry to its local origin before the final translation.
            IF (PrtlType(ithprtl) == "Bowl") THEN
                i = ndstaID(ithprtl)
                tp1 = xnd(i)
                tp2 = ynd(i)
                tp3 = znd(i)
!$OMP PARALLEL PRIVATE (i,rdbscl)
!$OMP DO
                DO i = ndstaID(ithprtl), ndendID(ithprtl)
                    xnd(i) = xnd(i)-tp1
                    ynd(i) = ynd(i)-tp2
                    znd(i) = znd(i)-tp3
                END DO
!$OMP END DO
!$OMP END PARALLEL
            END IF

        ! Apply the user-specified translation.
!$OMP PARALLEL PRIVATE (i,tpx,tpy,tpz,tp1,tp2,tp3)
!$OMP DO
            DO i = ndstaID(ithprtl), ndendID(ithprtl)

                xnd(i) = xnd(i) + xloctn(ithprtl)
                ynd(i) = ynd(i) + yloctn(ithprtl)
                znd(i) = znd(i) + zloctn(ithprtl)

            END DO
!$OMP END DO
!$OMP END PARALLEL

        END DO

        CALL Getndlnkelmnt
        CALL Getndlnknd

        IF (MeshType == 'L') THEN
            ttlddtnd = mxnmbrndlnknd2ndslf
        END IF

        IF (MeshType == 'Q') THEN
            ttlddtnd = mxnmbrndlnknd2ndslf
        END IF

        ALLOCATE(d_dt1(ttlnmbrnd,ttlddtnd))
        ALLOCATE(d_dt2(ttlnmbrnd,ttlddtnd))
        ALLOCATE(d2_dt1(ttlnmbrnd,ttlddtnd))
        ALLOCATE(d2_dt2(ttlnmbrnd,ttlddtnd))

!$OMP PARALLEL PRIVATE (i,j)
!$OMP DO
        DO i = 1, ttlnmbrnd
            DO j = 1, ttlddtnd
                d_dt1(i,j) = 0.0d0
                d_dt2(i,j) = 0.0d0
                d2_dt1(i,j) = 0.0d0
                d2_dt2(i,j) = 0.0d0
            END DO
        END DO
!$OMP END DO
!$OMP END PARALLEL

    END SUBROUTINE MeshEditing

    ! Build the list of surface elements incident on each global node.

    SUBROUTINE Getndlnkelmnt

        INTEGER :: ithprtl, rec_kmx, rec_k, k, i, j

        mxnmbrndlnkelmnt = 0

        rec_kmx = 0

        ! Determine the maximum number of incident elements per node.
        DO k = 1, ttlnmbrnd

            rec_k = 0

            DO i = 1, ttlnmbrelmnt

                IF (MeshType == "L") THEN
                    DO j = 1, 3
                        IF (k == elmntlnknd(i,j)) THEN
                            rec_k = rec_k+1
                        END IF
                    END DO
                END IF
                IF (MeshType == "Q") THEN
                    DO j = 1, 6
                        IF (k == elmntlnknd(i,j)) THEN
                            rec_k = rec_k+1
                        END IF
                    END DO
                END IF
            END DO

            IF (rec_k>rec_kmx) rec_kmx = rec_k

        END DO

        ! Record the elements incident on each node.
        IF (mxnmbrndlnkelmnt < rec_kmx) mxnmbrndlnkelmnt = rec_kmx

        ALLOCATE (ndlnkelmnt(ttlnmbrnd, mxnmbrndlnkelmnt))
        ndlnkelmnt(:,:) = 0

        DO k = 1, ttlnmbrnd

            rec_k = 0

            DO i = 1, ttlnmbrelmnt

                IF (MeshType == "L") THEN
                    DO j = 1, 3
                        IF (k == elmntlnknd(i,j)) THEN
                            rec_k = rec_k+1
                            ndlnkelmnt(k, rec_k) = i
                        END IF
                    END DO
                END IF
                IF (MeshType == "Q") THEN
                    DO j = 1, 6
                        IF (k == elmntlnknd(i,j)) THEN
                            rec_k = rec_k+1
                            ndlnkelmnt(k, rec_k) = i
                        END IF
                    END DO
                END IF
            END DO

        END DO

        IF (MeshType == "Q") THEN

            mxnmbrndlnkelmntlnr = 0

            rec_kmx = 0

            ! Determine the maximum number of incident elements per node.
            DO k = 1, ttlnmbrnd

                rec_k = 0

                DO i = 1, 4*ttlnmbrelmnt

                    DO j = 1, 3
                        IF (k == elmntlnkndlnr(i,j)) THEN
                            rec_k = rec_k+1
                        END IF
                    END DO

                END DO

                IF (rec_k>rec_kmx) rec_kmx = rec_k

            END DO

            ! Record the elements incident on each node.
            IF (mxnmbrndlnkelmntlnr < rec_kmx) mxnmbrndlnkelmntlnr = rec_kmx


            ALLOCATE (ndlnkelmntlnr(ttlnmbrnd,mxnmbrndlnkelmntlnr))
            ndlnkelmntlnr(:,:) = 0

            DO k = 1, ttlnmbrnd

                rec_k = 0

                DO i = 1, 4*ttlnmbrelmnt

                    DO j = 1, 3
                        IF (k == elmntlnkndlnr(i,j)) THEN
                            rec_k = rec_k+1
                            ndlnkelmntlnr(k, rec_k) = i
                        END IF
                    END DO

                END DO

            END DO

        END IF


    END SUBROUTINE Getndlnkelmnt

    ! Build unique node-to-node adjacency from the element connectivity.

    SUBROUTINE Getndlnknd

        INTEGER :: rec_i, rec_imx, rec_imxthd, k, i, j, rdj, rdmxnmbrndlnknd, ithprtl
        INTEGER, ALLOCATABLE,  DIMENSION (:,:) :: rdndlnknd

        rdmxnmbrndlnknd = 100
        ALLOCATE (rdndlnknd(ttlnmbrnd, rdmxnmbrndlnknd))
        rdndlnknd(:,:) = 0

        mxnmbrndlnknd1st = 0

!Get the connected nodes of the first ring of the host node

        rec_imx = 0
        rec_imxthd = 0
        DO i = 1, ttlnmbrnd

! Gather first-ring neighbouring nodes from incident elements.
            rec_i = 0
            IF (MeshType == "L") THEN
                DO k = 1, mxnmbrndlnkelmnt
                    IF (ndlnkelmnt(i,k) /= 0) THEN
                        IF (elmntlnknd(ndlnkelmnt(i,k),1) == i) THEN
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnknd(ndlnkelmnt(i,k),2)
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnknd(ndlnkelmnt(i,k),3)
                        END IF
                        IF (elmntlnknd(ndlnkelmnt(i,k),2) == i) THEN
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnknd(ndlnkelmnt(i,k),3)
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnknd(ndlnkelmnt(i,k),1)
                        END IF
                        IF (elmntlnknd(ndlnkelmnt(i,k),3) == i) THEN
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnknd(ndlnkelmnt(i,k),1)
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnknd(ndlnkelmnt(i,k),2)
                        END IF
                    END IF
                END DO
            END IF

            IF (MeshType == "Q") THEN
                DO k = 1, mxnmbrndlnkelmntlnr
                    IF (ndlnkelmntlnr(i,k) /= 0) THEN
                        IF (elmntlnkndlnr(ndlnkelmntlnr(i,k),1) == i) THEN
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnkndlnr(ndlnkelmntlnr(i,k),2)
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnkndlnr(ndlnkelmntlnr(i,k),3)
                        END IF
                        IF (elmntlnkndlnr(ndlnkelmntlnr(i,k),2) == i) THEN
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnkndlnr(ndlnkelmntlnr(i,k),3)
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnkndlnr(ndlnkelmntlnr(i,k),1)
                        END IF
                        IF (elmntlnkndlnr(ndlnkelmntlnr(i,k),3) == i) THEN
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnkndlnr(ndlnkelmntlnr(i,k),1)
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnkndlnr(ndlnkelmntlnr(i,k),2)
                        END IF
                    END IF
                END DO
            END IF

! Remove duplicate neighbours.
            DO k = 1, rdmxnmbrndlnknd-1
                DO j = k+1, rdmxnmbrndlnknd
                    IF (rdndlnknd(i,k) == rdndlnknd(i,j)) &
                    & rdndlnknd(i,j) = 0
                END DO
            END DO

! Determine the maximum first-ring neighbour count.
            rec_i = 0
            DO k = 1, rdmxnmbrndlnknd
                IF (rdndlnknd(i,k) /= 0) THEN
                    rec_i = rec_i + 1
                END IF
            END DO
            rec_imxthd=MAX(rec_imxthd,rec_i)

        END DO

        rec_imx=MAX(rec_imx,rec_imxthd)

        mxnmbrndlnknd1st = MAX(mxnmbrndlnknd1st,rec_imx)

! Store first-ring adjacency in a compact array.
        ALLOCATE (ndlnknd1st(ttlnmbrnd, mxnmbrndlnknd1st))
        ndlnknd1st(:,:) = 0

        DO i = 1, ttlnmbrnd

            rec_i = 0

            DO k = 1, rdmxnmbrndlnknd
                IF (rdndlnknd(i,k) /= 0) THEN
                    rec_i = rec_i + 1
                    ndlnknd1st(i,rec_i) = rdndlnknd(i,k)
                END IF
            END DO
        END DO


        mxnmbrndlnknd1stslf = mxnmbrndlnknd1st+1
        ALLOCATE (ndlnknd1stslf(ttlnmbrnd, mxnmbrndlnknd1stslf))
        ndlnknd1stslf(:,:) = 0

        DO i = 1, ttlnmbrnd

            ndlnknd1stslf(i,1) = i
            DO k = 1, mxnmbrndlnknd1st
                ndlnknd1stslf(i,k+1) = ndlnknd1st(i,k)
            END DO
        END DO

!==================
! Build the combined first- and second-ring stencil.

        rdndlnknd(:,:) = 0

        mxnmbrndlnknd2nd = 0

        rec_imx = 0
        rec_imxthd = 0
        DO i = 1, ttlnmbrnd

! Seed the stencil with the first ring.
            rdj = 0
            DO j = 1, mxnmbrndlnknd1st
                IF (ndlnknd1st(i,j)/=0) THEN
                    rdj = rdj + 1
                    rdndlnknd(i,rdj) = ndlnknd1st(i,j)
                END IF
            END DO

! Add neighbours of the first-ring nodes.
            DO j = 1, mxnmbrndlnknd1st
                IF (ndlnknd1st(i,j)/=0) THEN
                    DO k = 1, mxnmbrndlnknd1st
                        IF (ndlnknd1st(ndlnknd1st(i,j),k)/=0) THEN
                            rdj = rdj + 1
                            rdndlnknd(i,rdj) &
                        & = ndlnknd1st(ndlnknd1st(i,j),k)
                        END IF
                    END DO
                END IF
            END DO

! Remove the host node from its own stencil.
            DO k = 1, rdmxnmbrndlnknd
                IF (rdndlnknd(i,k)==i) rdndlnknd(i,k) = 0
            END DO

! Remove duplicate neighbours.
            DO k = 1, rdmxnmbrndlnknd-1
                DO j = k+1, rdmxnmbrndlnknd
                    IF (rdndlnknd(i,j)==rdndlnknd(i,k)) &
                    & rdndlnknd(i,j) = 0
                END DO
            END DO

! Determine the maximum two-ring stencil size.
            rec_i = 0
            DO k = 1, rdmxnmbrndlnknd
                IF (rdndlnknd(i,k)/=0) rec_i = rec_i+1
            END DO
            rec_imxthd=MAX(rec_imxthd,rec_i)

        END DO

        rec_imx=MAX(rec_imx,rec_imxthd)

        mxnmbrndlnknd2nd = MAX(mxnmbrndlnknd2nd, rec_imx)


! Store the two-ring adjacency in a compact array.
        ALLOCATE (ndlnknd2nd(ttlnmbrnd, mxnmbrndlnknd2nd))
        ndlnknd2nd(:,:) = 0

        DO i = 1, ttlnmbrnd
            rec_i = 0
            DO k = 1, rdmxnmbrndlnknd
                IF (rdndlnknd(i,k)/=0) THEN
                    rec_i = rec_i+1
                    ndlnknd2nd(i,rec_i) = rdndlnknd(i,k)
                END IF
            END DO
        END DO


        mxnmbrndlnknd2ndslf = mxnmbrndlnknd2nd+1
        ALLOCATE (ndlnknd2ndslf(ttlnmbrnd, mxnmbrndlnknd2ndslf))
        ndlnknd2ndslf(:,:) = 0

        DO i = 1, ttlnmbrnd

            ndlnknd2ndslf(i,1) = i
            DO k = 1, mxnmbrndlnknd2nd
                ndlnknd2ndslf(i,k+1) = ndlnknd2nd(i,k)
            END DO
        END DO

        DEALLOCATE (rdndlnknd)

    END SUBROUTINE Getndlnknd

END MODULE
