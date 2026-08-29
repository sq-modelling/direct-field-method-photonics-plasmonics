
! SPDX-FileCopyrightText: 2026 Qiang Sun
! SPDX-License-Identifier: BSD-3-Clause

MODULE EM_DmnCal_SliceOutput

! Sample incident/source, scattered/correction, and total electromagnetic fields on
! Cartesian planes, and write the plane plus surface traces as Tecplot ASCII zones.
! Field names and the ex/in, E1/E2/E3, and phasor conventions are defined in
! EM_SurfCal_GlobalData.

    USE omp_lib

    USE Pre_Constants
    USE Pre_csvformat

    USE Geom_GlobalData

    USE BRIEFGHReal

    USE EM_SurfCal_GlobalData
    USE EM_SurfCal_PhysBC
    USE EM_DmnCal

    IMPLICIT NONE

    CONTAINS


    SUBROUTINE GetDmn_SliceOutput_EM(dmksmcs, dmxyz_RMpln)

! Purpose: classify a rectangular plane by material domain, evaluate E1/H1 and E2/H2
!          in the appropriate domain at every grid point, and append plotting zones.
! Input  : dmksmcs is the case/time-step integer used in zone labels.  dmxyz_RMpln
!          must be xy, yz, or zx; the caller expands the Input_Phys_EM.dat value 3D
!          into three separate calls.  Slice centre, full extents, and base grid count
!          are supplied through PostOfst*, PostEgSZ*, and PostEgnd_RMpln.
! Contract: geometry, normals, topology, materials, sources, and solved surface traces
!           must already be initialised.  Plane extents must be positive and the base
!           grid count must exceed one.  Lengths remain in the mesh length unit.
! Output : RsltPos.dat stores x,y,z,pos, where pos=0 is domain 0 and pos=p is the
!          material domain owned by surface p.  Rslt_Dmn1_IncomingField_plot.dat is
!          the domain-0 incident field over the complete plane.  Rslt_Dmn1_plot.dat,
!          Rslt_Dmn2_plot.dat, and Rslt_Dmn3_plot.dat store E1/H1, E2/H2, and their
!          total E3/H3, respectively, as surface and classified-plane zones.
! Columns: E/H components are REAL phasor parts and |E|/|H| are complex-vector norms.
!          S = 0.5*E x CONJG(H); Sx/Sy/Sz are REAL(S), while |S| is the norm of the
!          complex S vector.  No field or coordinate unit conversion is applied.

        INTEGER, INTENT (IN) :: dmksmcs
        CHARACTER(LEN=2), INTENT (IN) :: dmxyz_RMpln
        DOUBLE PRECISION :: dmEgSZ_RMpln
        INTEGER :: dmEgnd_RMpln,dmphsOn

        LOGICAL :: filexist, filexistIn, filexistOutput

        DOUBLE PRECISION :: tpendx1,tpendx2,tpendx3,tpendx4,tpendx5,tpendx6
        DOUBLE PRECISION :: tpendy1,tpendy2,tpendy3,tpendy4,tpendy5,tpendy6
        DOUBLE PRECISION :: tpendz1,tpendz2,tpendz3,tpendz4,tpendz5,tpendz6

        DOUBLE PRECISION :: pr0x, pr0y, pr0z, pr0nnx, pr0nny, pr0nnz, Dmnx, Dmny, Dmnz

        DOUBLE PRECISION :: tp, tp1, tp2, tp3, tp4, tp5, tp_phase,tpttlphs
        DOUBLE PRECISION :: tpPoyntscl

        INTEGER :: ithprtl, i, j, k, icmpnt
        INTEGER :: NdA, NdB, NdC, NdD, NdE, NdF, ndst, nded, elst, eled, id_tp

        DOUBLE PRECISION :: g1, g2, g3, g4, g5, g6
        DOUBLE PRECISION :: h1, h2, h3, h4, h5, h6
        DOUBLE PRECISION :: tpsum

        COMPLEX(KIND=KIND(1.0D0)) :: ztp, ztp1, ztp2, ztp3, ztp4, ztp5, ztp6, ztp7, ztp8, ztp9

        INTEGER :: nmbrnd_RMpln, nmbrelmnt_RMpln
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: x_RMpln, y_RMpln, z_RMpln
        INTEGER, ALLOCATABLE, DIMENSION (:) :: pos_RMpln
        INTEGER, ALLOCATABLE, DIMENSION (:,:) :: ndid_RMpln
        COMPLEX(KIND=KIND(1.0D0)), ALLOCATABLE, DIMENSION (:) :: E1xdmn,E1ydmn,E1zdmn, &
        &                                                        H1xdmn,H1ydmn,H1zdmn, &
        &                                                        E2xdmn,E2ydmn,E2zdmn, &
        &                                                        H2xdmn,H2ydmn,H2zdmn

! Build a centred uniform plane grid.  PostEgnd_RMpln sets the node count along
! x for xy, y for yz, and z for zx; the other count is scaled by the aspect ratio.
        IF (dmxyz_RMpln == 'xy') THEN
            tp = PostEgSZy_RMpln/PostEgSZx_RMpln
            dmEgnd_RMpln = CEILING(PostEgnd_RMpln * tp)
        END IF
        IF (dmxyz_RMpln == 'yz') THEN
            tp = PostEgSZz_RMpln/PostEgSZy_RMpln
            dmEgnd_RMpln = CEILING(PostEgnd_RMpln * tp)
        END IF
        IF (dmxyz_RMpln == 'zx') THEN
            tp = PostEgSZx_RMpln/PostEgSZz_RMpln
            dmEgnd_RMpln = CEILING(PostEgnd_RMpln * tp)
        END IF

        nmbrnd_RMpln = PostEgnd_RMpln*dmEgnd_RMpln
        nmbrelmnt_RMpln = 2*(PostEgnd_RMpln-1)*(dmEgnd_RMpln-1)
        ALLOCATE (x_RMpln(nmbrnd_RMpln))
        ALLOCATE (y_RMpln(nmbrnd_RMpln))
        ALLOCATE (z_RMpln(nmbrnd_RMpln))
        ALLOCATE (pos_RMpln(nmbrnd_RMpln))
        ALLOCATE (ndid_RMpln(PostEgnd_RMpln,dmEgnd_RMpln))

        IF (dmxyz_RMpln == 'xy') THEN
            tp1 = PostEgSZx_RMpln / DBLE(PostEgnd_RMpln-1)
            tp2 = PostEgSZy_RMpln / DBLE(dmEgnd_RMpln-1)
        END IF
        IF (dmxyz_RMpln == 'yz') THEN
            tp2 = PostEgSZy_RMpln / DBLE(PostEgnd_RMpln-1)
            tp3 = PostEgSZz_RMpln / DBLE(dmEgnd_RMpln-1)
        END IF
        IF (dmxyz_RMpln == 'zx') THEN
            tp3 = PostEgSZz_RMpln / DBLE(PostEgnd_RMpln-1)
            tp1 = PostEgSZx_RMpln / DBLE(dmEgnd_RMpln-1)
        END IF

!$OMP PARALLEL PRIVATE (i,j,k)
!$OMP DO
        DO i = 1, PostEgnd_RMpln
            DO j = 1, dmEgnd_RMpln
                k = (i-1)*dmEgnd_RMpln + j
                IF (dmxyz_RMpln == 'xy') THEN
                    x_RMpln(k) = PostOfstx_RMpln-0.5d0*PostEgSZx_RMpln + tp1*DBLE(i-1)
                    y_RMpln(k) = PostOfsty_RMpln-0.5d0*PostEgSZy_RMpln + tp2*DBLE(j-1)
                    z_RMpln(k) = PostOfstz_RMpln+0.0d0
                END IF
                IF (dmxyz_RMpln == 'yz') THEN
                    y_RMpln(k) = PostOfsty_RMpln-0.5d0*PostEgSZy_RMpln + tp2*DBLE(i-1)
                    z_RMpln(k) = PostOfstz_RMpln-0.5d0*PostEgSZz_RMpln + tp3*DBLE(j-1)
                    x_RMpln(k) = PostOfstx_RMpln+0.0d0
                END IF
                IF (dmxyz_RMpln == 'zx') THEN
                    z_RMpln(k) = PostOfstz_RMpln-0.5d0*PostEgSZz_RMpln + tp3*DBLE(i-1)
                    x_RMpln(k) = PostOfstx_RMpln-0.5d0*PostEgSZx_RMpln + tp1*DBLE(j-1)
                    y_RMpln(k) = PostOfsty_RMpln+0.0d0
                END IF
                ndid_RMpln(i,j) = k
            END DO
        END DO
!$OMP END DO
!$OMP END PARALLEL

! Classify each grid point with a solid-angle test.  The corelnkshell condition lets
! a nested domain replace its parent while preventing unrelated surfaces from doing so.
!$OMP PARALLEL PRIVATE (i,j,k,ithprtl,id_tp) &
!$OMP & PRIVATE (pr0x,pr0y,pr0z) &
!$OMP & PRIVATE (g1,g2,g3,g4,g5,g6,h1,h2,h3,h4,h5,h6,tp1,tp2,tpsum)
!$OMP DO
        DO i = 1, nmbrnd_RMpln

            pos_RMpln(i) = 0

            pr0x = x_RMpln(i)
            pr0y = y_RMpln(i)
            pr0z = z_RMpln(i)

            DO ithprtl = 1, nmbrprtl

                tpsum = 0.0d0

                DO k = elstaID(ithprtl), elendID(ithprtl)

                    IF (MeshType == "L") THEN

                        CALL CalGHLnrBRIEFLnrREAL( 0.0d0, k, &
                        &       pr0x, pr0y, pr0z, 0.0d0, 0.0d0, 0.0d0, &
                        &       g1, g2, g3, h1, h2, h3, tp1, tp2)

                        tpsum = tpsum + h1 + h2 + h3

                   END IF


                    IF (MeshType == "Q") THEN

                        CALL CalGHQdrBRIEFLnrREAL( 0.0d0, k, &
                        &       pr0x, pr0y, pr0z, 0.0d0, 0.0d0, 0.0d0, &
                        &       g1, g2, g3, g4, g5, g6, &
                        &       h1, h2, h3, h4, h5, h6, tp1, tp2)

                        tpsum = tpsum + h1 + h2 + h3 + h4 + h5 + h6

                    END IF

                END DO

                IF (DABS(tpsum) > 12.0d0) THEN
                    IF (pos_RMpln(i) == corelnkshell(ithprtl)) THEN
                        pos_RMpln(i) = ithprtl
                    END IF
                END IF

            END DO

        END DO
!$OMP END DO
!$OMP END PARALLEL

! Write the plane classification before evaluating any fields.
        INQUIRE (FILE="RsltPos.dat", exist=filexistOutput)
        IF (filexistOutput) THEN
            OPEN (111, FILE="RsltPos.dat", STATUS="OLD", &
                & POSITION="APPEND", ACTION="WRITE")
        ELSE
            OPEN (111, FILE="RsltPos.dat", STATUS="NEW", &
                & ACTION="WRITE")
            CALL csv_write_char(111,'Variables="x", "y", "z",',.false.,'spc')
            CALL csv_write_char(111,'"pos"',.true.,'spc')
        END IF

        CALL csv_write_char(111,'Zone T ="pos_',.false.,'spc')
        CALL csv_write_char(111,'plane - ',.false.,'spc')
        CALL csv_write_char(111,dmxyz_RMpln,.false.,'spc')
        CALL csv_write_char(111,'", n=',.false.,'spc')
        CALL csv_write_integer(111,nmbrnd_RMpln,.false.,'spc')
        CALL csv_write_char(111,', e=',.false.,'spc')
        CALL csv_write_integer(111,nmbrelmnt_RMpln,.false.,'spc')
        CALL csv_write_char(111,', f=fepoint, et=triangle',.true.,'spc')

        DO i = 1, nmbrnd_RMpln
            Dmnx = x_RMpln(i)
            Dmny = y_RMpln(i)
            Dmnz = z_RMpln(i)
            CALL csv_write_dble(111,Dmnx,.false.,'cmr')
            CALL csv_write_dble(111,Dmny,.false.,'cmr')
            CALL csv_write_dble(111,Dmnz,.false.,'cmr')

            CALL csv_write_integer(111,pos_RMpln(i),.true.,'spc')

        END DO

        DO i = 1, PostEgnd_RMpln - 1
            DO j = 1, dmEgnd_RMpln - 1
                CALL csv_write_integer(111,ndid_RMpln(i  ,j  ),.false.,'cmr')
                CALL csv_write_integer(111,ndid_RMpln(i+1,j  ),.false.,'cmr')
                CALL csv_write_integer(111,ndid_RMpln(i+1,j+1),.true.,'spc')

                CALL csv_write_integer(111,ndid_RMpln(i+1,j+1),.false.,'cmr')
                CALL csv_write_integer(111,ndid_RMpln(i  ,j+1),.false.,'cmr')
                CALL csv_write_integer(111,ndid_RMpln(i  ,j  ),.true.,'spc')
            END DO
        END DO

        CLOSE (111)

! Evaluate the two field parts once per plane point in the classified material domain.
        ALLOCATE (E1xdmn(nmbrnd_RMpln))
        ALLOCATE (E1ydmn(nmbrnd_RMpln))
        ALLOCATE (E1zdmn(nmbrnd_RMpln))
        ALLOCATE (H1xdmn(nmbrnd_RMpln))
        ALLOCATE (H1ydmn(nmbrnd_RMpln))
        ALLOCATE (H1zdmn(nmbrnd_RMpln))
        ALLOCATE (E2xdmn(nmbrnd_RMpln))
        ALLOCATE (E2ydmn(nmbrnd_RMpln))
        ALLOCATE (E2zdmn(nmbrnd_RMpln))
        ALLOCATE (H2xdmn(nmbrnd_RMpln))
        ALLOCATE (H2ydmn(nmbrnd_RMpln))
        ALLOCATE (H2zdmn(nmbrnd_RMpln))

!$OMP PARALLEL PRIVATE (i,j) &
!$OMP & PRIVATE (Dmnx,Dmny,Dmnz)
!$OMP DO
        DO i = 1, nmbrnd_RMpln

            Dmnx = x_RMpln(i)
            Dmny = y_RMpln(i)
            Dmnz = z_RMpln(i)

            CALL GetDmnE1H1_EM (Dmnx,Dmny,Dmnz,pos_RMpln(i), &
            &                   E1xdmn(i),E1ydmn(i),E1zdmn(i), &
            &                   H1xdmn(i),H1ydmn(i),H1zdmn(i) )

            CALL GetDmnE2H2_EM ('E',Dmnx,Dmny,Dmnz,pos_RMpln(i), &
            &                   E2xdmn(i),E2ydmn(i),E2zdmn(i) )
            CALL GetDmnE2H2_EM ('H',Dmnx,Dmny,Dmnz,pos_RMpln(i), &
            &                   H2xdmn(i),H2ydmn(i),H2zdmn(i) )

        END DO
!$OMP END DO
!$OMP END PARALLEL

! Domain-0 incident/source field across the complete plane, without material
! classification.  This is a reference field distinct from Rslt_Dmn1_plot.dat.

        INQUIRE (FILE="Rslt_Dmn1_IncomingField_plot.dat", exist=filexistOutput)
        IF (filexistOutput) THEN
            OPEN (111, FILE="Rslt_Dmn1_IncomingField_plot.dat", STATUS="OLD", &
                & POSITION="APPEND", ACTION="WRITE")
        ELSE
            OPEN (111, FILE="Rslt_Dmn1_IncomingField_plot.dat", STATUS="NEW", &
                & ACTION="WRITE")
            CALL csv_write_char(111,'Variables="x", "y", "z",',.false.,'spc')
            CALL csv_write_char(111,'"E1x","E1y","E1z","|E1|",',.false.,'spc')
            CALL csv_write_char(111,'"H1x","H1y","H1z","|H1|"',.false.,'spc')
            CALL csv_write_char(111,'"S1x","S1y","S1z","|S1|"',.true.,'spc')
        END IF

        CALL csv_write_char(111,'Zone T ="time step',.false.,'spc')
        CALL csv_write_integer(111,dmksmcs,.false.,'spc')
        CALL csv_write_char(111,'plane - ',.false.,'spc')
        CALL csv_write_char(111,dmxyz_RMpln,.false.,'spc')
        CALL csv_write_char(111,'", n=',.false.,'spc')
        CALL csv_write_integer(111,nmbrnd_RMpln,.false.,'spc')
        CALL csv_write_char(111,', e=',.false.,'spc')
        CALL csv_write_integer(111,nmbrelmnt_RMpln,.false.,'spc')
        CALL csv_write_char(111,', f=fepoint, et=triangle',.true.,'spc')

        DO i = 1, nmbrnd_RMpln
            Dmnx = x_RMpln(i)
            Dmny = y_RMpln(i)
            Dmnz = z_RMpln(i)
            CALL GetDmnE1H1_EM (Dmnx,Dmny,Dmnz,0, &
            &                   ztp1, ztp2, ztp3, ztp4, ztp5, ztp6 )
            CALL csv_write_dble(111,Dmnx,.false.,'cmr')
            CALL csv_write_dble(111,Dmny,.false.,'cmr')
            CALL csv_write_dble(111,Dmnz,.false.,'cmr')
            CALL csv_write_dble(111,REAL(ztp1),.false.,'cmr')
            CALL csv_write_dble(111,REAL(ztp2),.false.,'cmr')
            CALL csv_write_dble(111,REAL(ztp3),.false.,'cmr')
            tp = DSQRT( CDABS(ztp1)**2 &
            &          +CDABS(ztp2)**2 &
            &          +CDABS(ztp3)**2 )
            CALL csv_write_dble(111,tp,.false.,'cmr')
            CALL csv_write_dble(111,REAL(ztp4),.false.,'cmr')
            CALL csv_write_dble(111,REAL(ztp5),.false.,'cmr')
            CALL csv_write_dble(111,REAL(ztp6),.false.,'cmr')
            tp = DSQRT( CDABS(ztp4)**2 &
            &          +CDABS(ztp5)**2 &
            &          +CDABS(ztp6)**2 )
            CALL csv_write_dble(111,tp,.false.,'cmr')
            ztp7 = 0.5d0*(ztp2*DCONJG(ztp6) - ztp3*DCONJG(ztp5))
            ztp8 = 0.5d0*(ztp3*DCONJG(ztp4) - ztp1*DCONJG(ztp6))
            ztp9 = 0.5d0*(ztp1*DCONJG(ztp5) - ztp2*DCONJG(ztp4))
            tp = DSQRT( CDABS(ztp7)**2 &
            &          +CDABS(ztp8)**2 &
            &          +CDABS(ztp9)**2 )
            CALL csv_write_dble(111,REAL(ztp7),.false.,'cmr')
            CALL csv_write_dble(111,REAL(ztp8),.false.,'cmr')
            CALL csv_write_dble(111,REAL(ztp9),.false.,'cmr')
            CALL csv_write_dble(111,tp,.true.,'spc')
        END DO

        DO i = 1, PostEgnd_RMpln - 1
            DO j = 1, dmEgnd_RMpln - 1
                CALL csv_write_integer(111,ndid_RMpln(i  ,j  ),.false.,'cmr')
                CALL csv_write_integer(111,ndid_RMpln(i+1,j  ),.false.,'cmr')
                CALL csv_write_integer(111,ndid_RMpln(i+1,j+1),.true.,'spc')

                CALL csv_write_integer(111,ndid_RMpln(i+1,j+1),.false.,'cmr')
                CALL csv_write_integer(111,ndid_RMpln(i  ,j+1),.false.,'cmr')
                CALL csv_write_integer(111,ndid_RMpln(i  ,j  ),.true.,'spc')
            END DO
        END DO

        CLOSE (111)

! Incident/source field (E1,H1): ex and in surface zones followed by the plane zone.

        INQUIRE (FILE="Rslt_Dmn1_plot.dat", exist=filexistOutput)
        IF (filexistOutput) THEN
            OPEN (111, FILE="Rslt_Dmn1_plot.dat", STATUS="OLD", &
                & POSITION="APPEND", ACTION="WRITE")
        ELSE
            OPEN (111, FILE="Rslt_Dmn1_plot.dat", STATUS="NEW", &
                & ACTION="WRITE")
            CALL csv_write_char(111,'Variables="x", "y", "z",',.false.,'spc')
            CALL csv_write_char(111,'"E1x","E1y","E1z","|E1|",',.false.,'spc')
            CALL csv_write_char(111,'"H1x","H1y","H1z","|H1|"',.false.,'spc')
            CALL csv_write_char(111,'"S1x","S1y","S1z","|S1|"',.true.,'spc')
        END IF

        DO ithprtl = 1, nmbrprtl
            IF (MeshType == "L") THEN
                CALL csv_write_char(111,'Zone T ="time step',.false.,'spc')
                CALL csv_write_integer(111,dmksmcs,.false.,'spc')
                CALL csv_write_char(111,'Particle ',.false.,'spc')
                CALL csv_write_integer(111,ithprtl,.false.,'spc')
                CALL csv_write_char(111,' Out - ',.false.,'spc')
                CALL csv_write_char(111,dmxyz_RMpln,.false.,'spc')
                CALL csv_write_char(111,'", n=',.false.,'spc')
                CALL csv_write_integer(111,nmbrnd(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', e=',.false.,'spc')
                CALL csv_write_integer(111,nmbrelmnt(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', f=fepoint, et=triangle',.true.,'spc')
            END IF
            IF (MeshType == "Q") THEN
                CALL csv_write_char(111,'Zone T ="time step',.false.,'spc')
                CALL csv_write_integer(111,dmksmcs,.false.,'spc')
                CALL csv_write_char(111,'Particle ',.false.,'spc')
                CALL csv_write_integer(111,ithprtl,.false.,'spc')
                CALL csv_write_char(111,' Out - ',.false.,'spc')
                CALL csv_write_char(111,dmxyz_RMpln,.false.,'spc')
                CALL csv_write_char(111,'", n=',.false.,'spc')
                CALL csv_write_integer(111,nmbrnd(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', e=',.false.,'spc')
                CALL csv_write_integer(111,4*nmbrelmnt(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', f=fepoint, et=triangle',.true.,'spc')
            END IF

            DO i = ndstaID(ithprtl), ndendID(ithprtl)
                CALL csv_write_dble(111,xnd(i),.false.,'cmr')
                CALL csv_write_dble(111,ynd(i),.false.,'cmr')
                CALL csv_write_dble(111,znd(i),.false.,'cmr')
                CALL csv_write_dble(111,REAL(exE1x_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(exE1y_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(exE1z_EM(i)),.false.,'cmr')
                tp = DSQRT( CDABS(exE1x_EM(i))**2 &
                &          +CDABS(exE1y_EM(i))**2 &
                &          +CDABS(exE1z_EM(i))**2 )
                CALL csv_write_dble(111,tp,.false.,'cmr')
                CALL csv_write_dble(111,REAL(exH1x_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(exH1y_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(exH1z_EM(i)),.false.,'cmr')
                tp = DSQRT( CDABS(exH1x_EM(i))**2 &
                &          +CDABS(exH1y_EM(i))**2 &
                &          +CDABS(exH1z_EM(i))**2 )
                CALL csv_write_dble(111,tp,.false.,'cmr')
                ztp7 = 0.5d0*(exE1y_EM(i)*DCONJG(exH1z_EM(i)) - exE1z_EM(i)*DCONJG(exH1y_EM(i)))
                ztp8 = 0.5d0*(exE1z_EM(i)*DCONJG(exH1x_EM(i)) - exE1x_EM(i)*DCONJG(exH1z_EM(i)))
                ztp9 = 0.5d0*(exE1x_EM(i)*DCONJG(exH1y_EM(i)) - exE1y_EM(i)*DCONJG(exH1x_EM(i)))
                tp = DSQRT( CDABS(ztp7)**2 &
                &          +CDABS(ztp8)**2 &
                &          +CDABS(ztp9)**2 )
                CALL csv_write_dble(111,REAL(ztp7),.false.,'cmr')
                CALL csv_write_dble(111,REAL(ztp8),.false.,'cmr')
                CALL csv_write_dble(111,REAL(ztp9),.false.,'cmr')
                CALL csv_write_dble(111,tp,.true.,'spc')
            END DO
            icmpnt = ndstaID(ithprtl) - 1
            IF (MeshType == "L") THEN
                DO k = elstaID(ithprtl), elendID(ithprtl)
                    CALL csv_write_integer(111,elmntlnknd(k,1)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,2)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,3)-icmpnt,.true.,'spc')
                    END DO
            END IF
            IF (MeshType == "Q") THEN
                DO k = elstaID(ithprtl), elendID(ithprtl)
                    CALL csv_write_integer(111,elmntlnknd(k,1)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,2)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,3)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.true.,'spc')
                END DO
            END IF
        END DO

        CLOSE (111)

        INQUIRE (FILE="Rslt_Dmn1_plot.dat", exist=filexistOutput)
        IF (filexistOutput) THEN
            OPEN (111, FILE="Rslt_Dmn1_plot.dat", STATUS="OLD", &
                & POSITION="APPEND", ACTION="WRITE")
        ELSE
            OPEN (111, FILE="Rslt_Dmn1_plot.dat", STATUS="NEW", &
                & ACTION="WRITE")
            CALL csv_write_char(111,'Variables="x", "y", "z",',.false.,'spc')
            CALL csv_write_char(111,'"E1x","E1y","E1z","|E1|",',.false.,'spc')
            CALL csv_write_char(111,'"H1x","H1y","H1z","|H1|"',.false.,'spc')
            CALL csv_write_char(111,'"S1x","S1y","S1z","|S1|"',.true.,'spc')
        END IF

        DO ithprtl = 1, nmbrprtl
            IF (MeshType == "L") THEN
                CALL csv_write_char(111,'Zone T ="time step',.false.,'spc')
                CALL csv_write_integer(111,dmksmcs,.false.,'spc')
                CALL csv_write_char(111,'Particle ',.false.,'spc')
                CALL csv_write_integer(111,ithprtl,.false.,'spc')
                CALL csv_write_char(111,' In - ',.false.,'spc')
                CALL csv_write_char(111,dmxyz_RMpln,.false.,'spc')
                CALL csv_write_char(111,'", n=',.false.,'spc')
                CALL csv_write_integer(111,nmbrnd(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', e=',.false.,'spc')
                CALL csv_write_integer(111,nmbrelmnt(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', f=fepoint, et=triangle',.true.,'spc')
            END IF
            IF (MeshType == "Q") THEN
                CALL csv_write_char(111,'Zone T ="time step',.false.,'spc')
                CALL csv_write_integer(111,dmksmcs,.false.,'spc')
                CALL csv_write_char(111,'Particle ',.false.,'spc')
                CALL csv_write_integer(111,ithprtl,.false.,'spc')
                CALL csv_write_char(111,' In - ',.false.,'spc')
                CALL csv_write_char(111,dmxyz_RMpln,.false.,'spc')
                CALL csv_write_char(111,'", n=',.false.,'spc')
                CALL csv_write_integer(111,nmbrnd(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', e=',.false.,'spc')
                CALL csv_write_integer(111,4*nmbrelmnt(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', f=fepoint, et=triangle',.true.,'spc')
            END IF
            DO i = ndstaID(ithprtl), ndendID(ithprtl)
                CALL csv_write_dble(111,xnd(i),.false.,'cmr')
                CALL csv_write_dble(111,ynd(i),.false.,'cmr')
                CALL csv_write_dble(111,znd(i),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inE1x_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inE1y_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inE1z_EM(i)),.false.,'cmr')
                tp = DSQRT( CDABS(inE1x_EM(i))**2 &
                &          +CDABS(inE1y_EM(i))**2 &
                &          +CDABS(inE1z_EM(i))**2 )
                CALL csv_write_dble(111,tp,.false.,'cmr')
                CALL csv_write_dble(111,REAL(inH1x_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inH1y_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inH1z_EM(i)),.false.,'cmr')
                tp = DSQRT( CDABS(inH1x_EM(i))**2 &
                &          +CDABS(inH1y_EM(i))**2 &
                &          +CDABS(inH1z_EM(i))**2 )
                CALL csv_write_dble(111,tp,.false.,'cmr')
                ztp7 = 0.5d0*(inE1y_EM(i)*DCONJG(inH1z_EM(i)) - inE1z_EM(i)*DCONJG(inH1y_EM(i)))
                ztp8 = 0.5d0*(inE1z_EM(i)*DCONJG(inH1x_EM(i)) - inE1x_EM(i)*DCONJG(inH1z_EM(i)))
                ztp9 = 0.5d0*(inE1x_EM(i)*DCONJG(inH1y_EM(i)) - inE1y_EM(i)*DCONJG(inH1x_EM(i)))
                tp = DSQRT( CDABS(ztp7)**2 &
                &          +CDABS(ztp8)**2 &
                &          +CDABS(ztp9)**2 )
                CALL csv_write_dble(111,REAL(ztp7),.false.,'cmr')
                CALL csv_write_dble(111,REAL(ztp8),.false.,'cmr')
                CALL csv_write_dble(111,REAL(ztp9),.false.,'cmr')
                CALL csv_write_dble(111,tp,.true.,'spc')
            END DO
            icmpnt = ndstaID(ithprtl) - 1
            IF (MeshType == "L") THEN
                DO k = elstaID(ithprtl), elendID(ithprtl)
                    CALL csv_write_integer(111,elmntlnknd(k,1)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,2)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,3)-icmpnt,.true.,'spc')
                    END DO
            END IF
            IF (MeshType == "Q") THEN
                DO k = elstaID(ithprtl), elendID(ithprtl)
                    CALL csv_write_integer(111,elmntlnknd(k,1)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,2)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,3)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.true.,'spc')
                END DO
            END IF
        END DO

        CLOSE (111)

        INQUIRE (FILE="Rslt_Dmn1_plot.dat", exist=filexistOutput)
        IF (filexistOutput) THEN
            OPEN (111, FILE="Rslt_Dmn1_plot.dat", STATUS="OLD", &
                & POSITION="APPEND", ACTION="WRITE")
        ELSE
            OPEN (111, FILE="Rslt_Dmn1_plot.dat", STATUS="NEW", &
                & ACTION="WRITE")
            CALL csv_write_char(111,'Variables="x", "y", "z",',.false.,'spc')
            CALL csv_write_char(111,'"E1x","E1y","E1z","|E1|",',.false.,'spc')
            CALL csv_write_char(111,'"H1x","H1y","H1z","|H1|"',.false.,'spc')
            CALL csv_write_char(111,'"S1x","S1y","S1z","|S1|"',.true.,'spc')
        END IF

        CALL csv_write_char(111,'Zone T ="time step',.false.,'spc')
        CALL csv_write_integer(111,dmksmcs,.false.,'spc')
        CALL csv_write_char(111,'plane - ',.false.,'spc')
        CALL csv_write_char(111,dmxyz_RMpln,.false.,'spc')
        CALL csv_write_char(111,'", n=',.false.,'spc')
        CALL csv_write_integer(111,nmbrnd_RMpln,.false.,'spc')
        CALL csv_write_char(111,', e=',.false.,'spc')
        CALL csv_write_integer(111,nmbrelmnt_RMpln,.false.,'spc')
        CALL csv_write_char(111,', f=fepoint, et=triangle',.true.,'spc')

        DO i = 1, nmbrnd_RMpln
            Dmnx = x_RMpln(i)
            Dmny = y_RMpln(i)
            Dmnz = z_RMpln(i)
            CALL csv_write_dble(111,Dmnx,.false.,'cmr')
            CALL csv_write_dble(111,Dmny,.false.,'cmr')
            CALL csv_write_dble(111,Dmnz,.false.,'cmr')
            CALL csv_write_dble(111,REAL(E1xdmn(i)),.false.,'cmr')
            CALL csv_write_dble(111,REAL(E1ydmn(i)),.false.,'cmr')
            CALL csv_write_dble(111,REAL(E1zdmn(i)),.false.,'cmr')
            tp = DSQRT( CDABS(E1xdmn(i))**2 &
            &          +CDABS(E1ydmn(i))**2 &
            &          +CDABS(E1zdmn(i))**2 )
            CALL csv_write_dble(111,tp,.false.,'cmr')
            CALL csv_write_dble(111,REAL(H1xdmn(i)),.false.,'cmr')
            CALL csv_write_dble(111,REAL(H1ydmn(i)),.false.,'cmr')
            CALL csv_write_dble(111,REAL(H1zdmn(i)),.false.,'cmr')
            tp = DSQRT( CDABS(H1xdmn(i))**2 &
            &          +CDABS(H1ydmn(i))**2 &
            &          +CDABS(H1zdmn(i))**2 )
            CALL csv_write_dble(111,tp,.false.,'cmr')
            ztp7 = 0.5d0*(E1ydmn(i)*DCONJG(H1zdmn(i)) - E1zdmn(i)*DCONJG(H1ydmn(i)))
            ztp8 = 0.5d0*(E1zdmn(i)*DCONJG(H1xdmn(i)) - E1xdmn(i)*DCONJG(H1zdmn(i)))
            ztp9 = 0.5d0*(E1xdmn(i)*DCONJG(H1ydmn(i)) - E1ydmn(i)*DCONJG(H1xdmn(i)))
            tp = DSQRT( CDABS(ztp7)**2 &
            &          +CDABS(ztp8)**2 &
            &          +CDABS(ztp9)**2 )
            CALL csv_write_dble(111,REAL(ztp7),.false.,'cmr')
            CALL csv_write_dble(111,REAL(ztp8),.false.,'cmr')
            CALL csv_write_dble(111,REAL(ztp9),.false.,'cmr')
            CALL csv_write_dble(111,tp,.true.,'spc')
        END DO

        DO i = 1, PostEgnd_RMpln - 1
            DO j = 1, dmEgnd_RMpln - 1
                CALL csv_write_integer(111,ndid_RMpln(i  ,j  ),.false.,'cmr')
                CALL csv_write_integer(111,ndid_RMpln(i+1,j  ),.false.,'cmr')
                CALL csv_write_integer(111,ndid_RMpln(i+1,j+1),.true.,'spc')

                CALL csv_write_integer(111,ndid_RMpln(i+1,j+1),.false.,'cmr')
                CALL csv_write_integer(111,ndid_RMpln(i  ,j+1),.false.,'cmr')
                CALL csv_write_integer(111,ndid_RMpln(i  ,j  ),.true.,'spc')
            END DO
        END DO

        CLOSE (111)

! Scattered/correction field (E2,H2): ex and in surface zones followed by the plane zone.

        INQUIRE (FILE="Rslt_Dmn2_plot.dat", exist=filexistOutput)
        IF (filexistOutput) THEN
            OPEN (111, FILE="Rslt_Dmn2_plot.dat", STATUS="OLD", &
                & POSITION="APPEND", ACTION="WRITE")
        ELSE
            OPEN (111, FILE="Rslt_Dmn2_plot.dat", STATUS="NEW", &
                & ACTION="WRITE")
            CALL csv_write_char(111,'Variables="x", "y", "z",',.false.,'spc')
            CALL csv_write_char(111,'"E2x","E2y","E2z","|E2|",',.false.,'spc')
            CALL csv_write_char(111,'"H2x","H2y","H2z","|H2|"',.false.,'spc')
            CALL csv_write_char(111,'"S2x","S2y","S2z","|S2|"',.true.,'spc')
        END IF

        DO ithprtl = 1, nmbrprtl
            IF (MeshType == "L") THEN
                CALL csv_write_char(111,'Zone T ="time step',.false.,'spc')
                CALL csv_write_integer(111,dmksmcs,.false.,'spc')
                CALL csv_write_char(111,'Particle ',.false.,'spc')
                CALL csv_write_integer(111,ithprtl,.false.,'spc')
                CALL csv_write_char(111,' Out - ',.false.,'spc')
                CALL csv_write_char(111,dmxyz_RMpln,.false.,'spc')
                CALL csv_write_char(111,'", n=',.false.,'spc')
                CALL csv_write_integer(111,nmbrnd(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', e=',.false.,'spc')
                CALL csv_write_integer(111,nmbrelmnt(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', f=fepoint, et=triangle',.true.,'spc')
            END IF
            IF (MeshType == "Q") THEN
                CALL csv_write_char(111,'Zone T ="time step',.false.,'spc')
                CALL csv_write_integer(111,dmksmcs,.false.,'spc')
                CALL csv_write_char(111,'Particle ',.false.,'spc')
                CALL csv_write_integer(111,ithprtl,.false.,'spc')
                CALL csv_write_char(111,' Out - ',.false.,'spc')
                CALL csv_write_char(111,dmxyz_RMpln,.false.,'spc')
                CALL csv_write_char(111,'", n=',.false.,'spc')
                CALL csv_write_integer(111,nmbrnd(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', e=',.false.,'spc')
                CALL csv_write_integer(111,4*nmbrelmnt(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', f=fepoint, et=triangle',.true.,'spc')
            END IF

            DO i = ndstaID(ithprtl), ndendID(ithprtl)
                CALL csv_write_dble(111,xnd(i),.false.,'cmr')
                CALL csv_write_dble(111,ynd(i),.false.,'cmr')
                CALL csv_write_dble(111,znd(i),.false.,'cmr')
                CALL csv_write_dble(111,REAL(exE2x_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(exE2y_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(exE2z_EM(i)),.false.,'cmr')
                tp = DSQRT( CDABS(exE2x_EM(i))**2 &
                &          +CDABS(exE2y_EM(i))**2 &
                &          +CDABS(exE2z_EM(i))**2 )
                CALL csv_write_dble(111,tp,.false.,'cmr')
                CALL csv_write_dble(111,REAL(exH2x_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(exH2y_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(exH2z_EM(i)),.false.,'cmr')
                tp = DSQRT( CDABS(exH2x_EM(i))**2 &
                &          +CDABS(exH2y_EM(i))**2 &
                &          +CDABS(exH2z_EM(i))**2 )
                CALL csv_write_dble(111,tp,.false.,'cmr')
                ztp7 = 0.5d0*(exE2y_EM(i)*DCONJG(exH2z_EM(i)) - exE2z_EM(i)*DCONJG(exH2y_EM(i)))
                ztp8 = 0.5d0*(exE2z_EM(i)*DCONJG(exH2x_EM(i)) - exE2x_EM(i)*DCONJG(exH2z_EM(i)))
                ztp9 = 0.5d0*(exE2x_EM(i)*DCONJG(exH2y_EM(i)) - exE2y_EM(i)*DCONJG(exH2x_EM(i)))
                tp = DSQRT( CDABS(ztp7)**2 &
                &          +CDABS(ztp8)**2 &
                &          +CDABS(ztp9)**2 )
                CALL csv_write_dble(111,REAL(ztp7),.false.,'cmr')
                CALL csv_write_dble(111,REAL(ztp8),.false.,'cmr')
                CALL csv_write_dble(111,REAL(ztp9),.false.,'cmr')
                CALL csv_write_dble(111,tp,.true.,'spc')
            END DO
            icmpnt = ndstaID(ithprtl) - 1
            IF (MeshType == "L") THEN
                DO k = elstaID(ithprtl), elendID(ithprtl)
                    CALL csv_write_integer(111,elmntlnknd(k,1)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,2)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,3)-icmpnt,.true.,'spc')
                    END DO
            END IF
            IF (MeshType == "Q") THEN
                DO k = elstaID(ithprtl), elendID(ithprtl)
                    CALL csv_write_integer(111,elmntlnknd(k,1)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,2)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,3)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.true.,'spc')
                END DO
            END IF
        END DO

        CLOSE (111)

        INQUIRE (FILE="Rslt_Dmn2_plot.dat", exist=filexistOutput)
        IF (filexistOutput) THEN
            OPEN (111, FILE="Rslt_Dmn2_plot.dat", STATUS="OLD", &
                & POSITION="APPEND", ACTION="WRITE")
        ELSE
            OPEN (111, FILE="Rslt_Dmn2_plot.dat", STATUS="NEW", &
                & ACTION="WRITE")
            CALL csv_write_char(111,'Variables="x", "y", "z",',.false.,'spc')
            CALL csv_write_char(111,'"E2x","E2y","E2z","|E2|",',.false.,'spc')
            CALL csv_write_char(111,'"H2x","H2y","H2z","|H2|"',.false.,'spc')
            CALL csv_write_char(111,'"S2x","S2y","S2z","|S2|"',.true.,'spc')
        END IF

        DO ithprtl = 1, nmbrprtl
            IF (MeshType == "L") THEN
                CALL csv_write_char(111,'Zone T ="time step',.false.,'spc')
                CALL csv_write_integer(111,dmksmcs,.false.,'spc')
                CALL csv_write_char(111,'Particle ',.false.,'spc')
                CALL csv_write_integer(111,ithprtl,.false.,'spc')
                CALL csv_write_char(111,' In - ',.false.,'spc')
                CALL csv_write_char(111,dmxyz_RMpln,.false.,'spc')
                CALL csv_write_char(111,'", n=',.false.,'spc')
                CALL csv_write_integer(111,nmbrnd(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', e=',.false.,'spc')
                CALL csv_write_integer(111,nmbrelmnt(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', f=fepoint, et=triangle',.true.,'spc')
            END IF
            IF (MeshType == "Q") THEN
                CALL csv_write_char(111,'Zone T ="time step',.false.,'spc')
                CALL csv_write_integer(111,dmksmcs,.false.,'spc')
                CALL csv_write_char(111,'Particle ',.false.,'spc')
                CALL csv_write_integer(111,ithprtl,.false.,'spc')
                CALL csv_write_char(111,' In - ',.false.,'spc')
                CALL csv_write_char(111,dmxyz_RMpln,.false.,'spc')
                CALL csv_write_char(111,'", n=',.false.,'spc')
                CALL csv_write_integer(111,nmbrnd(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', e=',.false.,'spc')
                CALL csv_write_integer(111,4*nmbrelmnt(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', f=fepoint, et=triangle',.true.,'spc')
            END IF
            DO i = ndstaID(ithprtl), ndendID(ithprtl)
                CALL csv_write_dble(111,xnd(i),.false.,'cmr')
                CALL csv_write_dble(111,ynd(i),.false.,'cmr')
                CALL csv_write_dble(111,znd(i),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inE2x_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inE2y_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inE2z_EM(i)),.false.,'cmr')
                tp = DSQRT( CDABS(inE2x_EM(i))**2 &
                &          +CDABS(inE2y_EM(i))**2 &
                &          +CDABS(inE2z_EM(i))**2 )
                CALL csv_write_dble(111,tp,.false.,'cmr')
                CALL csv_write_dble(111,REAL(inH2x_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inH2y_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inH2z_EM(i)),.false.,'cmr')
                tp = DSQRT( CDABS(inH2x_EM(i))**2 &
                &          +CDABS(inH2y_EM(i))**2 &
                &          +CDABS(inH2z_EM(i))**2 )
                CALL csv_write_dble(111,tp,.false.,'cmr')
                ztp7 = 0.5d0*(inE2y_EM(i)*DCONJG(inH2z_EM(i)) - inE2z_EM(i)*DCONJG(inH2y_EM(i)))
                ztp8 = 0.5d0*(inE2z_EM(i)*DCONJG(inH2x_EM(i)) - inE2x_EM(i)*DCONJG(inH2z_EM(i)))
                ztp9 = 0.5d0*(inE2x_EM(i)*DCONJG(inH2y_EM(i)) - inE2y_EM(i)*DCONJG(inH2x_EM(i)))
                tp = DSQRT( CDABS(ztp7)**2 &
                &          +CDABS(ztp8)**2 &
                &          +CDABS(ztp9)**2 )
                CALL csv_write_dble(111,REAL(ztp7),.false.,'cmr')
                CALL csv_write_dble(111,REAL(ztp8),.false.,'cmr')
                CALL csv_write_dble(111,REAL(ztp9),.false.,'cmr')
                CALL csv_write_dble(111,tp,.true.,'spc')
            END DO
            icmpnt = ndstaID(ithprtl) - 1
            IF (MeshType == "L") THEN
                DO k = elstaID(ithprtl), elendID(ithprtl)
                    CALL csv_write_integer(111,elmntlnknd(k,1)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,2)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,3)-icmpnt,.true.,'spc')
                    END DO
            END IF
            IF (MeshType == "Q") THEN
                DO k = elstaID(ithprtl), elendID(ithprtl)
                    CALL csv_write_integer(111,elmntlnknd(k,1)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,2)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,3)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.true.,'spc')
                END DO
            END IF
        END DO

        CLOSE (111)

        INQUIRE (FILE="Rslt_Dmn2_plot.dat", exist=filexistOutput)
        IF (filexistOutput) THEN
            OPEN (111, FILE="Rslt_Dmn2_plot.dat", STATUS="OLD", &
                & POSITION="APPEND", ACTION="WRITE")
        ELSE
            OPEN (111, FILE="Rslt_Dmn2_plot.dat", STATUS="NEW", &
                & ACTION="WRITE")
            CALL csv_write_char(111,'Variables="x", "y", "z",',.false.,'spc')
            CALL csv_write_char(111,'"E2x","E2y","E2z","|E2|",',.false.,'spc')
            CALL csv_write_char(111,'"H2x","H2y","H2z","|H2|"',.false.,'spc')
            CALL csv_write_char(111,'"S2x","S2y","S2z","|S2|"',.true.,'spc')
        END IF

        CALL csv_write_char(111,'Zone T ="time step',.false.,'spc')
        CALL csv_write_integer(111,dmksmcs,.false.,'spc')
        CALL csv_write_char(111,'plane - ',.false.,'spc')
        CALL csv_write_char(111,dmxyz_RMpln,.false.,'spc')
        CALL csv_write_char(111,'", n=',.false.,'spc')
        CALL csv_write_integer(111,nmbrnd_RMpln,.false.,'spc')
        CALL csv_write_char(111,', e=',.false.,'spc')
        CALL csv_write_integer(111,nmbrelmnt_RMpln,.false.,'spc')
        CALL csv_write_char(111,', f=fepoint, et=triangle',.true.,'spc')

        DO i = 1, nmbrnd_RMpln
            Dmnx = x_RMpln(i)
            Dmny = y_RMpln(i)
            Dmnz = z_RMpln(i)
            CALL csv_write_dble(111,Dmnx,.false.,'cmr')
            CALL csv_write_dble(111,Dmny,.false.,'cmr')
            CALL csv_write_dble(111,Dmnz,.false.,'cmr')
            CALL csv_write_dble(111,REAL(E2xdmn(i)),.false.,'cmr')
            CALL csv_write_dble(111,REAL(E2ydmn(i)),.false.,'cmr')
            CALL csv_write_dble(111,REAL(E2zdmn(i)),.false.,'cmr')
            tp = DSQRT( CDABS(E2xdmn(i))**2 &
            &          +CDABS(E2ydmn(i))**2 &
            &          +CDABS(E2zdmn(i))**2 )
            CALL csv_write_dble(111,tp,.false.,'cmr')
            CALL csv_write_dble(111,REAL(H2xdmn(i)),.false.,'cmr')
            CALL csv_write_dble(111,REAL(H2ydmn(i)),.false.,'cmr')
            CALL csv_write_dble(111,REAL(H2zdmn(i)),.false.,'cmr')
            tp = DSQRT( CDABS(H2xdmn(i))**2 &
            &          +CDABS(H2ydmn(i))**2 &
            &          +CDABS(H2zdmn(i))**2 )
            CALL csv_write_dble(111,tp,.false.,'cmr')
            ztp7 = 0.5d0*(E2ydmn(i)*DCONJG(H2zdmn(i)) - E2zdmn(i)*DCONJG(H2ydmn(i)))
            ztp8 = 0.5d0*(E2zdmn(i)*DCONJG(H2xdmn(i)) - E2xdmn(i)*DCONJG(H2zdmn(i)))
            ztp9 = 0.5d0*(E2xdmn(i)*DCONJG(H2ydmn(i)) - E2ydmn(i)*DCONJG(H2xdmn(i)))
            tp = DSQRT( CDABS(ztp7)**2 &
            &          +CDABS(ztp8)**2 &
            &          +CDABS(ztp9)**2 )
            CALL csv_write_dble(111,REAL(ztp7),.false.,'cmr')
            CALL csv_write_dble(111,REAL(ztp8),.false.,'cmr')
            CALL csv_write_dble(111,REAL(ztp9),.false.,'cmr')
            CALL csv_write_dble(111,tp,.true.,'spc')
        END DO

        DO i = 1, PostEgnd_RMpln - 1
            DO j = 1, dmEgnd_RMpln - 1
                CALL csv_write_integer(111,ndid_RMpln(i  ,j  ),.false.,'cmr')
                CALL csv_write_integer(111,ndid_RMpln(i+1,j  ),.false.,'cmr')
                CALL csv_write_integer(111,ndid_RMpln(i+1,j+1),.true.,'spc')

                CALL csv_write_integer(111,ndid_RMpln(i+1,j+1),.false.,'cmr')
                CALL csv_write_integer(111,ndid_RMpln(i  ,j+1),.false.,'cmr')
                CALL csv_write_integer(111,ndid_RMpln(i  ,j  ),.true.,'spc')
            END DO
        END DO

        CLOSE (111)

! Total field (E3,H3): ex and in surface zones followed by the plane zone.

        INQUIRE (FILE="Rslt_Dmn3_plot.dat", exist=filexistOutput)
        IF (filexistOutput) THEN
            OPEN (111, FILE="Rslt_Dmn3_plot.dat", STATUS="OLD", &
                & POSITION="APPEND", ACTION="WRITE")
        ELSE
            OPEN (111, FILE="Rslt_Dmn3_plot.dat", STATUS="NEW", &
                & ACTION="WRITE")
            CALL csv_write_char(111,'Variables="x", "y", "z",',.false.,'spc')
            CALL csv_write_char(111,'"E3x","E3y","E3z","|E3|",',.false.,'spc')
            CALL csv_write_char(111,'"H3x","H3y","H3z","|H3|"',.false.,'spc')
            CALL csv_write_char(111,'"S3x","S3y","S3z","|S3|"',.true.,'spc')
        END IF

        DO ithprtl = 1, nmbrprtl
            IF (MeshType == "L") THEN
                CALL csv_write_char(111,'Zone T ="time step',.false.,'spc')
                CALL csv_write_integer(111,dmksmcs,.false.,'spc')
                CALL csv_write_char(111,'Particle ',.false.,'spc')
                CALL csv_write_integer(111,ithprtl,.false.,'spc')
                CALL csv_write_char(111,' Out - ',.false.,'spc')
                CALL csv_write_char(111,dmxyz_RMpln,.false.,'spc')
                CALL csv_write_char(111,'", n=',.false.,'spc')
                CALL csv_write_integer(111,nmbrnd(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', e=',.false.,'spc')
                CALL csv_write_integer(111,nmbrelmnt(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', f=fepoint, et=triangle',.true.,'spc')
            END IF
            IF (MeshType == "Q") THEN
                CALL csv_write_char(111,'Zone T ="time step',.false.,'spc')
                CALL csv_write_integer(111,dmksmcs,.false.,'spc')
                CALL csv_write_char(111,'Particle ',.false.,'spc')
                CALL csv_write_integer(111,ithprtl,.false.,'spc')
                CALL csv_write_char(111,' Out - ',.false.,'spc')
                CALL csv_write_char(111,dmxyz_RMpln,.false.,'spc')
                CALL csv_write_char(111,'", n=',.false.,'spc')
                CALL csv_write_integer(111,nmbrnd(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', e=',.false.,'spc')
                CALL csv_write_integer(111,4*nmbrelmnt(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', f=fepoint, et=triangle',.true.,'spc')
            END IF

            DO i = ndstaID(ithprtl), ndendID(ithprtl)
                CALL csv_write_dble(111,xnd(i),.false.,'cmr')
                CALL csv_write_dble(111,ynd(i),.false.,'cmr')
                CALL csv_write_dble(111,znd(i),.false.,'cmr')
                CALL csv_write_dble(111,REAL(exE3x_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(exE3y_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(exE3z_EM(i)),.false.,'cmr')
                tp = DSQRT( CDABS(exE3x_EM(i))**2 &
                &          +CDABS(exE3y_EM(i))**2 &
                &          +CDABS(exE3z_EM(i))**2 )
                CALL csv_write_dble(111,tp,.false.,'cmr')
                CALL csv_write_dble(111,REAL(exH3x_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(exH3y_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(exH3z_EM(i)),.false.,'cmr')
                tp = DSQRT( CDABS(exH3x_EM(i))**2 &
                &          +CDABS(exH3y_EM(i))**2 &
                &          +CDABS(exH3z_EM(i))**2 )
                CALL csv_write_dble(111,tp,.false.,'cmr')
                ztp7 = 0.5d0*(exE3y_EM(i)*DCONJG(exH3z_EM(i)) - exE3z_EM(i)*DCONJG(exH3y_EM(i)))
                ztp8 = 0.5d0*(exE3z_EM(i)*DCONJG(exH3x_EM(i)) - exE3x_EM(i)*DCONJG(exH3z_EM(i)))
                ztp9 = 0.5d0*(exE3x_EM(i)*DCONJG(exH3y_EM(i)) - exE3y_EM(i)*DCONJG(exH3x_EM(i)))
                tp = DSQRT( CDABS(ztp7)**2 &
                &          +CDABS(ztp8)**2 &
                &          +CDABS(ztp9)**2 )
                CALL csv_write_dble(111,REAL(ztp7),.false.,'cmr')
                CALL csv_write_dble(111,REAL(ztp8),.false.,'cmr')
                CALL csv_write_dble(111,REAL(ztp9),.false.,'cmr')
                CALL csv_write_dble(111,tp,.true.,'spc')
            END DO
            icmpnt = ndstaID(ithprtl) - 1
            IF (MeshType == "L") THEN
                DO k = elstaID(ithprtl), elendID(ithprtl)
                    CALL csv_write_integer(111,elmntlnknd(k,1)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,2)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,3)-icmpnt,.true.,'spc')
                    END DO
            END IF
            IF (MeshType == "Q") THEN
                DO k = elstaID(ithprtl), elendID(ithprtl)
                    CALL csv_write_integer(111,elmntlnknd(k,1)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,2)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,3)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.true.,'spc')
                END DO
            END IF
        END DO

        CLOSE (111)

        INQUIRE (FILE="Rslt_Dmn3_plot.dat", exist=filexistOutput)
        IF (filexistOutput) THEN
            OPEN (111, FILE="Rslt_Dmn3_plot.dat", STATUS="OLD", &
                & POSITION="APPEND", ACTION="WRITE")
        ELSE
            OPEN (111, FILE="Rslt_Dmn3_plot.dat", STATUS="NEW", &
                & ACTION="WRITE")
            CALL csv_write_char(111,'Variables="x", "y", "z",',.false.,'spc')
            CALL csv_write_char(111,'"E3x","E3y","E3z","|E3|",',.false.,'spc')
            CALL csv_write_char(111,'"H3x","H3y","H3z","|H3|"',.false.,'spc')
            CALL csv_write_char(111,'"S3x","S3y","S3z","|S3|"',.true.,'spc')
        END IF

        DO ithprtl = 1, nmbrprtl
            IF (MeshType == "L") THEN
                CALL csv_write_char(111,'Zone T ="time step',.false.,'spc')
                CALL csv_write_integer(111,dmksmcs,.false.,'spc')
                CALL csv_write_char(111,'Particle ',.false.,'spc')
                CALL csv_write_integer(111,ithprtl,.false.,'spc')
                CALL csv_write_char(111,' In - ',.false.,'spc')
                CALL csv_write_char(111,dmxyz_RMpln,.false.,'spc')
                CALL csv_write_char(111,'", n=',.false.,'spc')
                CALL csv_write_integer(111,nmbrnd(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', e=',.false.,'spc')
                CALL csv_write_integer(111,nmbrelmnt(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', f=fepoint, et=triangle',.true.,'spc')
            END IF
            IF (MeshType == "Q") THEN
                CALL csv_write_char(111,'Zone T ="time step',.false.,'spc')
                CALL csv_write_integer(111,dmksmcs,.false.,'spc')
                CALL csv_write_char(111,'Particle ',.false.,'spc')
                CALL csv_write_integer(111,ithprtl,.false.,'spc')
                CALL csv_write_char(111,' In - ',.false.,'spc')
                CALL csv_write_char(111,dmxyz_RMpln,.false.,'spc')
                CALL csv_write_char(111,'", n=',.false.,'spc')
                CALL csv_write_integer(111,nmbrnd(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', e=',.false.,'spc')
                CALL csv_write_integer(111,4*nmbrelmnt(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', f=fepoint, et=triangle',.true.,'spc')
            END IF
            DO i = ndstaID(ithprtl), ndendID(ithprtl)
                CALL csv_write_dble(111,xnd(i),.false.,'cmr')
                CALL csv_write_dble(111,ynd(i),.false.,'cmr')
                CALL csv_write_dble(111,znd(i),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inE3x_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inE3y_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inE3z_EM(i)),.false.,'cmr')
                tp = DSQRT( CDABS(inE3x_EM(i))**2 &
                &          +CDABS(inE3y_EM(i))**2 &
                &          +CDABS(inE3z_EM(i))**2 )
                CALL csv_write_dble(111,tp,.false.,'cmr')
                CALL csv_write_dble(111,REAL(inH3x_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inH3y_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inH3z_EM(i)),.false.,'cmr')
                tp = DSQRT( CDABS(inH3x_EM(i))**2 &
                &          +CDABS(inH3y_EM(i))**2 &
                &          +CDABS(inH3z_EM(i))**2 )
                CALL csv_write_dble(111,tp,.false.,'cmr')
                ztp7 = 0.5d0*(inE3y_EM(i)*DCONJG(inH3z_EM(i)) - inE3z_EM(i)*DCONJG(inH3y_EM(i)))
                ztp8 = 0.5d0*(inE3z_EM(i)*DCONJG(inH3x_EM(i)) - inE3x_EM(i)*DCONJG(inH3z_EM(i)))
                ztp9 = 0.5d0*(inE3x_EM(i)*DCONJG(inH3y_EM(i)) - inE3y_EM(i)*DCONJG(inH3x_EM(i)))
                tp = DSQRT( CDABS(ztp7)**2 &
                &          +CDABS(ztp8)**2 &
                &          +CDABS(ztp9)**2 )
                CALL csv_write_dble(111,REAL(ztp7),.false.,'cmr')
                CALL csv_write_dble(111,REAL(ztp8),.false.,'cmr')
                CALL csv_write_dble(111,REAL(ztp9),.false.,'cmr')
                CALL csv_write_dble(111,tp,.true.,'spc')
            END DO
            icmpnt = ndstaID(ithprtl) - 1
            IF (MeshType == "L") THEN
                DO k = elstaID(ithprtl), elendID(ithprtl)
                    CALL csv_write_integer(111,elmntlnknd(k,1)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,2)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,3)-icmpnt,.true.,'spc')
                    END DO
            END IF
            IF (MeshType == "Q") THEN
                DO k = elstaID(ithprtl), elendID(ithprtl)
                    CALL csv_write_integer(111,elmntlnknd(k,1)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,2)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,3)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.false.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.true.,'spc')
                END DO
            END IF
        END DO

        CLOSE (111)

        INQUIRE (FILE="Rslt_Dmn3_plot.dat", exist=filexistOutput)
        IF (filexistOutput) THEN
            OPEN (111, FILE="Rslt_Dmn3_plot.dat", STATUS="OLD", &
                & POSITION="APPEND", ACTION="WRITE")
        ELSE
            OPEN (111, FILE="Rslt_Dmn3_plot.dat", STATUS="NEW", &
                & ACTION="WRITE")
            CALL csv_write_char(111,'Variables="x", "y", "z",',.false.,'spc')
            CALL csv_write_char(111,'"E3x","E3y","E3z","|E3|",',.false.,'spc')
            CALL csv_write_char(111,'"H3x","H3y","H3z","|H3|"',.false.,'spc')
            CALL csv_write_char(111,'"S3x","S3y","S3z","|S3|"',.true.,'spc')
        END IF

        CALL csv_write_char(111,'Zone T ="time step',.false.,'spc')
        CALL csv_write_integer(111,dmksmcs,.false.,'spc')
        CALL csv_write_char(111,'plane - ',.false.,'spc')
        CALL csv_write_char(111,dmxyz_RMpln,.false.,'spc')
        CALL csv_write_char(111,'", n=',.false.,'spc')
        CALL csv_write_integer(111,nmbrnd_RMpln,.false.,'spc')
        CALL csv_write_char(111,', e=',.false.,'spc')
        CALL csv_write_integer(111,nmbrelmnt_RMpln,.false.,'spc')
        CALL csv_write_char(111,', f=fepoint, et=triangle',.true.,'spc')

        DO i = 1, nmbrnd_RMpln
            Dmnx = x_RMpln(i)
            Dmny = y_RMpln(i)
            Dmnz = z_RMpln(i)
            CALL csv_write_dble(111,Dmnx,.false.,'cmr')
            CALL csv_write_dble(111,Dmny,.false.,'cmr')
            CALL csv_write_dble(111,Dmnz,.false.,'cmr')
            ztp1 = E1xdmn(i)+E2xdmn(i)
            ztp2 = E1ydmn(i)+E2ydmn(i)
            ztp3 = E1zdmn(i)+E2zdmn(i)
            CALL csv_write_dble(111,REAL(ztp1),.false.,'cmr')
            CALL csv_write_dble(111,REAL(ztp2),.false.,'cmr')
            CALL csv_write_dble(111,REAL(ztp3),.false.,'cmr')
            tp = DSQRT( CDABS(ztp1)**2 &
            &          +CDABS(ztp2)**2 &
            &          +CDABS(ztp3)**2 )
            CALL csv_write_dble(111,tp,.false.,'cmr')
            ztp4 = H1xdmn(i)+H2xdmn(i)
            ztp5 = H1ydmn(i)+H2ydmn(i)
            ztp6 = H1zdmn(i)+H2zdmn(i)
            CALL csv_write_dble(111,REAL(ztp4),.false.,'cmr')
            CALL csv_write_dble(111,REAL(ztp5),.false.,'cmr')
            CALL csv_write_dble(111,REAL(ztp6),.false.,'cmr')
            tp = DSQRT( CDABS(ztp4)**2 &
            &          +CDABS(ztp5)**2 &
            &          +CDABS(ztp6)**2 )
            CALL csv_write_dble(111,tp,.false.,'cmr')
            ztp7 = 0.5d0*(ztp2*DCONJG(ztp6) - ztp3*DCONJG(ztp5))
            ztp8 = 0.5d0*(ztp3*DCONJG(ztp4) - ztp1*DCONJG(ztp6))
            ztp9 = 0.5d0*(ztp1*DCONJG(ztp5) - ztp2*DCONJG(ztp4))
            tp = DSQRT( CDABS(ztp7)**2 &
            &          +CDABS(ztp8)**2 &
            &          +CDABS(ztp9)**2 )
            CALL csv_write_dble(111,REAL(ztp7),.false.,'cmr')
            CALL csv_write_dble(111,REAL(ztp8),.false.,'cmr')
            CALL csv_write_dble(111,REAL(ztp9),.false.,'cmr')
            CALL csv_write_dble(111,tp,.true.,'spc')
        END DO

        DO i = 1, PostEgnd_RMpln - 1
            DO j = 1, dmEgnd_RMpln - 1
                CALL csv_write_integer(111,ndid_RMpln(i  ,j  ),.false.,'cmr')
                CALL csv_write_integer(111,ndid_RMpln(i+1,j  ),.false.,'cmr')
                CALL csv_write_integer(111,ndid_RMpln(i+1,j+1),.true.,'spc')

                CALL csv_write_integer(111,ndid_RMpln(i+1,j+1),.false.,'cmr')
                CALL csv_write_integer(111,ndid_RMpln(i  ,j+1),.false.,'cmr')
                CALL csv_write_integer(111,ndid_RMpln(i  ,j  ),.true.,'spc')
            END DO
        END DO

        CLOSE (111)

        DEALLOCATE (x_RMpln, y_RMpln, z_RMpln, pos_RMpln, ndid_RMpln)
        DEALLOCATE (E1xdmn,E1ydmn,E1zdmn,H1xdmn,H1ydmn,H1zdmn, &
        &           E2xdmn,E2ydmn,E2zdmn,H2xdmn,H2ydmn,H2zdmn)

    END SUBROUTINE

END MODULE
