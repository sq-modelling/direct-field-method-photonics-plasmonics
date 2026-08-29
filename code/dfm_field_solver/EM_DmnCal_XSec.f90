
! SPDX-FileCopyrightText: 2026 Qiang Sun
! SPDX-License-Identifier: BSD-3-Clause

MODULE EM_DmnCal_XSec

    USE omp_lib

    USE Pre_Constants
    USE Pre_csvformat

    USE Geom_GlobalData
    USE Geom_MeshSphereCircle

    USE EM_SurfCal_GlobalData
    USE EM_SurfCal_PhysBC
    USE EM_DmnCal


    IMPLICIT NONE
    CONTAINS


    SUBROUTINE GetRCS_EM(dmkcase, dmnumpoint, dmsize, dmxoff, dmyoff, dmzoff)

        INTEGER, INTENT (IN) :: dmkcase, dmnumpoint

        DOUBLE PRECISION, INTENT (IN) :: dmsize, dmxoff, dmyoff, dmzoff

        LOGICAL :: filexist

        INTEGER :: icmpnt, i
        DOUBLE PRECISION :: tp,tp1,tp2,Dmnx,Dmny,Dmnz
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: tht_RMpln,x_RMpln,y_RMpln,z_RMpln,rcs_RMpln
        COMPLEX(KIND=KIND(1.0D0)) :: ztpEx,ztpEy,ztpEz,ztpHx,ztpHy,ztpHz,&
        &                            Exdmn,Eydmn,Ezdmn,Hxdmn,Hydmn,Hzdmn

        DO icmpnt = 1, 3

            ALLOCATE (tht_RMpln(dmnumpoint))
            ALLOCATE (x_RMpln(dmnumpoint))
            ALLOCATE (y_RMpln(dmnumpoint))
            ALLOCATE (z_RMpln(dmnumpoint))
            ALLOCATE (rcs_RMpln(dmnumpoint))

            tp = 2.0d0*pai/dble(dmnumpoint-1)

            DO i = 1, dmnumpoint

                tht_RMpln(i) = dble(i-1)*tp

                IF (icmpnt == 1) THEN
                    x_RMpln(i) = dmsize*dcos(dble(i-1)*tp) + dmxoff
                    y_RMpln(i) = dmsize*dsin(dble(i-1)*tp) + dmyoff
                    z_RMpln(i) = 0.0d0                     + dmzoff
                END IF

                IF (icmpnt == 2) THEN
                    x_RMpln(i) = 0.0d0                     + dmxoff
                    y_RMpln(i) = dmsize*dcos(dble(i-1)*tp) + dmyoff
                    z_RMpln(i) = dmsize*dsin(dble(i-1)*tp) + dmzoff
                END IF

                IF (icmpnt == 3) THEN
                    x_RMpln(i) = dmsize*dcos(dble(i-1)*tp) + dmxoff
                    y_RMpln(i) = 0.0d0                     + dmyoff
                    z_RMpln(i) = dmsize*dsin(dble(i-1)*tp) + dmzoff
                END IF

            END DO

            DO i = 1, dmnumpoint

                Dmnx = x_RMpln(i)
                Dmny = y_RMpln(i)
                Dmnz = z_RMpln(i)

                CALL GetDmnE1H1_EM (Dmnx,Dmny,Dmnz,0,ztpEx,ztpEy,ztpEz,ztpHx,ztpHy,ztpHz)
                CALL GetDmnE2H2_EM ('E',Dmnx,Dmny,Dmnz,0,Exdmn,Eydmn,Ezdmn)

                tp1 = REAL( Exdmn*DCONJG(Exdmn) + Eydmn*DCONJG(Eydmn) + Ezdmn*DCONJG(Ezdmn) )
                tp2 = REAL( ztpEx*DCONJG(ztpEx) + ztpEy*DCONJG(ztpEy) + ztpEz*DCONJG(ztpEz) )

                rcs_RMpln(i) = 4.0d0*pai * (dmsize**2) * tp1/tp2

            END DO

            IF (icmpnt == 1) THEN
                INQUIRE (FILE="Rslt_RCS_xy.dat", exist=filexist)
                IF (filexist) THEN
                    OPEN (111, FILE="Rslt_RCS_xy.dat", STATUS="OLD", &
                        & POSITION="APPEND", ACTION="WRITE")
                ELSE
                    OPEN (111, &
                    &FILE='Rslt_RCS_xy.dat', &
                    &STATUS="NEW", ACTION="WRITE")
                    CALL csv_write_char(111,&
                    &'Variables="theta", "x", "y", "z", "RCS","RCS(db_sm)"',&
                    &.true.,'spc')
                END IF
            END IF
            IF (icmpnt == 2) THEN
                INQUIRE (FILE="Rslt_RCS_yz.dat", exist=filexist)
                IF (filexist) THEN
                    OPEN (111, FILE="Rslt_RCS_yz.dat", STATUS="OLD", &
                        & POSITION="APPEND", ACTION="WRITE")
                ELSE
                    OPEN (111, &
                    &FILE='Rslt_RCS_yz.dat', &
                    &STATUS="NEW", ACTION="WRITE")
                    CALL csv_write_char(111,&
                    &'Variables="theta", "x", "y", "z", "RCS","RCS(db_sm)"',&
                    &.true.,'spc')
                END IF
            END IF
            IF (icmpnt == 3) THEN
                INQUIRE (FILE="Rslt_RCS_zx.dat", exist=filexist)
                IF (filexist) THEN
                    OPEN (111, FILE="Rslt_RCS_zx.dat", STATUS="OLD", &
                        & POSITION="APPEND", ACTION="WRITE")
                ELSE
                    OPEN (111, &
                    &FILE='Rslt_RCS_zx.dat', &
                    &STATUS="NEW", ACTION="WRITE")
                    CALL csv_write_char(111,&
                    &'Variables="theta", "x", "y", "z", "RCS","RCS(db_sm)"',&
                    &.true.,'spc')
                END IF
            END IF

            CALL csv_write_char(111,'Zone T ="Case ',.false.,'spc')
            CALL csv_write_integer(111,dmkcase,.false.,'spc')
            CALL csv_write_char(111,'"',.true.,'spc')

            DO i = 1, dmnumpoint
                Dmnx = x_RMpln(i)
                Dmny = y_RMpln(i)
                Dmnz = z_RMpln(i)
                tp = tht_RMpln(i) / pai * 180.0d0
                CALL csv_write_dble(111,tp,.false.,'cmr')
                CALL csv_write_dble(111,Dmnx,.false.,'cmr')
                CALL csv_write_dble(111,Dmny,.false.,'cmr')
                CALL csv_write_dble(111,Dmnz,.false.,'cmr')
                CALL csv_write_dble(111,rcs_RMpln(i),.false.,'cmr')

                tp = 10.0d0*dLog10(rcs_RMpln(i))

                CALL csv_write_dble(111,tp,.true.,'spc')

            END DO

            CLOSE (111)

            DEALLOCATE (tht_RMpln, x_RMpln, y_RMpln, z_RMpln,rcs_RMpln)

        END DO

    END SUBROUTINE

END MODULE
