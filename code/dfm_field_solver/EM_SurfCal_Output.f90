
! SPDX-FileCopyrightText: 2026 Qiang Sun
! SPDX-License-Identifier: BSD-3-Clause

MODULE EM_SurfCal_Output

! Write solved surface traces as Tecplot ASCII FEPOINT zones.  Field names and the
! ex/in, E1/E2/E3, phasor, and normal conventions are defined in
! EM_SurfCal_GlobalData.

    USE omp_lib

    USE Pre_Constants
    USE Pre_csvformat

    USE Geom_GlobalData

    USE EM_SurfCal_GlobalData


    IMPLICIT NONE

    CONTAINS

    SUBROUTINE RsltSurf_fieldCase(dmkcase)

! Purpose: append one zone per surface to each of the incident/source (1), solved
!          scattered/correction (2), and total (3) surface-field files.
! Input  : dmkcase labels the appended zones; the solved mesh and surface traces must
!          already be available in Geom_GlobalData and EM_SurfCal_GlobalData.
! Output : Rslt_SurfCal1_EM.dat, Rslt_SurfCal2_EM.dat, and Rslt_SurfCal3_EM.dat.
!          Existing files are appended; otherwise a Tecplot variable header is written.
! Columns: x,y,z use the mesh length unit.  For each ex/in E or H trace, x/y/z are
!          REAL parts of the complex phasor, |F| is sqrt(sum_j |F_j|**2), and the
!          header labels Se/Sh denote REAL(F.n), not a Poynting vector.  Field units
!          are the solver's native units; no output conversion is applied.
! Mesh   : linear triangles are written directly; each six-node quadratic triangle is
!          represented by four three-node plotting triangles without changing its data.

        INTEGER, INTENT (IN) :: dmkcase

        LOGICAL :: filexist
        INTEGER :: i, j, k, ithprtl, jthprtl, kthprtl, icmpnt, id_tp
        DOUBLE PRECISION :: tp, tp0, tp1, tp2, tp3, tp4, tp5, tp6, tp7, tp8, tp9
        COMPLEX(KIND=KIND(1.0D0)) :: ztp,ztp1,ztp2,ztp3,ztp4,ztp5,ztp6,ztp7,ztp8,ztp9

! Incident/source-field surface traces (E1 and H1).

        INQUIRE (FILE="Rslt_SurfCal1_EM.dat", exist=filexist)
        IF (filexist) THEN
            OPEN (111, FILE="Rslt_SurfCal1_EM.dat", STATUS="OLD", &
                & POSITION="APPEND", ACTION="WRITE")
        ELSE
            OPEN (111, FILE="Rslt_SurfCal1_EM.dat", STATUS="NEW", &
                & ACTION="WRITE")
            CALL csv_write_char(111,'Variables="x", "y", "z",',.false.,'spc')
            CALL csv_write_char(111,'"exE1x", "exE1y", "exE1z", "|exE1|", "exSe1", ',.false.,'spc')
            CALL csv_write_char(111,'"inE1x", "inE1y", "inE1z", "|inE1|", "inSe1", ',.false.,'spc')
            CALL csv_write_char(111,'"exH1x", "exH1y", "exH1z", "|exH1|", "exSh1", ',.false.,'spc')
            CALL csv_write_char(111,'"inH1x", "inH1y", "inH1z", "|inH1|", "inSh1"  ',.true., 'spc')
        END IF

        DO ithprtl = 1, nmbrprtl
            IF (MeshType == "L") THEN
                CALL csv_write_char(111,'Zone T ="Input - case ',.false.,'spc')
                CALL csv_write_integer(111,dmkcase,.false.,'spc')
                CALL csv_write_char(111,' Particle ',.false.,'spc')
                CALL csv_write_integer(111,ithprtl,.false.,'spc')
                CALL csv_write_char(111,'", n=',.false.,'spc')
                CALL csv_write_integer(111,nmbrnd(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', e=',.false.,'spc')
                CALL csv_write_integer(111,nmbrelmnt(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', f=fepoint, et=triangle',.true.,'spc')
            END IF
            IF (MeshType == "Q") THEN
                CALL csv_write_char(111,'Zone T ="Input - case ',.false.,'spc')
                CALL csv_write_integer(111,dmkcase,.false.,'spc')
                CALL csv_write_char(111,' Particle ',.false.,'spc')
                CALL csv_write_integer(111,ithprtl,.false.,'spc')
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
                tp = REAL(  exE1x_EM(i)*nnx(i) &
                &          +exE1y_EM(i)*nny(i) &
                &          +exE1z_EM(i)*nnz(i) )
                CALL csv_write_dble(111,tp,.false.,'cmr')

                CALL csv_write_dble(111,REAL(inE1x_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inE1y_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inE1z_EM(i)),.false.,'cmr')
                tp = DSQRT( CDABS(inE1x_EM(i))**2 &
                &          +CDABS(inE1y_EM(i))**2 &
                &          +CDABS(inE1z_EM(i))**2 )
                CALL csv_write_dble(111,tp,.false.,'cmr')
                tp = REAL(  inE1x_EM(i)*nnx(i) &
                &          +inE1y_EM(i)*nny(i) &
                &          +inE1z_EM(i)*nnz(i) )
                CALL csv_write_dble(111,tp,.false.,'cmr')

                CALL csv_write_dble(111,REAL(exH1x_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(exH1y_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(exH1z_EM(i)),.false.,'cmr')
                tp = DSQRT( CDABS(exH1x_EM(i))**2 &
                &          +CDABS(exH1y_EM(i))**2 &
                &          +CDABS(exH1z_EM(i))**2 )
                CALL csv_write_dble(111,tp,.false.,'cmr')
                tp = REAL(  exH1x_EM(i)*nnx(i) &
                &          +exH1y_EM(i)*nny(i) &
                &          +exH1z_EM(i)*nnz(i) )
                CALL csv_write_dble(111,tp,.false.,'cmr')

                CALL csv_write_dble(111,REAL(inH1x_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inH1y_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inH1z_EM(i)),.false.,'cmr')
                tp = DSQRT( CDABS(inH1x_EM(i))**2 &
                &          +CDABS(inH1y_EM(i))**2 &
                &          +CDABS(inH1z_EM(i))**2 )
                CALL csv_write_dble(111,tp,.false.,'cmr')
                tp = REAL(  inH1x_EM(i)*nnx(i) &
                &          +inH1y_EM(i)*nny(i) &
                &          +inH1z_EM(i)*nnz(i) )
                CALL csv_write_dble(111,tp,.true., 'spc')

            END DO
            icmpnt = ndstaID(ithprtl) - 1
            IF (MeshType == "L") THEN
                DO k = elstaID(ithprtl), elendID(ithprtl)
                    CALL csv_write_integer(111,elmntlnknd(k,1)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,2)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,3)-icmpnt,.true.,'spc')
                END DO
            END IF
            IF (MeshType == "Q") THEN
                DO k = elstaID(ithprtl), elendID(ithprtl)
                    CALL csv_write_integer(111,elmntlnknd(k,1)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,2)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,3)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.true.,'spc')
                END DO
            END IF
        END DO

        CLOSE (111)

! Solved scattered/correction-field surface traces (E2 and H2).

        INQUIRE (FILE="Rslt_SurfCal2_EM.dat", exist=filexist)
        IF (filexist) THEN
            OPEN (111, FILE="Rslt_SurfCal2_EM.dat", STATUS="OLD", &
                & POSITION="APPEND", ACTION="WRITE")
        ELSE
            OPEN (111, FILE="Rslt_SurfCal2_EM.dat", STATUS="NEW", &
                & ACTION="WRITE")
            CALL csv_write_char(111,'Variables="x", "y", "z",',.false.,'spc')
            CALL csv_write_char(111,'"exE2x", "exE2y", "exE2z", "|exE2|", "exSe2", ',.false.,'spc')
            CALL csv_write_char(111,'"inE2x", "inE2y", "inE2z", "|inE2|", "inSe2", ',.false.,'spc')
            CALL csv_write_char(111,'"exH2x", "exH2y", "exH2z", "|exH2|", "exSh2", ',.false.,'spc')
            CALL csv_write_char(111,'"inH2x", "inH2y", "inH2z", "|inH2|", "inSh2"  ',.true., 'spc')
        END IF

        DO ithprtl = 1, nmbrprtl
            IF (MeshType == "L") THEN
                CALL csv_write_char(111,'Zone T ="Input - case ',.false.,'spc')
                CALL csv_write_integer(111,dmkcase,.false.,'spc')
                CALL csv_write_char(111,' Particle ',.false.,'spc')
                CALL csv_write_integer(111,ithprtl,.false.,'spc')
                CALL csv_write_char(111,'", n=',.false.,'spc')
                CALL csv_write_integer(111,nmbrnd(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', e=',.false.,'spc')
                CALL csv_write_integer(111,nmbrelmnt(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', f=fepoint, et=triangle',.true.,'spc')
            END IF
            IF (MeshType == "Q") THEN
                CALL csv_write_char(111,'Zone T ="Input - case ',.false.,'spc')
                CALL csv_write_integer(111,dmkcase,.false.,'spc')
                CALL csv_write_char(111,' Particle ',.false.,'spc')
                CALL csv_write_integer(111,ithprtl,.false.,'spc')
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
                tp = REAL(  exE2x_EM(i)*nnx(i) &
                &          +exE2y_EM(i)*nny(i) &
                &          +exE2z_EM(i)*nnz(i) )
                CALL csv_write_dble(111,tp,.false.,'cmr')

                CALL csv_write_dble(111,REAL(inE2x_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inE2y_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inE2z_EM(i)),.false.,'cmr')
                tp = DSQRT( CDABS(inE2x_EM(i))**2 &
                &          +CDABS(inE2y_EM(i))**2 &
                &          +CDABS(inE2z_EM(i))**2 )
                CALL csv_write_dble(111,tp,.false.,'cmr')
                tp = REAL(  inE2x_EM(i)*nnx(i) &
                &          +inE2y_EM(i)*nny(i) &
                &          +inE2z_EM(i)*nnz(i) )
                CALL csv_write_dble(111,tp,.false.,'cmr')

                CALL csv_write_dble(111,REAL(exH2x_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(exH2y_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(exH2z_EM(i)),.false.,'cmr')
                tp = DSQRT( CDABS(exH2x_EM(i))**2 &
                &          +CDABS(exH2y_EM(i))**2 &
                &          +CDABS(exH2z_EM(i))**2 )
                CALL csv_write_dble(111,tp,.false.,'cmr')
                tp = REAL(  exH2x_EM(i)*nnx(i) &
                &          +exH2y_EM(i)*nny(i) &
                &          +exH2z_EM(i)*nnz(i) )
                CALL csv_write_dble(111,tp,.false.,'cmr')

                CALL csv_write_dble(111,REAL(inH2x_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inH2y_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inH2z_EM(i)),.false.,'cmr')
                tp = DSQRT( CDABS(inH2x_EM(i))**2 &
                &          +CDABS(inH2y_EM(i))**2 &
                &          +CDABS(inH2z_EM(i))**2 )
                CALL csv_write_dble(111,tp,.false.,'cmr')
                tp = REAL(  inH2x_EM(i)*nnx(i) &
                &          +inH2y_EM(i)*nny(i) &
                &          +inH2z_EM(i)*nnz(i) )
                CALL csv_write_dble(111,tp,.true., 'spc')

            END DO
            icmpnt = ndstaID(ithprtl) - 1
            IF (MeshType == "L") THEN
                DO k = elstaID(ithprtl), elendID(ithprtl)
                    CALL csv_write_integer(111,elmntlnknd(k,1)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,2)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,3)-icmpnt,.true.,'spc')
                END DO
            END IF
            IF (MeshType == "Q") THEN
                DO k = elstaID(ithprtl), elendID(ithprtl)
                    CALL csv_write_integer(111,elmntlnknd(k,1)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,2)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,3)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.true.,'spc')
                END DO
            END IF
        END DO

        CLOSE (111)

! Total-field surface traces (E3 and H3).

        INQUIRE (FILE="Rslt_SurfCal3_EM.dat", exist=filexist)
        IF (filexist) THEN
            OPEN (111, FILE="Rslt_SurfCal3_EM.dat", STATUS="OLD", &
                & POSITION="APPEND", ACTION="WRITE")
        ELSE
            OPEN (111, FILE="Rslt_SurfCal3_EM.dat", STATUS="NEW", &
                & ACTION="WRITE")
            CALL csv_write_char(111,'Variables="x", "y", "z",',.false.,'spc')
            CALL csv_write_char(111,'"exE3x", "exE3y", "exE3z", "|exE3|", "exSe3", ',.false.,'spc')
            CALL csv_write_char(111,'"inE3x", "inE3y", "inE3z", "|inE3|", "inSe3", ',.false.,'spc')
            CALL csv_write_char(111,'"exH3x", "exH3y", "exH3z", "|exH3|", "exSh3", ',.false.,'spc')
            CALL csv_write_char(111,'"inH3x", "inH3y", "inH3z", "|inH3|", "inSh3"  ',.true., 'spc')
        END IF

        DO ithprtl = 1, nmbrprtl
            IF (MeshType == "L") THEN
                CALL csv_write_char(111,'Zone T ="Input - case ',.false.,'spc')
                CALL csv_write_integer(111,dmkcase,.false.,'spc')
                CALL csv_write_char(111,' Particle ',.false.,'spc')
                CALL csv_write_integer(111,ithprtl,.false.,'spc')
                CALL csv_write_char(111,'", n=',.false.,'spc')
                CALL csv_write_integer(111,nmbrnd(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', e=',.false.,'spc')
                CALL csv_write_integer(111,nmbrelmnt(ithprtl),.false.,'spc')
                CALL csv_write_char(111,', f=fepoint, et=triangle',.true.,'spc')
            END IF
            IF (MeshType == "Q") THEN
                CALL csv_write_char(111,'Zone T ="Input - case ',.false.,'spc')
                CALL csv_write_integer(111,dmkcase,.false.,'spc')
                CALL csv_write_char(111,' Particle ',.false.,'spc')
                CALL csv_write_integer(111,ithprtl,.false.,'spc')
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
                tp = REAL(  exE3x_EM(i)*nnx(i) &
                &          +exE3y_EM(i)*nny(i) &
                &          +exE3z_EM(i)*nnz(i) )
! For a PEM boundary, the existing Se3/Sh3 columns both receive the real normal
! component of eps*E3 - mu*BCValue*H3, as required by this output convention.
                IF (BCType_EM(ithprtl) == 'PEM') THEN
                    id_tp = corelnkshell(ithprtl)
                    IF (id_tp == 0) THEN
                        ztp1 = vcm_eps0*exeps_EM
                        ztp2 = vcm_mu0*exmiu_EM
                    ELSE
                        ztp1 = vcm_eps0*ineps_EM(id_tp)
                        ztp2 = vcm_mu0*inmiu_EM(id_tp)
                    END IF
                    tp = REAL(  (ztp1*exE3x_EM(i) - ztp2*exH3x_EM(i) * BCValue_EM(ithprtl))*nnx(i) &
                    &          +(ztp1*exE3y_EM(i) - ztp2*exH3y_EM(i) * BCValue_EM(ithprtl))*nny(i) &
                    &          +(ztp1*exE3z_EM(i) - ztp2*exH3z_EM(i) * BCValue_EM(ithprtl))*nnz(i) )
                END IF
                CALL csv_write_dble(111,tp,.false.,'cmr')

                CALL csv_write_dble(111,REAL(inE3x_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inE3y_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inE3z_EM(i)),.false.,'cmr')
                tp = DSQRT( CDABS(inE3x_EM(i))**2 &
                &          +CDABS(inE3y_EM(i))**2 &
                &          +CDABS(inE3z_EM(i))**2 )
                CALL csv_write_dble(111,tp,.false.,'cmr')
                tp = REAL(  inE3x_EM(i)*nnx(i) &
                &          +inE3y_EM(i)*nny(i) &
                &          +inE3z_EM(i)*nnz(i) )
                CALL csv_write_dble(111,tp,.false.,'cmr')

                CALL csv_write_dble(111,REAL(exH3x_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(exH3y_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(exH3z_EM(i)),.false.,'cmr')
                tp = DSQRT( CDABS(exH3x_EM(i))**2 &
                &          +CDABS(exH3y_EM(i))**2 &
                &          +CDABS(exH3z_EM(i))**2 )
                CALL csv_write_dble(111,tp,.false.,'cmr')
                tp = REAL(  exH3x_EM(i)*nnx(i) &
                &          +exH3y_EM(i)*nny(i) &
                &          +exH3z_EM(i)*nnz(i) )
                IF (BCType_EM(ithprtl) == 'PEM') THEN
                    id_tp = corelnkshell(ithprtl)
                    IF (id_tp == 0) THEN
                        ztp1 = vcm_eps0*exeps_EM
                        ztp2 = vcm_mu0*exmiu_EM
                    ELSE
                        ztp1 = vcm_eps0*ineps_EM(id_tp)
                        ztp2 = vcm_mu0*inmiu_EM(id_tp)
                    END IF
                    tp = REAL(  (ztp1*exE3x_EM(i) - ztp2*exH3x_EM(i) * BCValue_EM(ithprtl))*nnx(i) &
                    &          +(ztp1*exE3y_EM(i) - ztp2*exH3y_EM(i) * BCValue_EM(ithprtl))*nny(i) &
                    &          +(ztp1*exE3z_EM(i) - ztp2*exH3z_EM(i) * BCValue_EM(ithprtl))*nnz(i) )
                END IF
                CALL csv_write_dble(111,tp,.false.,'cmr')

                CALL csv_write_dble(111,REAL(inH3x_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inH3y_EM(i)),.false.,'cmr')
                CALL csv_write_dble(111,REAL(inH3z_EM(i)),.false.,'cmr')
                tp = DSQRT( CDABS(inH3x_EM(i))**2 &
                &          +CDABS(inH3y_EM(i))**2 &
                &          +CDABS(inH3z_EM(i))**2 )
                CALL csv_write_dble(111,tp,.false.,'cmr')
                tp = REAL(  inH3x_EM(i)*nnx(i) &
                &          +inH3y_EM(i)*nny(i) &
                &          +inH3z_EM(i)*nnz(i) )
                CALL csv_write_dble(111,tp,.true., 'spc')

            END DO
            icmpnt = ndstaID(ithprtl) - 1
            IF (MeshType == "L") THEN
                DO k = elstaID(ithprtl), elendID(ithprtl)
                    CALL csv_write_integer(111,elmntlnknd(k,1)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,2)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,3)-icmpnt,.true.,'spc')
                END DO
            END IF
            IF (MeshType == "Q") THEN
                DO k = elstaID(ithprtl), elendID(ithprtl)
                    CALL csv_write_integer(111,elmntlnknd(k,1)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,2)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,3)-icmpnt,.true.,'spc')
                    CALL csv_write_integer(111,elmntlnknd(k,4)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,5)-icmpnt,.false.,'cmr')
                    CALL csv_write_integer(111,elmntlnknd(k,6)-icmpnt,.true.,'spc')
                END DO
            END IF
        END DO

        CLOSE (111)

    END SUBROUTINE

END MODULE
