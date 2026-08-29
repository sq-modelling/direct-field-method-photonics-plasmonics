
! SPDX-FileCopyrightText: 2026 Qiang Sun
! SPDX-License-Identifier: BSD-3-Clause

MODULE EM_SurfCal_VariableInt

    USE omp_lib

    USE Pre_Constants
    USE Pre_csvformat

    USE Geom_GlobalData

    USE EM_SurfCal_GlobalData

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE GetVariableInt_EM

        IMPLICIT NONE

        INTEGER :: i, j, ithprtl, IOS
        DOUBLE PRECISION :: tp

        IF (ttlnmbrsrc_EM > 0) THEN
            ALLOCATE(srcType_EM(ttlnmbrsrc_EM))
            ALLOCATE(srcStrength_EM(ttlnmbrsrc_EM))
            ALLOCATE(xsrc_EM(ttlnmbrsrc_EM))
            ALLOCATE(ysrc_EM(ttlnmbrsrc_EM))
            ALLOCATE(zsrc_EM(ttlnmbrsrc_EM))
            ALLOCATE(polxsrc_EM(ttlnmbrsrc_EM))
            ALLOCATE(polysrc_EM(ttlnmbrsrc_EM))
            ALLOCATE(polzsrc_EM(ttlnmbrsrc_EM))
            DO i = 1, ttlnmbrsrc_EM
                srcStrength_EM(i) = 0.0d0
                xsrc_EM(i) = 0.0d0
                ysrc_EM(i) = 0.0d0
                zsrc_EM(i) = 0.0d0
            END DO
        ELSE
            ALLOCATE(srcType_EM(1))
            ALLOCATE(srcStrength_EM(1))
            ALLOCATE(xsrc_EM(1))
            ALLOCATE(ysrc_EM(1))
            ALLOCATE(zsrc_EM(1))
            ALLOCATE(polxsrc_EM(1))
            ALLOCATE(polysrc_EM(1))
            ALLOCATE(polzsrc_EM(1))
            srcStrength_EM(1) = 0.0d0
            xsrc_EM(1) = 1.0d18
            ysrc_EM(1) = 1.0d18
            zsrc_EM(1) = 1.0d18
            polxsrc_EM(1) = 0.0d0
            polysrc_EM(1) = 0.0d0
            polzsrc_EM(1) = 0.0d0
            srcType_EM(1) = 2
        END IF

        IF (ttlnmbrsrc_EM > 0) THEN
            OPEN (111, file = 'Input_Source_EM.dat', STATUS = 'OLD', IOSTAT = IOS)
            IF (IOS /= 0) THEN
                PRINT*, "'Input_Source_EM.dat' does not exist! Please check!"
                STOP
            END IF
            DO ithprtl = 0, nmbrprtl
                IF (srcendID_EM(ithprtl) > 0) THEN
                    READ (111, *)
                    DO i = srcstaID_EM(ithprtl), srcendID_EM(ithprtl)
                        READ(111, *) j, &
                        &            xsrc_EM(i), ysrc_EM(i), zsrc_EM(i), &
                        &            polxsrc_EM(i), polysrc_EM(i), polzsrc_EM(i), &
                        &            srcStrength_EM(i), srcType_EM(i)
                        IF (ithprtl > 0) THEN
                            xsrc_EM(i) = xsrc_EM(i) + xloctn(ithprtl)
                            ysrc_EM(i) = ysrc_EM(i) + yloctn(ithprtl)
                            zsrc_EM(i) = zsrc_EM(i) + zloctn(ithprtl)
                            tp = polxsrc_EM(i)**2+polysrc_EM(i)**2+polzsrc_EM(i)**2
                            tp = 1.0d0/DSQRT(tp)
                            polxsrc_EM(i) = polxsrc_EM(i) * tp
                            polysrc_EM(i) = polysrc_EM(i) * tp
                            polzsrc_EM(i) = polzsrc_EM(i) * tp
                            polxsrc_EM(i) = polxsrc_EM(i) * srcStrength_EM(i)
                            polysrc_EM(i) = polysrc_EM(i) * srcStrength_EM(i)
                            polzsrc_EM(i) = polzsrc_EM(i) * srcStrength_EM(i)
                        END IF
                    END DO
                    READ (111, *)
                END IF
            END DO
            CLOSE (111)
        END IF

        ALLOCATE(exE1x_EM(ttlnmbrnd))
        ALLOCATE(exE1y_EM(ttlnmbrnd))
        ALLOCATE(exE1z_EM(ttlnmbrnd))
        ALLOCATE(inE1x_EM(ttlnmbrnd))
        ALLOCATE(inE1y_EM(ttlnmbrnd))
        ALLOCATE(inE1z_EM(ttlnmbrnd))
        ALLOCATE(exE1xdnn_EM(ttlnmbrnd))
        ALLOCATE(exE1ydnn_EM(ttlnmbrnd))
        ALLOCATE(exE1zdnn_EM(ttlnmbrnd))
        ALLOCATE(inE1xdnn_EM(ttlnmbrnd))
        ALLOCATE(inE1ydnn_EM(ttlnmbrnd))
        ALLOCATE(inE1zdnn_EM(ttlnmbrnd))
        ALLOCATE(exE1xdt1_EM(ttlnmbrnd))
        ALLOCATE(exE1ydt1_EM(ttlnmbrnd))
        ALLOCATE(exE1zdt1_EM(ttlnmbrnd))
        ALLOCATE(inE1xdt1_EM(ttlnmbrnd))
        ALLOCATE(inE1ydt1_EM(ttlnmbrnd))
        ALLOCATE(inE1zdt1_EM(ttlnmbrnd))
        ALLOCATE(exE1xdt2_EM(ttlnmbrnd))
        ALLOCATE(exE1ydt2_EM(ttlnmbrnd))
        ALLOCATE(exE1zdt2_EM(ttlnmbrnd))
        ALLOCATE(inE1xdt2_EM(ttlnmbrnd))
        ALLOCATE(inE1ydt2_EM(ttlnmbrnd))
        ALLOCATE(inE1zdt2_EM(ttlnmbrnd))

        ALLOCATE(exE2x_EM(ttlnmbrnd))
        ALLOCATE(exE2y_EM(ttlnmbrnd))
        ALLOCATE(exE2z_EM(ttlnmbrnd))
        ALLOCATE(inE2x_EM(ttlnmbrnd))
        ALLOCATE(inE2y_EM(ttlnmbrnd))
        ALLOCATE(inE2z_EM(ttlnmbrnd))
        ALLOCATE(exE2xdnn_EM(ttlnmbrnd))
        ALLOCATE(exE2ydnn_EM(ttlnmbrnd))
        ALLOCATE(exE2zdnn_EM(ttlnmbrnd))
        ALLOCATE(inE2xdnn_EM(ttlnmbrnd))
        ALLOCATE(inE2ydnn_EM(ttlnmbrnd))
        ALLOCATE(inE2zdnn_EM(ttlnmbrnd))

        ALLOCATE(exE3x_EM(ttlnmbrnd))
        ALLOCATE(exE3y_EM(ttlnmbrnd))
        ALLOCATE(exE3z_EM(ttlnmbrnd))
        ALLOCATE(inE3x_EM(ttlnmbrnd))
        ALLOCATE(inE3y_EM(ttlnmbrnd))
        ALLOCATE(inE3z_EM(ttlnmbrnd))
        ALLOCATE(exE3xdnn_EM(ttlnmbrnd))
        ALLOCATE(exE3ydnn_EM(ttlnmbrnd))
        ALLOCATE(exE3zdnn_EM(ttlnmbrnd))
        ALLOCATE(inE3xdnn_EM(ttlnmbrnd))
        ALLOCATE(inE3ydnn_EM(ttlnmbrnd))
        ALLOCATE(inE3zdnn_EM(ttlnmbrnd))

        ALLOCATE(exH1x_EM(ttlnmbrnd))
        ALLOCATE(exH1y_EM(ttlnmbrnd))
        ALLOCATE(exH1z_EM(ttlnmbrnd))
        ALLOCATE(inH1x_EM(ttlnmbrnd))
        ALLOCATE(inH1y_EM(ttlnmbrnd))
        ALLOCATE(inH1z_EM(ttlnmbrnd))
        ALLOCATE(exH1xdnn_EM(ttlnmbrnd))
        ALLOCATE(exH1ydnn_EM(ttlnmbrnd))
        ALLOCATE(exH1zdnn_EM(ttlnmbrnd))
        ALLOCATE(inH1xdnn_EM(ttlnmbrnd))
        ALLOCATE(inH1ydnn_EM(ttlnmbrnd))
        ALLOCATE(inH1zdnn_EM(ttlnmbrnd))
        ALLOCATE(exH1xdt1_EM(ttlnmbrnd))
        ALLOCATE(exH1ydt1_EM(ttlnmbrnd))
        ALLOCATE(exH1zdt1_EM(ttlnmbrnd))
        ALLOCATE(inH1xdt1_EM(ttlnmbrnd))
        ALLOCATE(inH1ydt1_EM(ttlnmbrnd))
        ALLOCATE(inH1zdt1_EM(ttlnmbrnd))
        ALLOCATE(exH1xdt2_EM(ttlnmbrnd))
        ALLOCATE(exH1ydt2_EM(ttlnmbrnd))
        ALLOCATE(exH1zdt2_EM(ttlnmbrnd))
        ALLOCATE(inH1xdt2_EM(ttlnmbrnd))
        ALLOCATE(inH1ydt2_EM(ttlnmbrnd))
        ALLOCATE(inH1zdt2_EM(ttlnmbrnd))

        ALLOCATE(exH2x_EM(ttlnmbrnd))
        ALLOCATE(exH2y_EM(ttlnmbrnd))
        ALLOCATE(exH2z_EM(ttlnmbrnd))
        ALLOCATE(inH2x_EM(ttlnmbrnd))
        ALLOCATE(inH2y_EM(ttlnmbrnd))
        ALLOCATE(inH2z_EM(ttlnmbrnd))
        ALLOCATE(exH2xdnn_EM(ttlnmbrnd))
        ALLOCATE(exH2ydnn_EM(ttlnmbrnd))
        ALLOCATE(exH2zdnn_EM(ttlnmbrnd))
        ALLOCATE(inH2xdnn_EM(ttlnmbrnd))
        ALLOCATE(inH2ydnn_EM(ttlnmbrnd))
        ALLOCATE(inH2zdnn_EM(ttlnmbrnd))

        ALLOCATE(exH3x_EM(ttlnmbrnd))
        ALLOCATE(exH3y_EM(ttlnmbrnd))
        ALLOCATE(exH3z_EM(ttlnmbrnd))
        ALLOCATE(inH3x_EM(ttlnmbrnd))
        ALLOCATE(inH3y_EM(ttlnmbrnd))
        ALLOCATE(inH3z_EM(ttlnmbrnd))
        ALLOCATE(exH3xdnn_EM(ttlnmbrnd))
        ALLOCATE(exH3ydnn_EM(ttlnmbrnd))
        ALLOCATE(exH3zdnn_EM(ttlnmbrnd))
        ALLOCATE(inH3xdnn_EM(ttlnmbrnd))
        ALLOCATE(inH3ydnn_EM(ttlnmbrnd))
        ALLOCATE(inH3zdnn_EM(ttlnmbrnd))

!$OMP PARALLEL PRIVATE (i)
!$OMP DO
        DO i = 1, ttlnmbrnd

            exE1x_EM(i) = ztpzero
            exE1y_EM(i) = ztpzero
            exE1z_EM(i) = ztpzero
            inE1x_EM(i) = ztpzero
            inE1y_EM(i) = ztpzero
            inE1z_EM(i) = ztpzero
            exE1xdnn_EM(i) = ztpzero
            exE1ydnn_EM(i) = ztpzero
            exE1zdnn_EM(i) = ztpzero
            inE1xdnn_EM(i) = ztpzero
            inE1ydnn_EM(i) = ztpzero
            inE1zdnn_EM(i) = ztpzero
            exE1xdt1_EM(i) = ztpzero
            exE1ydt1_EM(i) = ztpzero
            exE1zdt1_EM(i) = ztpzero
            inE1xdt1_EM(i) = ztpzero
            inE1ydt1_EM(i) = ztpzero
            inE1zdt1_EM(i) = ztpzero
            exE1xdt2_EM(i) = ztpzero
            exE1ydt2_EM(i) = ztpzero
            exE1zdt2_EM(i) = ztpzero
            inE1xdt2_EM(i) = ztpzero
            inE1ydt2_EM(i) = ztpzero
            inE1zdt2_EM(i) = ztpzero

            exE2x_EM(i) = ztpzero
            exE2y_EM(i) = ztpzero
            exE2z_EM(i) = ztpzero
            inE2x_EM(i) = ztpzero
            inE2y_EM(i) = ztpzero
            inE2z_EM(i) = ztpzero
            exE2xdnn_EM(i) = ztpzero
            exE2ydnn_EM(i) = ztpzero
            exE2zdnn_EM(i) = ztpzero
            inE2xdnn_EM(i) = ztpzero
            inE2ydnn_EM(i) = ztpzero
            inE2zdnn_EM(i) = ztpzero

            exE3x_EM(i) = ztpzero
            exE3y_EM(i) = ztpzero
            exE3z_EM(i) = ztpzero
            inE3x_EM(i) = ztpzero
            inE3y_EM(i) = ztpzero
            inE3z_EM(i) = ztpzero
            exE3xdnn_EM(i) = ztpzero
            exE3ydnn_EM(i) = ztpzero
            exE3zdnn_EM(i) = ztpzero
            inE3xdnn_EM(i) = ztpzero
            inE3ydnn_EM(i) = ztpzero
            inE3zdnn_EM(i) = ztpzero

            exH1x_EM(i) = ztpzero
            exH1y_EM(i) = ztpzero
            exH1z_EM(i) = ztpzero
            inH1x_EM(i) = ztpzero
            inH1y_EM(i) = ztpzero
            inH1z_EM(i) = ztpzero
            exH1xdnn_EM(i) = ztpzero
            exH1ydnn_EM(i) = ztpzero
            exH1zdnn_EM(i) = ztpzero
            inH1xdnn_EM(i) = ztpzero
            inH1ydnn_EM(i) = ztpzero
            inH1zdnn_EM(i) = ztpzero
            exH1xdt1_EM(i) = ztpzero
            exH1ydt1_EM(i) = ztpzero
            exH1zdt1_EM(i) = ztpzero
            inH1xdt1_EM(i) = ztpzero
            inH1ydt1_EM(i) = ztpzero
            inH1zdt1_EM(i) = ztpzero
            exH1xdt2_EM(i) = ztpzero
            exH1ydt2_EM(i) = ztpzero
            exH1zdt2_EM(i) = ztpzero
            inH1xdt2_EM(i) = ztpzero
            inH1ydt2_EM(i) = ztpzero
            inH1zdt2_EM(i) = ztpzero

            exH2x_EM(i) = ztpzero
            exH2y_EM(i) = ztpzero
            exH2z_EM(i) = ztpzero
            inH2x_EM(i) = ztpzero
            inH2y_EM(i) = ztpzero
            inH2z_EM(i) = ztpzero
            exH2xdnn_EM(i) = ztpzero
            exH2ydnn_EM(i) = ztpzero
            exH2zdnn_EM(i) = ztpzero
            inH2xdnn_EM(i) = ztpzero
            inH2ydnn_EM(i) = ztpzero
            inH2zdnn_EM(i) = ztpzero

            exH3x_EM(i) = ztpzero
            exH3y_EM(i) = ztpzero
            exH3z_EM(i) = ztpzero
            inH3x_EM(i) = ztpzero
            inH3y_EM(i) = ztpzero
            inH3z_EM(i) = ztpzero
            exH3xdnn_EM(i) = ztpzero
            exH3ydnn_EM(i) = ztpzero
            exH3zdnn_EM(i) = ztpzero
            inH3xdnn_EM(i) = ztpzero
            inH3ydnn_EM(i) = ztpzero
            inH3zdnn_EM(i) = ztpzero

        END DO
!$OMP END DO
!$OMP END PARALLEL

    END SUBROUTINE


END MODULE
