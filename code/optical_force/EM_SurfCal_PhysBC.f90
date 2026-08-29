
! SPDX-FileCopyrightText: 2026 Qiang Sun
! SPDX-License-Identifier: BSD-3-Clause

! Construct incident fields and material boundary data at surface nodes.
! Complex amplitudes use the exp(-i*omega*t) convention; plane-wave spatial
! factors therefore use exp(+i*k dot r). Geometry, wavelength, beam width,
! and penetration depth must be expressed in one consistent length unit.
!
MODULE EM_SurfCal_PhysBC

    USE omp_lib

    USE Pre_Constants
    USE Pre_csvformat

    USE Geom_GlobalData

    USE EM_SurfCal_GlobalData

    IMPLICIT NONE

    COMPLEX(KIND=KIND(1.0D0)), PRIVATE :: poralz_c1 = ztpone, poralz_c2 = ztponei
    DOUBLE PRECISION, PRIVATE :: fdm_dlh = 1.d-3
    INTEGER, PRIVATE :: fdm_nmx = 7

    CONTAINS

    ! Populate incident E/H values and their normal derivatives at surface
    ! nodes, then initialise material-side boundary data for the solver.
    ! Input_Source_EM.dat is read only when impressed sources are declared.

    SUBROUTINE GetPhysBC_EM

        IMPLICIT NONE

        INTEGER :: i, j, k, ithprtl, ii, jj, kk, id_tp, IOS

        DOUBLE PRECISION :: uMedium = 1.0d0, uPrism = 1.0d0
        ! For evanescent-wave excitation, incFeature_EM is the 1/e intensity
        ! penetration depth; field amplitude decays as exp[-x/(2*depth)].
        DOUBLE PRECISION :: penDepthDef
        DOUBLE PRECISION :: betaDef, ampE0, thetaI

        COMPLEX(KIND=KIND(1.0D0)) :: thetaT, sinT, cosT, ampE_P, ampE_S

        DOUBLE PRECISION :: tp, tp1, tp2, tp3, tp4, tp5, tp6, tp_eva, &
        &                   tpomega_0, tpomega_z, tpzR, tpRRz, tpr2, &
        &                   tpx, tpy, tpz, tpxhat, tpyhat, tpzhat, tpr, &
        &                   rhat_dot_p, rhat_crs_p_x, rhat_crs_p_y, rhat_crs_p_z, &
        &                   rhat_crs_p_crs_rhad_x, rhat_crs_p_crs_rhad_y, rhat_crs_p_crs_rhad_z

        DOUBLE PRECISION :: tpEin_read(27)

        COMPLEX(KIND=KIND(1.0D0)) :: ztp,ztp1,ztp2,ztp3,ztp4,ztp5,ztp6,ztp7,ztp8,ztp9
        COMPLEX(KIND=KIND(1.0D0)) :: ztp_eva, ztp1x,ztp1y,ztp1z,ztp2x,ztp2y,ztp2z
        COMPLEX(KIND=KIND(1.0D0)) :: ztpeps,ztpmiu,ztpk,ztpexp

        DOUBLE PRECISION :: dlh=1.d-7,dlh_1,dmx,dmy,dmz
        COMPLEX(KIND=KIND(1.0D0)) :: ztpdExdx,ztpdExdy,ztpdExdz,ztpdEydx,ztpdEydy,ztpdEydz,&
        &             ztpdEzdx,ztpdEzdy,ztpdEzdz


    ! Plane-wave excitation.

        IF (excitetype_EM == 'spe') THEN

            DO i = 1, ttlnmbrnd

                ithprtl = 1
                IF (nmbrprtl > 1 .AND. i > ndendID(1)) THEN
                    DO id_tp = 2, nmbrprtl
                        IF (i >= ndstaID(id_tp) .AND. i <= ndendID(id_tp)) THEN
                            ithprtl = id_tp
                            EXIT
                        END IF
                    END DO
                END IF

                IF (corelnkshell(ithprtl) == 0) THEN
                    tp1 = inckx_EM*xnd(i)&
                    &    +incky_EM*ynd(i)&
                    &    +inckz_EM*znd(i)
                    ztp1x = incFieldx_EM*CDEXP( ztponei*exk_EM*tp1 )
                    ztp1y = incFieldy_EM*CDEXP( ztponei*exk_EM*tp1 )
                    ztp1z = incFieldz_EM*CDEXP( ztponei*exk_EM*tp1 )
                    exE1x_EM(i) = exE1x_EM(i) + ztp1x
                    exE1y_EM(i) = exE1y_EM(i) + ztp1y
                    exE1z_EM(i) = exE1z_EM(i) + ztp1z
                    tp2 = inckx_EM*nnx(i)&
                    &    +incky_EM*nny(i)&
                    &    +inckz_EM*nnz(i)
                    exE1xdnn_EM(i) = exE1xdnn_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exE1ydnn_EM(i) = exE1ydnn_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exE1zdnn_EM(i) = exE1zdnn_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                    tp2 = inckx_EM*t1x(i)&
                    &    +incky_EM*t1y(i)&
                    &    +inckz_EM*t1z(i)
                    exE1xdt1_EM(i) = exE1xdt1_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exE1ydt1_EM(i) = exE1ydt1_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exE1zdt1_EM(i) = exE1zdt1_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                    tp2 = inckx_EM*t2x(i)&
                    &    +incky_EM*t2y(i)&
                    &    +inckz_EM*t2z(i)
                    exE1xdt2_EM(i) = exE1xdt2_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exE1ydt2_EM(i) = exE1ydt2_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exE1zdt2_EM(i) = exE1zdt2_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                    ztp1 = ztp1x*(ztponei*exk_EM*inckx_EM)
                    ztp2 = ztp1x*(ztponei*exk_EM*incky_EM)
                    ztp3 = ztp1x*(ztponei*exk_EM*inckz_EM)
                    ztp4 = ztp1y*(ztponei*exk_EM*inckx_EM)
                    ztp5 = ztp1y*(ztponei*exk_EM*incky_EM)
                    ztp6 = ztp1y*(ztponei*exk_EM*inckz_EM)
                    ztp7 = ztp1z*(ztponei*exk_EM*inckx_EM)
                    ztp8 = ztp1z*(ztponei*exk_EM*incky_EM)
                    ztp9 = ztp1z*(ztponei*exk_EM*inckz_EM)

                    ztp1x = ztpzero
                    ztp1y = ztpzero
                    ztp1z = ztpzero
                    exH1x_EM(i) = exH1x_EM(i) + ztp1x
                    exH1y_EM(i) = exH1y_EM(i) + ztp1y
                    exH1z_EM(i) = exH1z_EM(i) + ztp1z
                    tp2 = inckx_EM*nnx(i)&
                    &    +incky_EM*nny(i)&
                    &    +inckz_EM*nnz(i)
                    exH1xdnn_EM(i) = exH1xdnn_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exH1ydnn_EM(i) = exH1ydnn_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exH1zdnn_EM(i) = exH1zdnn_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                    tp2 = inckx_EM*t1x(i)&
                    &    +incky_EM*t1y(i)&
                    &    +inckz_EM*t1z(i)
                    exH1xdt1_EM(i) = exH1xdt1_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exH1ydt1_EM(i) = exH1ydt1_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exH1zdt1_EM(i) = exH1zdt1_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                    tp2 = inckx_EM*t2x(i)&
                    &    +incky_EM*t2y(i)&
                    &    +inckz_EM*t2z(i)
                    exH1xdt2_EM(i) = exH1xdt2_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exH1ydt2_EM(i) = exH1ydt2_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exH1zdt2_EM(i) = exH1zdt2_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                END IF

            END DO

        END IF

        IF (excitetype_EM == 'sph') THEN

            DO i = 1, ttlnmbrnd

                ithprtl = 1
                IF (nmbrprtl > 1 .AND. i > ndendID(1)) THEN
                    DO id_tp = 2, nmbrprtl
                        IF (i >= ndstaID(id_tp) .AND. i <= ndendID(id_tp)) THEN
                            ithprtl = id_tp
                            EXIT
                        END IF
                    END DO
                END IF

                IF (corelnkshell(ithprtl) == 0) THEN
                    tp1 = inckx_EM*xnd(i)&
                    &    +incky_EM*ynd(i)&
                    &    +inckz_EM*znd(i)
                    ztp1x = incFieldx_EM*CDEXP( ztponei*exk_EM*tp1 )
                    ztp1y = incFieldy_EM*CDEXP( ztponei*exk_EM*tp1 )
                    ztp1z = incFieldz_EM*CDEXP( ztponei*exk_EM*tp1 )
                    exH1x_EM(i) = exH1x_EM(i) + ztp1x
                    exH1y_EM(i) = exH1y_EM(i) + ztp1y
                    exH1z_EM(i) = exH1z_EM(i) + ztp1z
                    tp2 = inckx_EM*nnx(i)&
                    &    +incky_EM*nny(i)&
                    &    +inckz_EM*nnz(i)
                    exH1xdnn_EM(i) = exH1xdnn_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exH1ydnn_EM(i) = exH1ydnn_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exH1zdnn_EM(i) = exH1zdnn_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                    tp2 = inckx_EM*t1x(i)&
                    &    +incky_EM*t1y(i)&
                    &    +inckz_EM*t1z(i)
                    exH1xdt1_EM(i) = exH1xdt1_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exH1ydt1_EM(i) = exH1ydt1_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exH1zdt1_EM(i) = exH1zdt1_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                    tp2 = inckx_EM*t2x(i)&
                    &    +incky_EM*t2y(i)&
                    &    +inckz_EM*t2z(i)
                    exH1xdt2_EM(i) = exH1xdt2_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exH1ydt2_EM(i) = exH1ydt2_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exH1zdt2_EM(i) = exH1zdt2_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                    ztp1 = ztp1x*(ztponei*exk_EM*inckx_EM)
                    ztp2 = ztp1x*(ztponei*exk_EM*incky_EM)
                    ztp3 = ztp1x*(ztponei*exk_EM*inckz_EM)
                    ztp4 = ztp1y*(ztponei*exk_EM*inckx_EM)
                    ztp5 = ztp1y*(ztponei*exk_EM*incky_EM)
                    ztp6 = ztp1y*(ztponei*exk_EM*inckz_EM)
                    ztp7 = ztp1z*(ztponei*exk_EM*inckx_EM)
                    ztp8 = ztp1z*(ztponei*exk_EM*incky_EM)
                    ztp9 = ztp1z*(ztponei*exk_EM*inckz_EM)

                    ztp1x = ztpzero
                    ztp1y = ztpzero
                    ztp1z = ztpzero
                    exE1x_EM(i) = exE1x_EM(i) + ztp1x
                    exE1y_EM(i) = exE1y_EM(i) + ztp1y
                    exE1z_EM(i) = exE1z_EM(i) + ztp1z
                    tp2 = inckx_EM*nnx(i)&
                    &    +incky_EM*nny(i)&
                    &    +inckz_EM*nnz(i)
                    exE1xdnn_EM(i) = exE1xdnn_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exE1ydnn_EM(i) = exE1ydnn_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exE1zdnn_EM(i) = exE1zdnn_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                    tp2 = inckx_EM*t1x(i)&
                    &    +incky_EM*t1y(i)&
                    &    +inckz_EM*t1z(i)
                    exE1xdt1_EM(i) = exE1xdt1_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exE1ydt1_EM(i) = exE1ydt1_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exE1zdt1_EM(i) = exE1zdt1_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                    tp2 = inckx_EM*t2x(i)&
                    &    +incky_EM*t2y(i)&
                    &    +inckz_EM*t2z(i)
                    exE1xdt2_EM(i) = exE1xdt2_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exE1ydt2_EM(i) = exE1ydt2_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exE1zdt2_EM(i) = exE1zdt2_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                END IF

            END DO

        END IF

        IF (excitetype_EM == 'pwe') THEN

            DO i = 1, ttlnmbrnd

                ithprtl = 1
                IF (nmbrprtl > 1 .AND. i > ndendID(1)) THEN
                    DO id_tp = 2, nmbrprtl
                        IF (i >= ndstaID(id_tp) .AND. i <= ndendID(id_tp)) THEN
                            ithprtl = id_tp
                            EXIT
                        END IF
                    END DO
                END IF

                IF (corelnkshell(ithprtl) == 0) THEN
                    tp1 = inckx_EM*xnd(i)&
                    &    +incky_EM*ynd(i)&
                    &    +inckz_EM*znd(i)
                    IF (incOrder_EM == 1) THEN
                        ztp1x = incFieldx_EM*CDEXP( ztponei*exk_EM*tp1 )
                        ztp1y = incFieldy_EM*CDEXP( ztponei*exk_EM*tp1 + ztponei*incFeature_EM)
                        ztp1z = incFieldz_EM*CDEXP( ztponei*exk_EM*tp1 )
                    END IF
                    IF (incOrder_EM == 2) THEN
                        ztp1x = incFieldx_EM*CDEXP( ztponei*exk_EM*tp1 )
                        ztp1y = incFieldy_EM*CDEXP( ztponei*exk_EM*tp1 )
                        ztp1z = incFieldz_EM*CDEXP( ztponei*exk_EM*tp1 + ztponei*incFeature_EM )
                    END IF
                    IF (incOrder_EM == 3) THEN
                        ztp1x = incFieldx_EM*CDEXP( ztponei*exk_EM*tp1 + ztponei*incFeature_EM )
                        ztp1y = incFieldy_EM*CDEXP( ztponei*exk_EM*tp1 )
                        ztp1z = incFieldz_EM*CDEXP( ztponei*exk_EM*tp1 )
                    END IF
                    exE1x_EM(i) = exE1x_EM(i) + ztp1x
                    exE1y_EM(i) = exE1y_EM(i) + ztp1y
                    exE1z_EM(i) = exE1z_EM(i) + ztp1z
                    tp2 = inckx_EM*nnx(i)&
                    &    +incky_EM*nny(i)&
                    &    +inckz_EM*nnz(i)
                    exE1xdnn_EM(i) = exE1xdnn_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exE1ydnn_EM(i) = exE1ydnn_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exE1zdnn_EM(i) = exE1zdnn_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                    tp2 = inckx_EM*t1x(i)&
                    &    +incky_EM*t1y(i)&
                    &    +inckz_EM*t1z(i)
                    exE1xdt1_EM(i) = exE1xdt1_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exE1ydt1_EM(i) = exE1ydt1_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exE1zdt1_EM(i) = exE1zdt1_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                    tp2 = inckx_EM*t2x(i)&
                    &    +incky_EM*t2y(i)&
                    &    +inckz_EM*t2z(i)
                    exE1xdt2_EM(i) = exE1xdt2_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exE1ydt2_EM(i) = exE1ydt2_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exE1zdt2_EM(i) = exE1zdt2_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                    ztp1 = ztp1x*(ztponei*exk_EM*inckx_EM)
                    ztp2 = ztp1x*(ztponei*exk_EM*incky_EM)
                    ztp3 = ztp1x*(ztponei*exk_EM*inckz_EM)
                    ztp4 = ztp1y*(ztponei*exk_EM*inckx_EM)
                    ztp5 = ztp1y*(ztponei*exk_EM*incky_EM)
                    ztp6 = ztp1y*(ztponei*exk_EM*inckz_EM)
                    ztp7 = ztp1z*(ztponei*exk_EM*inckx_EM)
                    ztp8 = ztp1z*(ztponei*exk_EM*incky_EM)
                    ztp9 = ztp1z*(ztponei*exk_EM*inckz_EM)

                    ztp = 1.0d0/(ztponei*vcm_mu0*exmiu_EM*AngFrqnc_EM)
                    ztp1x = ztp*(ztp8 - ztp6)
                    ztp1y = ztp*(ztp3 - ztp7)
                    ztp1z = ztp*(ztp4 - ztp2)
                    exH1x_EM(i) = exH1x_EM(i) + ztp1x
                    exH1y_EM(i) = exH1y_EM(i) + ztp1y
                    exH1z_EM(i) = exH1z_EM(i) + ztp1z
                    tp2 = inckx_EM*nnx(i)&
                    &    +incky_EM*nny(i)&
                    &    +inckz_EM*nnz(i)
                    exH1xdnn_EM(i) = exH1xdnn_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exH1ydnn_EM(i) = exH1ydnn_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exH1zdnn_EM(i) = exH1zdnn_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                    tp2 = inckx_EM*t1x(i)&
                    &    +incky_EM*t1y(i)&
                    &    +inckz_EM*t1z(i)
                    exH1xdt1_EM(i) = exH1xdt1_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exH1ydt1_EM(i) = exH1ydt1_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exH1zdt1_EM(i) = exH1zdt1_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                    tp2 = inckx_EM*t2x(i)&
                    &    +incky_EM*t2y(i)&
                    &    +inckz_EM*t2z(i)
                    exH1xdt2_EM(i) = exH1xdt2_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exH1ydt2_EM(i) = exH1ydt2_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exH1zdt2_EM(i) = exH1zdt2_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                END IF

            END DO

        END IF

        IF (excitetype_EM == 'pwh') THEN

            DO i = 1, ttlnmbrnd

                ithprtl = 1
                IF (nmbrprtl > 1 .AND. i > ndendID(1)) THEN
                    DO id_tp = 2, nmbrprtl
                        IF (i >= ndstaID(id_tp) .AND. i <= ndendID(id_tp)) THEN
                            ithprtl = id_tp
                            EXIT
                        END IF
                    END DO
                END IF

                IF (corelnkshell(ithprtl) == 0) THEN
                    tp1 = inckx_EM*xnd(i)&
                    &    +incky_EM*ynd(i)&
                    &    +inckz_EM*znd(i)
                    IF (incOrder_EM == 1) THEN
                        ztp1x = incFieldx_EM*CDEXP( ztponei*exk_EM*tp1 )
                        ztp1y = incFieldy_EM*CDEXP( ztponei*exk_EM*tp1 + ztponei*incFeature_EM)
                        ztp1z = incFieldz_EM*CDEXP( ztponei*exk_EM*tp1 )
                    END IF
                    IF (incOrder_EM == 2) THEN
                        ztp1x = incFieldx_EM*CDEXP( ztponei*exk_EM*tp1 )
                        ztp1y = incFieldy_EM*CDEXP( ztponei*exk_EM*tp1 )
                        ztp1z = incFieldz_EM*CDEXP( ztponei*exk_EM*tp1 + ztponei*incFeature_EM )
                    END IF
                    IF (incOrder_EM == 3) THEN
                        ztp1x = incFieldx_EM*CDEXP( ztponei*exk_EM*tp1 + ztponei*incFeature_EM )
                        ztp1y = incFieldy_EM*CDEXP( ztponei*exk_EM*tp1 )
                        ztp1z = incFieldz_EM*CDEXP( ztponei*exk_EM*tp1 )
                    END IF
                    exH1x_EM(i) = exH1x_EM(i) + ztp1x
                    exH1y_EM(i) = exH1y_EM(i) + ztp1y
                    exH1z_EM(i) = exH1z_EM(i) + ztp1z
                    tp2 = inckx_EM*nnx(i)&
                    &    +incky_EM*nny(i)&
                    &    +inckz_EM*nnz(i)
                    exH1xdnn_EM(i) = exH1xdnn_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exH1ydnn_EM(i) = exH1ydnn_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exH1zdnn_EM(i) = exH1zdnn_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                    tp2 = inckx_EM*t1x(i)&
                    &    +incky_EM*t1y(i)&
                    &    +inckz_EM*t1z(i)
                    exH1xdt1_EM(i) = exH1xdt1_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exH1ydt1_EM(i) = exH1ydt1_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exH1zdt1_EM(i) = exH1zdt1_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                    tp2 = inckx_EM*t2x(i)&
                    &    +incky_EM*t2y(i)&
                    &    +inckz_EM*t2z(i)
                    exH1xdt2_EM(i) = exH1xdt2_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exH1ydt2_EM(i) = exH1ydt2_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exH1zdt2_EM(i) = exH1zdt2_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                    ztp1 = ztp1x*(ztponei*exk_EM*inckx_EM)
                    ztp2 = ztp1x*(ztponei*exk_EM*incky_EM)
                    ztp3 = ztp1x*(ztponei*exk_EM*inckz_EM)
                    ztp4 = ztp1y*(ztponei*exk_EM*inckx_EM)
                    ztp5 = ztp1y*(ztponei*exk_EM*incky_EM)
                    ztp6 = ztp1y*(ztponei*exk_EM*inckz_EM)
                    ztp7 = ztp1z*(ztponei*exk_EM*inckx_EM)
                    ztp8 = ztp1z*(ztponei*exk_EM*incky_EM)
                    ztp9 = ztp1z*(ztponei*exk_EM*inckz_EM)

                    ztp =-1.0d0/(ztponei*vcm_eps0*exeps_EM*AngFrqnc_EM)
                    ztp1x = ztp*(ztp8 - ztp6)
                    ztp1y = ztp*(ztp3 - ztp7)
                    ztp1z = ztp*(ztp4 - ztp2)
                    exE1x_EM(i) = exE1x_EM(i) + ztp1x
                    exE1y_EM(i) = exE1y_EM(i) + ztp1y
                    exE1z_EM(i) = exE1z_EM(i) + ztp1z
                    tp2 = inckx_EM*nnx(i)&
                    &    +incky_EM*nny(i)&
                    &    +inckz_EM*nnz(i)
                    exE1xdnn_EM(i) = exE1xdnn_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exE1ydnn_EM(i) = exE1ydnn_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exE1zdnn_EM(i) = exE1zdnn_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                    tp2 = inckx_EM*t1x(i)&
                    &    +incky_EM*t1y(i)&
                    &    +inckz_EM*t1z(i)
                    exE1xdt1_EM(i) = exE1xdt1_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exE1ydt1_EM(i) = exE1ydt1_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exE1zdt1_EM(i) = exE1zdt1_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                    tp2 = inckx_EM*t2x(i)&
                    &    +incky_EM*t2y(i)&
                    &    +inckz_EM*t2z(i)
                    exE1xdt2_EM(i) = exE1xdt2_EM(i) &
                    &               +ztp1x*(ztponei*exk_EM*tp2)
                    exE1ydt2_EM(i) = exE1ydt2_EM(i) &
                    &               +ztp1y*(ztponei*exk_EM*tp2)
                    exE1zdt2_EM(i) = exE1zdt2_EM(i) &
                    &               +ztp1z*(ztponei*exk_EM*tp2)
                END IF

            END DO

        END IF



    ! Evanescent-wave excitation.
        IF (excitetype_EM == 'eva') THEN

            penDepthDef = incFeature_EM

            ! betaDef is the field-amplitude decay constant implied by the
            ! 1/e intensity penetration depth.
            betaDef = 1.0d0/(2.0d0*penDepthDef)
            tp = REAL(exk_EM)      !kM <------- tp
            thetaI = ASIN(DSQRT((((betaDef/tp)**2) + 1)*(exepsn_EM**2)/(prismn_EM**2)))

            tp1 = DSQRT((1.0d0 / (2.0d0 * penDepthDef * tp))**2 + 1.0d0)    !cosh(alfa)
            tp2 = DSQRT((prismn_EM * DSIN(thetaI) / exepsn_EM)**2 - 1.0d0)       !sinh(alfa)

            tp3 = exepsn_EM/prismn_EM

            tp4 = (DCMPLX(exepsn_EM,exepsk_EM)**2)/(DCMPLX(prismn_EM,prismk_EM)**2)   !eps/eps_1
            ampE_P = 2.0d0 * tp3 * DCOS(thetaI) / (tp4 * DCOS(thetaI)+ ztponei * tp3 * tp2)

            tp4 = 1.0d0     !mu/mu_1
            ampE_S = 2.0d0 * tp4 * DCOS(thetaI) / (tp4 * DCOS(thetaI)+ ztponei * tp3 * tp2)


            IF (poltype_EM == 'S') THEN

                tp_eva =-REAL(exk_EM)*tp2
                ztp_eva= ztponei*REAL(exk_EM)*tp1

                DO i = 1, ttlnmbrnd

                    ithprtl = 1
                    IF (nmbrprtl > 1 .AND. i > ndendID(1)) THEN
                        DO id_tp = 2, nmbrprtl
                            IF (i >= ndstaID(id_tp) .AND. i <= ndendID(id_tp)) THEN
                                ithprtl = id_tp
                                EXIT
                            END IF
                        END DO
                    END IF

                    IF (corelnkshell(ithprtl) == 0) THEN

                        ztp1x = ztpzero
                        ztp1y = ampE_S*DEXP(tp_eva*xnd(i))*CDEXP(ztp_eva*znd(i))
                        ztp1z = ztpzero
                        exE1x_EM(i) = exE1x_EM(i) + ztp1x
                        exE1y_EM(i) = exE1y_EM(i) + ztp1y
                        exE1z_EM(i) = exE1z_EM(i) + ztp1z
                        exE1xdnn_EM(i) = exE1xdnn_EM(i) &
                        &               +ztp1x*(tp_eva*nnx(i) + ztp_eva*nnz(i))
                        exE1ydnn_EM(i) = exE1ydnn_EM(i) &
                        &               +ztp1y*(tp_eva*nnx(i) + ztp_eva*nnz(i))
                        exE1zdnn_EM(i) = exE1zdnn_EM(i) &
                        &               +ztp1z*(tp_eva*nnx(i) + ztp_eva*nnz(i))
                        exE1xdt1_EM(i) = exE1xdt1_EM(i) &
                        &               +ztp1x*(tp_eva*t1x(i) + ztp_eva*t1z(i))
                        exE1ydt1_EM(i) = exE1ydt1_EM(i) &
                        &               +ztp1y*(tp_eva*t1x(i) + ztp_eva*t1z(i))
                        exE1zdt1_EM(i) = exE1zdt1_EM(i) &
                        &               +ztp1z*(tp_eva*t1x(i) + ztp_eva*t1z(i))
                        exE1xdt2_EM(i) = exE1xdt2_EM(i) &
                        &               +ztp1x*(tp_eva*t2x(i) + ztp_eva*t2z(i))
                        exE1ydt2_EM(i) = exE1ydt2_EM(i) &
                        &               +ztp1y*(tp_eva*t2x(i) + ztp_eva*t2z(i))
                        exE1zdt2_EM(i) = exE1zdt2_EM(i) &
                        &               +ztp1z*(tp_eva*t2x(i) + ztp_eva*t2z(i))

                        ztp1 = ztpzero
                        ztp2 = ztpzero
                        ztp3 = ztpzero
                        ztp4 = ztp1y * tp_eva
                        ztp5 = ztpzero
                        ztp6 = ztp1y * ztp_eva
                        ztp7 = ztpzero
                        ztp8 = ztpzero
                        ztp9 = ztpzero

                        ztp = 1.0d0/(ztponei*vcm_mu0*exmiu_EM*AngFrqnc_EM)
                        ztp1x = ztp*(ztp8 - ztp6)
                        ztp1y = ztp*(ztp3 - ztp7)
                        ztp1z = ztp*(ztp4 - ztp2)
                        exH1x_EM(i) = exH1x_EM(i) + ztp1x
                        exH1y_EM(i) = exH1y_EM(i) + ztp1y
                        exH1z_EM(i) = exH1z_EM(i) + ztp1z
                        exH1xdnn_EM(i) = exH1xdnn_EM(i) &
                        &               +ztp1x*(tp_eva*nnx(i) + ztp_eva*nnz(i))
                        exH1ydnn_EM(i) = exH1ydnn_EM(i) &
                        &               +ztp1y*(tp_eva*nnx(i) + ztp_eva*nnz(i))
                        exH1zdnn_EM(i) = exH1zdnn_EM(i) &
                        &               +ztp1z*(tp_eva*nnx(i) + ztp_eva*nnz(i))
                        exH1xdt1_EM(i) = exH1xdt1_EM(i) &
                        &               +ztp1x*(tp_eva*t1x(i) + ztp_eva*t1z(i))
                        exH1ydt1_EM(i) = exH1ydt1_EM(i) &
                        &               +ztp1y*(tp_eva*t1x(i) + ztp_eva*t1z(i))
                        exH1zdt1_EM(i) = exH1zdt1_EM(i) &
                        &               +ztp1z*(tp_eva*t1x(i) + ztp_eva*t1z(i))
                        exH1xdt2_EM(i) = exH1xdt2_EM(i) &
                        &               +ztp1x*(tp_eva*t2x(i) + ztp_eva*t2z(i))
                        exH1ydt2_EM(i) = exH1ydt2_EM(i) &
                        &               +ztp1y*(tp_eva*t2x(i) + ztp_eva*t2z(i))
                        exH1zdt2_EM(i) = exH1zdt2_EM(i) &
                        &               +ztp1z*(tp_eva*t2x(i) + ztp_eva*t2z(i))

                    END IF

                END DO

            ELSE IF (poltype_EM == 'P') THEN !% <----- we're interested in this one!

                sinT = tp1
                cosT = ztponei*tp2

                tp_eva =-REAL(exk_EM)*tp2
                ztp_eva= ztponei*REAL(exk_EM)*tp1

                DO i = 1, ttlnmbrnd

                    ithprtl = 1
                    IF (nmbrprtl > 1 .AND. i > ndendID(1)) THEN
                        DO id_tp = 2, nmbrprtl
                            IF (i >= ndstaID(id_tp) .AND. i <= ndendID(id_tp)) THEN
                                ithprtl = id_tp
                                EXIT
                            END IF
                        END DO
                    END IF

                    IF (corelnkshell(ithprtl) == 0) THEN

                        ztp1x = ampE_P*sinT*DEXP(tp_eva*xnd(i))*CDEXP(ztp_eva*znd(i))
                        ztp1y = ztpzero
                        ztp1z =-ampE_P*cosT*DEXP(tp_eva*xnd(i))*CDEXP(ztp_eva*znd(i))
                        exE1x_EM(i) = exE1x_EM(i) + ztp1x
                        exE1y_EM(i) = exE1y_EM(i) + ztp1y
                        exE1z_EM(i) = exE1z_EM(i) + ztp1z
                        exE1xdnn_EM(i) = exE1xdnn_EM(i) &
                        &               +ztp1x*(tp_eva*nnx(i) + ztp_eva*nnz(i))
                        exE1ydnn_EM(i) = exE1ydnn_EM(i) &
                        &               +ztp1y*(tp_eva*nnx(i) + ztp_eva*nnz(i))
                        exE1zdnn_EM(i) = exE1zdnn_EM(i) &
                        &               +ztp1z*(tp_eva*nnx(i) + ztp_eva*nnz(i))
                        exE1xdt1_EM(i) = exE1xdt1_EM(i) &
                        &               +ztp1x*(tp_eva*t1x(i) + ztp_eva*t1z(i))
                        exE1ydt1_EM(i) = exE1ydt1_EM(i) &
                        &               +ztp1y*(tp_eva*t1x(i) + ztp_eva*t1z(i))
                        exE1zdt1_EM(i) = exE1zdt1_EM(i) &
                        &               +ztp1z*(tp_eva*t1x(i) + ztp_eva*t1z(i))
                        exE1xdt2_EM(i) = exE1xdt2_EM(i) &
                        &               +ztp1x*(tp_eva*t2x(i) + ztp_eva*t2z(i))
                        exE1ydt2_EM(i) = exE1ydt2_EM(i) &
                        &               +ztp1y*(tp_eva*t2x(i) + ztp_eva*t2z(i))
                        exE1zdt2_EM(i) = exE1zdt2_EM(i) &
                        &               +ztp1z*(tp_eva*t2x(i) + ztp_eva*t2z(i))

                        ztp1 = ztp1x * tp_eva
                        ztp2 = ztpzero
                        ztp3 = ztp1x * ztp_eva
                        ztp4 = ztpzero
                        ztp5 = ztpzero
                        ztp6 = ztpzero
                        ztp7 = ztp1z * tp_eva
                        ztp8 = ztpzero
                        ztp9 = ztp1z * ztp_eva

                        ztp = 1.0d0/(ztponei*vcm_mu0*exmiu_EM*AngFrqnc_EM)
                        ztp1x = ztp*(ztp8 - ztp6)
                        ztp1y = ztp*(ztp3 - ztp7)
                        ztp1z = ztp*(ztp4 - ztp2)
                        exH1x_EM(i) = exH1x_EM(i) + ztp1x
                        exH1y_EM(i) = exH1y_EM(i) + ztp1y
                        exH1z_EM(i) = exH1z_EM(i) + ztp1z
                        exH1xdnn_EM(i) = exH1xdnn_EM(i) &
                        &               +ztp1x*(tp_eva*nnx(i) + ztp_eva*nnz(i))
                        exH1ydnn_EM(i) = exH1ydnn_EM(i) &
                        &               +ztp1y*(tp_eva*nnx(i) + ztp_eva*nnz(i))
                        exH1zdnn_EM(i) = exH1zdnn_EM(i) &
                        &               +ztp1z*(tp_eva*nnx(i) + ztp_eva*nnz(i))
                        exH1xdt1_EM(i) = exH1xdt1_EM(i) &
                        &               +ztp1x*(tp_eva*t1x(i) + ztp_eva*t1z(i))
                        exH1ydt1_EM(i) = exH1ydt1_EM(i) &
                        &               +ztp1y*(tp_eva*t1x(i) + ztp_eva*t1z(i))
                        exH1zdt1_EM(i) = exH1zdt1_EM(i) &
                        &               +ztp1z*(tp_eva*t1x(i) + ztp_eva*t1z(i))
                        exH1xdt2_EM(i) = exH1xdt2_EM(i) &
                        &               +ztp1x*(tp_eva*t2x(i) + ztp_eva*t2z(i))
                        exH1ydt2_EM(i) = exH1ydt2_EM(i) &
                        &               +ztp1y*(tp_eva*t2x(i) + ztp_eva*t2z(i))
                        exH1zdt2_EM(i) = exH1zdt2_EM(i) &
                        &               +ztp1z*(tp_eva*t2x(i) + ztp_eva*t2z(i))

                    END IF

                END DO

            END IF

        END IF



    !Gaussian Beam

        IF (excitetype_EM == 'gau') THEN

            DO i = 1, ttlnmbrnd

                ithprtl = 1
                IF (nmbrprtl > 1 .AND. i > ndendID(1)) THEN
                    DO id_tp = 2, nmbrprtl
                        IF (i >= ndstaID(id_tp) .AND. i <= ndendID(id_tp)) THEN
                            ithprtl = id_tp
                            EXIT
                        END IF
                    END DO
                END IF

                IF (corelnkshell(ithprtl) == 0) THEN

                    dmx = xnd(i)
                    dmy = ynd(i)
                    dmz = znd(i)
                    dlh_1 = 1.0d0/dlh

                    CALL GauBM_Einc(i,ztp1,ztp2,ztp3)

                    exE1x_EM(i) = exE1x_EM(i) + ztp1
                    exE1y_EM(i) = exE1y_EM(i) + ztp2
                    exE1z_EM(i) = exE1z_EM(i) + ztp3

!                   CALL FDMGauBM_dEdninc(dmx,dmy,dmz,&
                    CALL GauBM_dEdninc(i,ztp1,ztp2,ztp3,ztp4,ztp5,ztp6,ztp7,ztp8,ztp9)

                    exE1xdnn_EM(i) = exE1xdnn_EM(i) + ztp1
                    exE1ydnn_EM(i) = exE1ydnn_EM(i) + ztp2
                    exE1zdnn_EM(i) = exE1zdnn_EM(i) + ztp3
                    exE1xdt1_EM(i) = exE1xdt1_EM(i) + ztp4
                    exE1ydt1_EM(i) = exE1ydt1_EM(i) + ztp5
                    exE1zdt1_EM(i) = exE1zdt1_EM(i) + ztp6
                    exE1xdt2_EM(i) = exE1xdt2_EM(i) + ztp7
                    exE1ydt2_EM(i) = exE1ydt2_EM(i) + ztp8
                    exE1zdt2_EM(i) = exE1zdt2_EM(i) + ztp9

                    CALL GauBM_Hinc(i,ztp1,ztp2,ztp3)

                    exH1x_EM(i) = exH1x_EM(i) + ztp1
                    exH1y_EM(i) = exH1y_EM(i) + ztp2
                    exH1z_EM(i) = exH1z_EM(i) + ztp3

!                   CALL FDMGauBM_dHdninc(dmx,dmy,dmz,&
                    CALL GauBM_dHdninc(i,ztp1,ztp2,ztp3,ztp4,ztp5,ztp6,ztp7,ztp8,ztp9)

                    exH1xdnn_EM(i) = exH1xdnn_EM(i) + ztp1
                    exH1ydnn_EM(i) = exH1ydnn_EM(i) + ztp2
                    exH1zdnn_EM(i) = exH1zdnn_EM(i) + ztp3
                    exH1xdt1_EM(i) = exH1xdt1_EM(i) + ztp4
                    exH1ydt1_EM(i) = exH1ydt1_EM(i) + ztp5
                    exH1zdt1_EM(i) = exH1zdt1_EM(i) + ztp6
                    exH1xdt2_EM(i) = exH1xdt2_EM(i) + ztp7
                    exH1ydt2_EM(i) = exH1ydt2_EM(i) + ztp8
                    exH1zdt2_EM(i) = exH1zdt2_EM(i) + ztp9

                END IF

            END DO

        END IF



    !Gaussian Beam 5th order

        IF (excitetype_EM == 'gb5') THEN

            IF (poltype_EM == 'e') THEN

                DO i = 1, ttlnmbrnd

                    ithprtl = 1
                    IF (nmbrprtl > 1 .AND. i > ndendID(1)) THEN
                        DO id_tp = 2, nmbrprtl
                            IF (i >= ndstaID(id_tp) .AND. i <= ndendID(id_tp)) THEN
                                ithprtl = id_tp
                                EXIT
                            END IF
                        END DO
                    END IF

                    IF (corelnkshell(ithprtl) == 0) THEN

                        dmx = xnd(i)
                        dmy = ynd(i)
                        dmz = znd(i)
                        dlh_1 = 1.0d0/dlh

                        tp = DBLE(incOrder_EM)/180.0d0*pai
                        ztp = exk_EM/(vcm_mu0*exmiu_EM*AngFrqnc_EM)

                        CALL GauBM5th_Einc(i,ztp1,ztp2,ztp3)    !TEM^(x)_00~E

                        exE1x_EM(i) = exE1x_EM(i) + ztp1 * DCOS(tp)
                        exE1y_EM(i) = exE1y_EM(i) + ztp2 * DCOS(tp)
                        exE1z_EM(i) = exE1z_EM(i) + ztp3 * DCOS(tp)

                        CALL GauBM5th_Hinc(i,ztp1,ztp2,ztp3)    !TEM^(y)_00~E = TEM^(x)_00~H / ztp

                        exE1x_EM(i) = exE1x_EM(i) + ztp1 / ztp * DSIN(tp) * ztponei
                        exE1y_EM(i) = exE1y_EM(i) + ztp2 / ztp * DSIN(tp) * ztponei
                        exE1z_EM(i) = exE1z_EM(i) + ztp3 / ztp * DSIN(tp) * ztponei

                        ! CALL FDMGauBM5th_dEdninc(dmx,dmy,dmz,&
                        CALL GauBM5th_dEdninc(i,ztp1,ztp2,ztp3,ztp4,ztp5,ztp6,ztp7,ztp8,ztp9)    !TEM^(x)_00~E

                        exE1xdnn_EM(i) = exE1xdnn_EM(i) + ztp1 * DCOS(tp)
                        exE1ydnn_EM(i) = exE1ydnn_EM(i) + ztp2 * DCOS(tp)
                        exE1zdnn_EM(i) = exE1zdnn_EM(i) + ztp3 * DCOS(tp)
                        exE1xdt1_EM(i) = exE1xdt1_EM(i) + ztp4 * DCOS(tp)
                        exE1ydt1_EM(i) = exE1ydt1_EM(i) + ztp5 * DCOS(tp)
                        exE1zdt1_EM(i) = exE1zdt1_EM(i) + ztp6 * DCOS(tp)
                        exE1xdt2_EM(i) = exE1xdt2_EM(i) + ztp7 * DCOS(tp)
                        exE1ydt2_EM(i) = exE1ydt2_EM(i) + ztp8 * DCOS(tp)
                        exE1zdt2_EM(i) = exE1zdt2_EM(i) + ztp9 * DCOS(tp)

                        ! CALL FDMGauBM5th_dHdninc(dmx,dmy,dmz,&
                        CALL GauBM5th_dHdninc(i,ztp1,ztp2,ztp3,ztp4,ztp5,ztp6,ztp7,ztp8,ztp9)   !TEM^(y)_00~E = TEM^(x)_00~H / ztp

                        exE1xdnn_EM(i) = exE1xdnn_EM(i) + ztp1 / ztp * DSIN(tp) * ztponei
                        exE1ydnn_EM(i) = exE1ydnn_EM(i) + ztp2 / ztp * DSIN(tp) * ztponei
                        exE1zdnn_EM(i) = exE1zdnn_EM(i) + ztp3 / ztp * DSIN(tp) * ztponei
                        exE1xdt1_EM(i) = exE1xdt1_EM(i) + ztp4 / ztp * DSIN(tp) * ztponei
                        exE1ydt1_EM(i) = exE1ydt1_EM(i) + ztp5 / ztp * DSIN(tp) * ztponei
                        exE1zdt1_EM(i) = exE1zdt1_EM(i) + ztp6 / ztp * DSIN(tp) * ztponei
                        exE1xdt2_EM(i) = exE1xdt2_EM(i) + ztp7 / ztp * DSIN(tp) * ztponei
                        exE1ydt2_EM(i) = exE1ydt2_EM(i) + ztp8 / ztp * DSIN(tp) * ztponei
                        exE1zdt2_EM(i) = exE1zdt2_EM(i) + ztp9 / ztp * DSIN(tp) * ztponei

                        CALL GauBM5th_Hinc(i,ztp1,ztp2,ztp3)    !TEM^(x)_00~H

                        exH1x_EM(i) = exH1x_EM(i) + ztp1 * DCOS(tp)
                        exH1y_EM(i) = exH1y_EM(i) + ztp2 * DCOS(tp)
                        exH1z_EM(i) = exH1z_EM(i) + ztp3 * DCOS(tp)

                        CALL GauBM5th_Einc(i,ztp1,ztp2,ztp3)    !TEM^(y)_00~H =-TEM^(x)_00~E * ztp

                        exH1x_EM(i) = exH1x_EM(i) - ztp1 * ztp * DSIN(tp) * ztponei
                        exH1y_EM(i) = exH1y_EM(i) - ztp2 * ztp * DSIN(tp) * ztponei
                        exH1z_EM(i) = exH1z_EM(i) - ztp3 * ztp * DSIN(tp) * ztponei

                        ! CALL FDMGauBM5th_dHdninc(dmx,dmy,dmz,&
                        CALL GauBM5th_dHdninc(i,ztp1,ztp2,ztp3,ztp4,ztp5,ztp6,ztp7,ztp8,ztp9)    !TEM^(x)_00~H

                        exH1xdnn_EM(i) = exH1xdnn_EM(i) + ztp1 * DCOS(tp)
                        exH1ydnn_EM(i) = exH1ydnn_EM(i) + ztp2 * DCOS(tp)
                        exH1zdnn_EM(i) = exH1zdnn_EM(i) + ztp3 * DCOS(tp)
                        exH1xdt1_EM(i) = exH1xdt1_EM(i) + ztp4 * DCOS(tp)
                        exH1ydt1_EM(i) = exH1ydt1_EM(i) + ztp5 * DCOS(tp)
                        exH1zdt1_EM(i) = exH1zdt1_EM(i) + ztp6 * DCOS(tp)
                        exH1xdt2_EM(i) = exH1xdt2_EM(i) + ztp7 * DCOS(tp)
                        exH1ydt2_EM(i) = exH1ydt2_EM(i) + ztp8 * DCOS(tp)
                        exH1zdt2_EM(i) = exH1zdt2_EM(i) + ztp9 * DCOS(tp)

                        ! CALL FDMGauBM5th_dEdninc(dmx,dmy,dmz,&
                        CALL GauBM5th_dEdninc(i,ztp1,ztp2,ztp3,ztp4,ztp5,ztp6,ztp7,ztp8,ztp9)    !TEM^(y)_00~H =-TEM^(x)_00~E * ztp

                        exH1xdnn_EM(i) = exH1xdnn_EM(i) - ztp1 * ztp * DSIN(tp) * ztponei
                        exH1ydnn_EM(i) = exH1ydnn_EM(i) - ztp2 * ztp * DSIN(tp) * ztponei
                        exH1zdnn_EM(i) = exH1zdnn_EM(i) - ztp3 * ztp * DSIN(tp) * ztponei
                        exH1xdt1_EM(i) = exH1xdt1_EM(i) - ztp4 * ztp * DSIN(tp) * ztponei
                        exH1ydt1_EM(i) = exH1ydt1_EM(i) - ztp5 * ztp * DSIN(tp) * ztponei
                        exH1zdt1_EM(i) = exH1zdt1_EM(i) - ztp6 * ztp * DSIN(tp) * ztponei
                        exH1xdt2_EM(i) = exH1xdt2_EM(i) - ztp7 * ztp * DSIN(tp) * ztponei
                        exH1ydt2_EM(i) = exH1ydt2_EM(i) - ztp8 * ztp * DSIN(tp) * ztponei
                        exH1zdt2_EM(i) = exH1zdt2_EM(i) - ztp9 * ztp * DSIN(tp) * ztponei

                    END IF

                END DO

            ELSE
            
                DO i = 1, ttlnmbrnd

                    ithprtl = 1
                    IF (nmbrprtl > 1 .AND. i > ndendID(1)) THEN
                        DO id_tp = 2, nmbrprtl
                            IF (i >= ndstaID(id_tp) .AND. i <= ndendID(id_tp)) THEN
                                ithprtl = id_tp
                                EXIT
                            END IF
                        END DO
                    END IF

                    IF (corelnkshell(ithprtl) == 0) THEN

                        dmx = xnd(i)
                        dmy = ynd(i)
                        dmz = znd(i)
                        dlh_1 = 1.0d0/dlh

                        CALL GauBM5th_Einc(i,ztp1,ztp2,ztp3)

                        exE1x_EM(i) = exE1x_EM(i) + ztp1
                        exE1y_EM(i) = exE1y_EM(i) + ztp2
                        exE1z_EM(i) = exE1z_EM(i) + ztp3

                        ! CALL FDMGauBM5th_dEdninc(dmx,dmy,dmz,&
                        CALL GauBM5th_dEdninc(i,ztp1,ztp2,ztp3,ztp4,ztp5,ztp6,ztp7,ztp8,ztp9)

                        exE1xdnn_EM(i) = exE1xdnn_EM(i) + ztp1
                        exE1ydnn_EM(i) = exE1ydnn_EM(i) + ztp2
                        exE1zdnn_EM(i) = exE1zdnn_EM(i) + ztp3
                        exE1xdt1_EM(i) = exE1xdt1_EM(i) + ztp4
                        exE1ydt1_EM(i) = exE1ydt1_EM(i) + ztp5
                        exE1zdt1_EM(i) = exE1zdt1_EM(i) + ztp6
                        exE1xdt2_EM(i) = exE1xdt2_EM(i) + ztp7
                        exE1ydt2_EM(i) = exE1ydt2_EM(i) + ztp8
                        exE1zdt2_EM(i) = exE1zdt2_EM(i) + ztp9

                        CALL GauBM5th_Hinc(i,ztp1,ztp2,ztp3)

                        exH1x_EM(i) = exH1x_EM(i) + ztp1
                        exH1y_EM(i) = exH1y_EM(i) + ztp2
                        exH1z_EM(i) = exH1z_EM(i) + ztp3

                        ! CALL FDMGauBM5th_dHdninc(dmx,dmy,dmz,&
                        CALL GauBM5th_dHdninc(i,ztp1,ztp2,ztp3,ztp4,ztp5,ztp6,ztp7,ztp8,ztp9)

                        exH1xdnn_EM(i) = exH1xdnn_EM(i) + ztp1
                        exH1ydnn_EM(i) = exH1ydnn_EM(i) + ztp2
                        exH1zdnn_EM(i) = exH1zdnn_EM(i) + ztp3
                        exH1xdt1_EM(i) = exH1xdt1_EM(i) + ztp4
                        exH1ydt1_EM(i) = exH1ydt1_EM(i) + ztp5
                        exH1zdt1_EM(i) = exH1zdt1_EM(i) + ztp6
                        exH1xdt2_EM(i) = exH1xdt2_EM(i) + ztp7
                        exH1ydt2_EM(i) = exH1ydt2_EM(i) + ztp8
                        exH1zdt2_EM(i) = exH1zdt2_EM(i) + ztp9

                    END IF

                END DO

            END IF

        END IF



    !Bessle Beam

        IF (excitetype_EM == 'bsl') THEN

            DO i = 1, ttlnmbrnd

                ithprtl = 1
                IF (nmbrprtl > 1 .AND. i > ndendID(1)) THEN
                    DO id_tp = 2, nmbrprtl
                        IF (i >= ndstaID(id_tp) .AND. i <= ndendID(id_tp)) THEN
                            ithprtl = id_tp
                            EXIT
                        END IF
                    END DO
                END IF

                IF (corelnkshell(ithprtl) == 0) THEN

                    dmx = xnd(i)
                    dmy = ynd(i)
                    dmz = znd(i)
                    dlh_1 = 1.0d0/dlh

                    CALL BsslBM_Einc(i,ztp1,ztp2,ztp3)

                    exE1x_EM(i) = exE1x_EM(i) + ztp1
                    exE1y_EM(i) = exE1y_EM(i) + ztp2
                    exE1z_EM(i) = exE1z_EM(i) + ztp3

                    CALL FDMBsslBM_dEdninc(i,ztp1,ztp2,ztp3,ztp4,ztp5,ztp6,ztp7,ztp8,ztp9)

                    exE1xdnn_EM(i) = exE1xdnn_EM(i) + ztp1
                    exE1ydnn_EM(i) = exE1ydnn_EM(i) + ztp2
                    exE1zdnn_EM(i) = exE1zdnn_EM(i) + ztp3
                    exE1xdt1_EM(i) = exE1xdt1_EM(i) + ztp4
                    exE1ydt1_EM(i) = exE1ydt1_EM(i) + ztp5
                    exE1zdt1_EM(i) = exE1zdt1_EM(i) + ztp6
                    exE1xdt2_EM(i) = exE1xdt2_EM(i) + ztp7
                    exE1ydt2_EM(i) = exE1ydt2_EM(i) + ztp8
                    exE1zdt2_EM(i) = exE1zdt2_EM(i) + ztp9

                    CALL BsslBM_Hinc(i,ztp1,ztp2,ztp3)

                    exH1x_EM(i) = exH1x_EM(i) + ztp1
                    exH1y_EM(i) = exH1y_EM(i) + ztp2
                    exH1z_EM(i) = exH1z_EM(i) + ztp3

                    CALL FDMBsslBM_dHdninc(i,ztp1,ztp2,ztp3,ztp4,ztp5,ztp6,ztp7,ztp8,ztp9)

                    exH1xdnn_EM(i) = exH1xdnn_EM(i) + ztp1
                    exH1ydnn_EM(i) = exH1ydnn_EM(i) + ztp2
                    exH1zdnn_EM(i) = exH1zdnn_EM(i) + ztp3
                    exH1xdt1_EM(i) = exH1xdt1_EM(i) + ztp4
                    exH1ydt1_EM(i) = exH1ydt1_EM(i) + ztp5
                    exH1zdt1_EM(i) = exH1zdt1_EM(i) + ztp6
                    exH1xdt2_EM(i) = exH1xdt2_EM(i) + ztp7
                    exH1ydt2_EM(i) = exH1ydt2_EM(i) + ztp8
                    exH1zdt2_EM(i) = exH1zdt2_EM(i) + ztp9

                END IF

            END DO

        END IF



    !Bessle Beam

        IF (excitetype_EM == 'mch') THEN

            DO i = 1, ttlnmbrnd

                ithprtl = 1
                IF (nmbrprtl > 1 .AND. i > ndendID(1)) THEN
                    DO id_tp = 2, nmbrprtl
                        IF (i >= ndstaID(id_tp) .AND. i <= ndendID(id_tp)) THEN
                            ithprtl = id_tp
                            EXIT
                        END IF
                    END DO
                END IF

                IF (corelnkshell(ithprtl) == 0) THEN

                    dmx = xnd(i)
                    dmy = ynd(i)
                    dmz = znd(i)
                    dlh_1 = 1.0d0/dlh

                    CALL monochromaticBM_Einc(i,ztp1,ztp2,ztp3)

                    exE1x_EM(i) = exE1x_EM(i) + ztp1
                    exE1y_EM(i) = exE1y_EM(i) + ztp2
                    exE1z_EM(i) = exE1z_EM(i) + ztp3

                    !CALL FDMmonochromaticBM_dEdninc(dmx,dmy,dmz,&
                    CALL monochromaticBM_dEdninc(i,ztp1,ztp2,ztp3,ztp4,ztp5,ztp6,&
                    &                            ztp7,ztp8,ztp9)

                    exE1xdnn_EM(i) = exE1xdnn_EM(i) + ztp1
                    exE1ydnn_EM(i) = exE1ydnn_EM(i) + ztp2
                    exE1zdnn_EM(i) = exE1zdnn_EM(i) + ztp3
                    exE1xdt1_EM(i) = exE1xdt1_EM(i) + ztp4
                    exE1ydt1_EM(i) = exE1ydt1_EM(i) + ztp5
                    exE1zdt1_EM(i) = exE1zdt1_EM(i) + ztp6
                    exE1xdt2_EM(i) = exE1xdt2_EM(i) + ztp7
                    exE1ydt2_EM(i) = exE1ydt2_EM(i) + ztp8
                    exE1zdt2_EM(i) = exE1zdt2_EM(i) + ztp9

                    CALL monochromaticBM_Hinc(i,ztp1,ztp2,ztp3)

                    exH1x_EM(i) = exH1x_EM(i) + ztp1
                    exH1y_EM(i) = exH1y_EM(i) + ztp2
                    exH1z_EM(i) = exH1z_EM(i) + ztp3

                    !CALL FDMmonochromaticBM_dHdninc(dmx,dmy,dmz,&
                    CALL monochromaticBM_dHdninc(i,ztp1,ztp2,ztp3,ztp4,ztp5,ztp6,&
                    &                            ztp7,ztp8,ztp9)

                    exH1xdnn_EM(i) = exH1xdnn_EM(i) + ztp1
                    exH1ydnn_EM(i) = exH1ydnn_EM(i) + ztp2
                    exH1zdnn_EM(i) = exH1zdnn_EM(i) + ztp3
                    exH1xdt1_EM(i) = exH1xdt1_EM(i) + ztp4
                    exH1ydt1_EM(i) = exH1ydt1_EM(i) + ztp5
                    exH1zdt1_EM(i) = exH1zdt1_EM(i) + ztp6
                    exH1xdt2_EM(i) = exH1xdt2_EM(i) + ztp7
                    exH1ydt2_EM(i) = exH1ydt2_EM(i) + ztp8
                    exH1zdt2_EM(i) = exH1zdt2_EM(i) + ztp9

                END IF

            END DO

        END IF


    END SUBROUTINE


    SUBROUTINE GauBM_Einc(dmi,dmexE1x_EM,dmexE1y_EM,dmexE1z_EM)

! Fifth-order Gaussian-beam expansion following Barton and Alexander, J. Appl. Phys. 66 (1989) 2800.
! The symbol i in the cited expressions is represented by tponej = (0,-1).
! Sign adjusted for the exp(-i*omega*t) convention used here.

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: dmi
        COMPLEX(KIND=KIND(1.0D0)), INTENT (OUT) :: dmexE1x_EM,dmexE1y_EM,dmexE1z_EM

        DOUBLE PRECISION :: dmx,dmy,dmz

        COMPLEX(KIND=KIND(1.0D0)) :: ztp, ztp1, ztp2, ztp3

        DOUBLE PRECISION :: omgo, xi, eta, zeta, ss, ll, rho, zR
        COMPLEX(KIND=KIND(1.0D0)) :: QQ, psi_0

        dmx = xnd(dmi)
        dmy = ynd(dmi)
        dmz = znd(dmi)

        omgo = incFeature_EM
        ll = ((2.0d0*pai)/(vcmwl_EM/exepsn_EM))*omgo**2
        ss = omgo/ll

        xi   = dmx/omgo
        eta  = dmy/omgo
        zeta = dmz/ll

        rho = DSQRT(xi**2 + eta**2)
        QQ = 1.0d0/(ztponej + 2.0d0*zeta)
        psi_0 = ztponej*QQ*CDEXP(-ztponej*rho**2*QQ)

        zR = (pai*(omgo**2)*exepsn_EM)/vcmwl_EM

        ztp1 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*xi**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi**2

        ztp2 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

        ztp3 = ss*(-2.0d0*QQ*xi)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*xi&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*xi

        ztp = -ztponej*(exk_EM*dmz)

        dmexE1x_EM = ztp1 * psi_0 * CDEXP(ztp)
        dmexE1y_EM = ztp2 * psi_0 * CDEXP(ztp)
        dmexE1z_EM = ztp3 * psi_0 * CDEXP(ztp)

    END SUBROUTINE

    SUBROUTINE GauBM_Hinc(dmi,dmexH1x_EM,dmexH1y_EM,dmexH1z_EM)

! Fifth-order Gaussian-beam expansion following Barton and Alexander, J. Appl. Phys. 66 (1989) 2800.
! The symbol i in the cited expressions is represented by tponej = (0,-1).
! Sign adjusted for the exp(-i*omega*t) convention used here.

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: dmi
        COMPLEX(KIND=KIND(1.0D0)), INTENT (OUT) :: dmexH1x_EM,dmexH1y_EM,dmexH1z_EM

        DOUBLE PRECISION :: dmx,dmy,dmz

        COMPLEX(KIND=KIND(1.0D0)) :: ztp, ztp1, ztp2, ztp3, ztp4
        DOUBLE PRECISION :: tp1,tp2,tp3,tp4

        DOUBLE PRECISION :: omgo, xi, eta, zeta, ss, ll, rho, zR
        COMPLEX(KIND=KIND(1.0D0)) :: QQ, psi_0

        dmx = xnd(dmi)
        dmy = ynd(dmi)
        dmz = znd(dmi)

        omgo = incFeature_EM
        ll = ((2.0d0*pai)/(vcmwl_EM/exepsn_EM))*omgo**2
        ss = omgo/ll

        xi   = dmx/omgo
        eta  = dmy/omgo
        zeta = dmz/ll

        rho = DSQRT(xi**2 + eta**2)
        QQ = 1.0d0/(ztponej + 2.0d0*zeta)
        psi_0 = ztponej*QQ*CDEXP(-ztponej*rho**2*QQ)

        zR = (pai*(omgo**2)*exepsn_EM)/vcmwl_EM

        ztp2 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*eta**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*eta**2

        ztp1 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

        ztp3 = ss*(-2.0d0*QQ*eta)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*eta&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*eta

        ztp = -ztponej*(exk_EM*dmz)

        ztp4 = (vcm_eps0*exeps_EM)/(vcm_mu0*exmiu_EM)
        tp1 =  REAL(ztp4)
        tp2 = AIMAG(ztp4)
        tp3 = (1.0d0/DSQRT(2.0d0)) * DSQRT( tp1 + DSQRT(tp1**2+tp2**2) )
        tp4 = SIGN(1.0d0,tp2)*(1.0d0/DSQRT(2.0d0)) * DSQRT( -tp1 + DSQRT(tp1**2+tp2**2) )

        ztp4 = exk_EM/(vcm_mu0*exmiu_EM*AngFrqnc_EM)
        dmexH1x_EM = ztp1 * psi_0 * CDEXP(ztp) * ztp4
        dmexH1y_EM = ztp2 * psi_0 * CDEXP(ztp) * ztp4
        dmexH1z_EM = ztp3 * psi_0 * CDEXP(ztp) * ztp4

    END SUBROUTINE





    SUBROUTINE GauBM_dEdninc(dmi,dmdexE1xdnn_EM,dmdexE1ydnn_EM,dmdexE1zdnn_EM, &
    &                           dmdexE1xdt1_EM,dmdexE1ydt1_EM,dmdexE1zdt1_EM, &
    &                           dmdexE1xdt2_EM,dmdexE1ydt2_EM,dmdexE1zdt2_EM)

! Fifth-order Gaussian-beam expansion following Barton and Alexander, J. Appl. Phys. 66 (1989) 2800.
! The symbol i in the cited expressions is represented by tponej = (0,-1).
! Sign adjusted for the exp(-i*omega*t) convention used here.

!***
!   Note:
!   drho/dxi = xi/rho;
!   drho/deta = eta/rho
!   dQQ/dzeta = -2*QQ^2
!   dpsi_0/dxi = -2*I*xi*psi_0*QQ;
!   dpsi_0/deta = -2*I*eta*psi_0*QQ
!   dpsi_0/dzeta = 2*I*(-1+rho**2+2*I*zeta)*psi_0*QQ**2
!   dxi/dx = 1/omgo
!   deta/dy = 1/omgo
!   dzeta/dz = 1/ll
!***

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: dmi
        COMPLEX(KIND=KIND(1.0D0)), INTENT (OUT) :: &
        &                           dmdexE1xdnn_EM,dmdexE1ydnn_EM,dmdexE1zdnn_EM, &
        &                           dmdexE1xdt1_EM,dmdexE1ydt1_EM,dmdexE1zdt1_EM, &
        &                           dmdexE1xdt2_EM,dmdexE1ydt2_EM,dmdexE1zdt2_EM

        DOUBLE PRECISION :: dmx,dmy,dmz,dmnnx,dmnny,dmnnz,&
        &                   dmt1x,dmt1y,dmt1z,dmt2x,dmt2y,dmt2z

        COMPLEX(KIND=KIND(1.0D0)) :: ztp, ztp1, ztp2, ztp3, ztp4
        COMPLEX(KIND=KIND(1.0D0)) :: ztpdEdx, ztpdEdy, ztpdEdz

        DOUBLE PRECISION :: omgo, xi, eta, zeta, ss, ll, rho, zR
        COMPLEX(KIND=KIND(1.0D0)) :: QQ, psi_0

        dmx = xnd(dmi)
        dmy = ynd(dmi)
        dmz = znd(dmi)
        dmnnx = nnx(dmi)
        dmnny = nny(dmi)
        dmnnz = nnz(dmi)
        dmt1x = t1x(dmi)
        dmt1y = t1y(dmi)
        dmt1z = t1z(dmi)
        dmt2x = t2x(dmi)
        dmt2y = t2y(dmi)
        dmt2z = t2z(dmi)

        omgo = incFeature_EM
        ll = ((2.0d0*pai)/(vcmwl_EM/exepsn_EM))*omgo**2
        ss = omgo/ll

        xi   = dmx/omgo
        eta  = dmy/omgo
        zeta = dmz/ll

        rho = DSQRT(xi**2 + eta**2)
        QQ = 1.0d0/(ztponej + 2.0d0*zeta)
        psi_0 = ztponej*QQ*CDEXP(-ztponej*rho**2*QQ)

        zR = (pai*(omgo**2)*exepsn_EM)/vcmwl_EM

        ztp = -ztponej*(exk_EM*dmz)

        ztpdEdx = ztpzero
        ztpdEdy = ztpzero
        ztpdEdz = ztpzero

!       d[...]/dxi in Eq.(25)
        ztp1 = (ss**2)*(-2.0d0*rho**0*QQ**2*xi+4.0d0*ztponej*rho**2*QQ**3*xi-4.0d0*QQ**2*xi**1)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4*xi-ztponej*18.0d0*rho**4*QQ**5*xi-4.0d0*rho**6*QQ**6*xi)&
        &     +(ss**4)*(16.0d0*rho**0*QQ**4*xi-ztponej*8.0d0*rho**2*QQ**5*xi)*xi**2 &
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*2.0d0*xi**1

!       d[...]/dx in Eq.(25)
        ztpdEdx = ztpdEdx + ztp1 * psi_0 * CDEXP(ztp) / omgo

        ztp1 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*xi**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi**2

!       dpsi_0/dx in Eq.(25)
        ztpdEdx = ztpdEdx + ztp1 * (-2.0d0*ztponej*xi*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/deta in Eq.(25)
        ztp1 = (ss**2)*(-2.0d0*rho**0*QQ**2*eta+4.0d0*ztponej*rho**2*QQ**3*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4*eta-ztponej*18.0d0*rho**4*QQ**5*eta-4.0d0*rho**6*QQ**6*eta)&
        &     +(ss**4)*(16.0d0*rho**0*QQ**4*eta-ztponej*8.0d0*rho**2*QQ**5*eta)*xi**2

!       d[...]/dy in Eq.(25)
        ztpdEdy = ztpdEdy + ztp1 * psi_0 * CDEXP(ztp) / omgo

        ztp1 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*xi**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi**2

!       dpsi_0/dy in Eq.(25)
        ztpdEdy = ztpdEdy + ztp1 * (-2.0d0*ztponej*eta*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/dzeta in Eq.(25)
        ztp1 = (ss**2)*(-2.0d0*rho**2*QQ**1+3.0d0*ztponej*rho**4*QQ**2-4.0d0*QQ**1*xi**2)&
        &     +(ss**4)*(8.0d0*rho**4*QQ**3-ztponej*15.0d0*rho**6*QQ**4-3.0d0*rho**8*QQ**5)&
        &     +(ss**4)*(32.0d0*rho**2*QQ**3-ztponej*10.0d0*rho**4*QQ**4)*xi**2

!       d[...]/dz in Eq.(25)
        ztpdEdz = ztpdEdz + ztp1 * (-2.0d0*QQ**2) * psi_0 * CDEXP(ztp) / ll

        ztp1 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*xi**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi**2

!       dpsi_0/dz in Eq.(25)
        ztpdEdz = ztpdEdz + ztp1 &
        &                  *2.0d0*ztponej*(-1.0d0+rho**2+2.0d0*ztponej*zeta)*psi_0*QQ**2 &
        &                  *CDEXP(ztp) / ll

        ztp1 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*xi**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi**2

!       dExp(ztp)/dz in Eq.(25)
        ztpdEdz = ztpdEdz + ztp1 * psi_0 * (-ztponej*exk_EM) * CDEXP(ztp)

!       dexE1x_EM/dn
        dmdexE1xdnn_EM = ztpdEdx * dmnnx + ztpdEdy * dmnny + ztpdEdz * dmnnz
        dmdexE1xdt1_EM = ztpdEdx * dmt1x + ztpdEdy * dmt1y + ztpdEdz * dmt1z
        dmdexE1xdt2_EM = ztpdEdx * dmt2x + ztpdEdy * dmt2y + ztpdEdz * dmt2z


        ztpdEdx = ztpzero
        ztpdEdy = ztpzero
        ztpdEdz = ztpzero

!       d[...]/dxi in Eq.(26)
        ztp2 = (ss**2)*(-2.0d0*QQ**2*eta)&
        &     +(ss**4)*(16.0d0*rho**0*QQ**4*xi-ztponej*8.0d0*rho**2*QQ**5*xi)*xi*eta&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*eta

!       d[...]/dx in Eq.(26)
        ztpdEdx = ztpdEdx + ztp2 * psi_0 * CDEXP(ztp) / omgo

        ztp2 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

!       dpsi_0/dx in Eq.(26)
        ztpdEdx = ztpdEdx + ztp2 * (-2.0d0*ztponej*xi*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/deta in Eq.(26)
        ztp2 = (ss**2)*(-2.0d0*QQ**2*xi)&
        &     +(ss**4)*(16.0d0*rho**0*QQ**4*eta-ztponej*8.0d0*rho**2*QQ**5*eta)*xi*eta&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi

!       d[...]/dy in Eq.(26)
        ztpdEdy = ztpdEdy + ztp2 * psi_0 * CDEXP(ztp) / omgo

        ztp2 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

!       dpsi_0/dy in Eq.(26)
        ztpdEdy = ztpdEdy + ztp2 * (-2.0d0*ztponej*eta*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/dzeta in Eq.(26)
        ztp2 = (ss**2)*(-4.0d0*QQ**1*xi*eta)&
        &     +(ss**4)*(32.0d0*rho**2*QQ**3-ztponej*10.0d0*rho**4*QQ**4)*xi*eta

!       d[...]/dz in Eq.(26)
        ztpdEdz = ztpdEdz + ztp2 * (-2.0d0*QQ**2) * psi_0 * CDEXP(ztp) / ll

        ztp2 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

!       dpsi_0/dz in Eq.(26)
        ztpdEdz = ztpdEdz + ztp2 &
        &                  *2.0d0*ztponej*(-1.0d0+rho**2+2.0d0*ztponej*zeta)*psi_0*QQ**2 &
        &                  *CDEXP(ztp) / ll

        ztp2 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

!       dExp(ztp)/dz in Eq.(26)
        ztpdEdz = ztpdEdz + ztp2 * psi_0 * (-ztponej*exk_EM) * CDEXP(ztp)

!       dexE1y_EM/dn
        dmdexE1ydnn_EM = ztpdEdx * dmnnx + ztpdEdy * dmnny + ztpdEdz * dmnnz
        dmdexE1ydt1_EM = ztpdEdx * dmt1x + ztpdEdy * dmt1y + ztpdEdz * dmt1z
        dmdexE1ydt2_EM = ztpdEdx * dmt2x + ztpdEdy * dmt2y + ztpdEdz * dmt2z


        ztpdEdx = ztpzero
        ztpdEdy = ztpzero
        ztpdEdz = ztpzero

!       d[...]/dxi in Eq.(27)
        ztp3 = ss*(-2.0d0*QQ)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)&
        &     +(ss**3)*(12.0d0*rho**0*QQ**3*xi-ztponej*8.0d0*rho**2*QQ**4*xi)*xi&
        &     +(ss**5)*(-80.0d0*rho**2*QQ**5*xi+ztponej*60.0d0*rho**4*QQ**6*xi+8.0d0*rho**6*QQ**7*xi)*xi

!       d[...]/dx in Eq.(27)
        ztpdEdx = ztpdEdx + ztp3 * psi_0 * CDEXP(ztp) / omgo

        ztp3 = ss*(-2.0d0*QQ*xi)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*xi&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*xi

!       dpsi_0/dx in Eq.(27)
        ztpdEdx = ztpdEdx + ztp3 * (-2.0d0*ztponej*xi*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/deta in Eq.(27)
        ztp3 = (ss**3)*(12.0d0*rho**0*QQ**3*eta-ztponej*8.0d0*rho**2*QQ**4*eta)*xi&
        &     +(ss**5)*(-80.0d0*rho**2*QQ**5*eta+ztponej*60.0d0*rho**4*QQ**6*eta+8.0d0*rho**6*QQ**7*eta)*xi

!       d[...]/dy in Eq.(27)
        ztpdEdy = ztpdEdy + ztp3 * psi_0 * CDEXP(ztp) / omgo

        ztp3 = ss*(-2.0d0*QQ*xi)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*xi&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*xi

!       dpsi_0/dy in Eq.(27)
        ztpdEdy = ztpdEdy + ztp3 * (-2.0d0*ztponej*eta*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/dzeta in Eq.(27)
        ztp3 = ss*(-2.0d0*xi)&
        &     +(ss**3)*(18.0d0*rho**2*QQ**2-ztponej*8.0d0*rho**4*QQ**3)*xi&
        &     +(ss**5)*(-100.0d0*rho**4*QQ**4+ztponej*60.0d0*rho**6*QQ**5+7.0d0*rho**8*QQ**6)*xi

!       d[...]/dz in Eq.(27)
        ztpdEdz = ztpdEdz + ztp3 * (-2.0d0*QQ**2) * psi_0 * CDEXP(ztp) / ll

        ztp3 = ss*(-2.0d0*QQ*xi)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*xi&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*xi

!       dpsi_0/dz in Eq.(26)
        ztpdEdz = ztpdEdz + ztp3 &
        &                  *2.0d0*ztponej*(-1.0d0+rho**2+2.0d0*ztponej*zeta)*psi_0*QQ**2 &
        &                  *CDEXP(ztp) / ll

        ztp3 = ss*(-2.0d0*QQ*xi)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*xi&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*xi

!       dExp(ztp)/dz in Eq.(26)
        ztpdEdz = ztpdEdz + ztp3 * psi_0 * (-ztponej*exk_EM) * CDEXP(ztp)

!       dexE1z_EM/dn
        dmdexE1zdnn_EM = ztpdEdx * dmnnx + ztpdEdy * dmnny + ztpdEdz * dmnnz
        dmdexE1zdt1_EM = ztpdEdx * dmt1x + ztpdEdy * dmt1y + ztpdEdz * dmt1z
        dmdexE1zdt2_EM = ztpdEdx * dmt2x + ztpdEdy * dmt2y + ztpdEdz * dmt2z

!----

        ztp1 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*xi**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi**2

        ztp2 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

        ztp3 = ss*(-2.0d0*QQ*xi)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*xi&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*xi

!----

    END SUBROUTINE

    SUBROUTINE GauBM_dHdninc(dmi,dmdexE1xdnn_EM,dmdexE1ydnn_EM,dmdexE1zdnn_EM, &
    &                           dmdexE1xdt1_EM,dmdexE1ydt1_EM,dmdexE1zdt1_EM, &
    &                           dmdexE1xdt2_EM,dmdexE1ydt2_EM,dmdexE1zdt2_EM)

! Fifth-order Gaussian-beam expansion following Barton and Alexander, J. Appl. Phys. 66 (1989) 2800.
! The symbol i in the cited expressions is represented by tponej = (0,-1).
! Sign adjusted for the exp(-i*omega*t) convention used here.

!***
!   Note:
!   drho/dxi = xi/rho;
!   drho/deta = eta/rho
!   dQQ/dzeta = -2*QQ^2
!   dpsi_0/dxi = -2*I*xi*psi_0*QQ;
!   dpsi_0/deta = -2*I*eta*psi_0*QQ
!   dpsi_0/dzeta = 2*I*(-1+rho**2+2*I*zeta)*psi_0*QQ**2
!   dxi/dx = 1/omgo
!   deta/dy = 1/omgo
!   dzeta/dz = 1/ll
!***

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: dmi
        COMPLEX(KIND=KIND(1.0D0)), INTENT (OUT) :: &
        &                           dmdexE1xdnn_EM,dmdexE1ydnn_EM,dmdexE1zdnn_EM, &
        &                           dmdexE1xdt1_EM,dmdexE1ydt1_EM,dmdexE1zdt1_EM, &
        &                           dmdexE1xdt2_EM,dmdexE1ydt2_EM,dmdexE1zdt2_EM

        DOUBLE PRECISION :: dmx,dmy,dmz,dmnnx,dmnny,dmnnz,&
        &                   dmt1x,dmt1y,dmt1z,dmt2x,dmt2y,dmt2z

        COMPLEX(KIND=KIND(1.0D0)) :: ztp, ztp1, ztp2, ztp3, ztp4
        COMPLEX(KIND=KIND(1.0D0)) :: ztpdEdx, ztpdEdy, ztpdEdz

        DOUBLE PRECISION :: omgo, xi, eta, zeta, ss, ll, rho, zR
        COMPLEX(KIND=KIND(1.0D0)) :: QQ, psi_0

        dmx = xnd(dmi)
        dmy = ynd(dmi)
        dmz = znd(dmi)
        dmnnx = nnx(dmi)
        dmnny = nny(dmi)
        dmnnz = nnz(dmi)
        dmt1x = t1x(dmi)
        dmt1y = t1y(dmi)
        dmt1z = t1z(dmi)
        dmt2x = t2x(dmi)
        dmt2y = t2y(dmi)
        dmt2z = t2z(dmi)

        omgo = incFeature_EM
        ll = ((2.0d0*pai)/(vcmwl_EM/exepsn_EM))*omgo**2
        ss = omgo/ll

        xi   = dmx/omgo
        eta  = dmy/omgo
        zeta = dmz/ll

        rho = DSQRT(xi**2 + eta**2)
        QQ = 1.0d0/(ztponej + 2.0d0*zeta)
        psi_0 = ztponej*QQ*CDEXP(-ztponej*rho**2*QQ)

        zR = (pai*(omgo**2)*exepsn_EM)/vcmwl_EM

        ztp = -ztponej*(exk_EM*dmz)

        ztpdEdx = ztpzero
        ztpdEdy = ztpzero
        ztpdEdz = ztpzero

!       d[...]/dxi in Eq.(28)
        ztp1 = (ss**2)*(-2.0d0*QQ**2*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*eta&
        &     +(ss**4)*(16.0d0*rho**0*QQ**4*xi-ztponej*8.0d0*rho**2*QQ**5*xi)*xi*eta

!       d[...]/dx in Eq.(28)
        ztpdEdx = ztpdEdx + ztp1 * psi_0 * CDEXP(ztp) / omgo

        ztp1 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

!       dpsi_0/dx in Eq.(28)
        ztpdEdx = ztpdEdx + ztp1 * (-2.0d0*ztponej*xi*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/deta in Eq.(28)
        ztp1 = (ss**2)*(-2.0d0*QQ**2*xi)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi&
        &     +(ss**4)*(16.0d0*rho**0*QQ**4*eta-ztponej*8.0d0*rho**2*QQ**5*eta)*xi*eta

!       d[...]/dy in Eq.(28)
        ztpdEdy = ztpdEdy + ztp1 * psi_0 * CDEXP(ztp) / omgo

        ztp1 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

!       dpsi_0/dy in Eq.(28)
        ztpdEdy = ztpdEdy + ztp1 * (-2.0d0*ztponej*eta*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/dzeta in Eq.(28)
        ztp1 = (ss**2)*(-4.0d0*QQ**1*xi*eta)&
        &     +(ss**4)*(32.0d0*rho**2*QQ**3-ztponej*10.0d0*rho**4*QQ**4)*xi*eta

!       d[...]/dz in Eq.(28)
        ztpdEdz = ztpdEdz + ztp1 * (-2.0d0*QQ**2) * psi_0 * CDEXP(ztp) / ll

        ztp1 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

!       dpsi_0/dz in Eq.(28)
        ztpdEdz = ztpdEdz + ztp1 &
        &                  *2.0d0*ztponej*(-1.0d0+rho**2+2.0d0*ztponej*zeta)*psi_0*QQ**2 &
        &                  *CDEXP(ztp) / ll

        ztp1 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

!       dExp(ztp)/dz in Eq.(28)
        ztpdEdz = ztpdEdz + ztp1 * psi_0 * (-ztponej*exk_EM) * CDEXP(ztp)

!       dexE1x_EM/dn
        dmdexE1xdnn_EM = ztpdEdx * dmnnx + ztpdEdy * dmnny + ztpdEdz * dmnnz
        dmdexE1xdt1_EM = ztpdEdx * dmt1x + ztpdEdy * dmt1y + ztpdEdz * dmt1z
        dmdexE1xdt2_EM = ztpdEdx * dmt2x + ztpdEdy * dmt2y + ztpdEdz * dmt2z


        ztpdEdx = ztpzero
        ztpdEdy = ztpzero
        ztpdEdz = ztpzero

!       d[...]/dxi in Eq.(29)
        ztp2 = (ss**2)*(-2.0d0*rho**0*QQ**2*xi+4.0d0*ztponej*rho**2*QQ**3*xi)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4*xi-ztponej*18.0d0*rho**4*QQ**5*xi-4.0d0*rho**6*QQ**6*xi)&
        &     +(ss**4)*(16.0d0*rho**0*QQ**4*xi-ztponej*8.0d0*rho**2*QQ**5*xi)*eta**2

!       d[...]/dx in Eq.(29)
        ztpdEdx = ztpdEdx + ztp2 * psi_0 * CDEXP(ztp) / omgo

        ztp2 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*eta**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*eta**2

!       dpsi_0/dx in Eq.(29)
        ztpdEdx = ztpdEdx + ztp2 * (-2.0d0*ztponej*xi*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/deta in Eq.(29)
        ztp2 = (ss**2)*(-2.0d0*rho**0*QQ**2*eta+4.0d0*ztponej*rho**2*QQ**3*eta-4.0d0*QQ**2*eta**1)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4*eta-ztponej*18.0d0*rho**4*QQ**5*eta-4.0d0*rho**6*QQ**6*eta)&
        &     +(ss**4)*(16.0d0*rho**0*QQ**4*eta-ztponej*8.0d0*rho**2*QQ**5*eta)*eta**2&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*2.0d0*eta**1

!       d[...]/dy in Eq.(29)
        ztpdEdy = ztpdEdy + ztp2 * psi_0 * CDEXP(ztp) / omgo

        ztp2 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*eta**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*eta**2

!       dpsi_0/dy in Eq.(29)
        ztpdEdy = ztpdEdy + ztp2 * (-2.0d0*ztponej*eta*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/dzeta in Eq.(26)
        ztp2 = (ss**2)*(-2.0d0*rho**2*QQ**1+3.0d0*ztponej*rho**4*QQ**2-4.0d0*QQ**1*eta**2)&
        &     +(ss**4)*(8.0d0*rho**4*QQ**3-ztponej*15.0d0*rho**6*QQ**4-3.0d0*rho**8*QQ**5)&
        &     +(ss**4)*(32.0d0*rho**2*QQ**3-ztponej*10.0d0*rho**4*QQ**4)*eta**2

!       d[...]/dz in Eq.(29)
        ztpdEdz = ztpdEdz + ztp2 * (-2.0d0*QQ**2) * psi_0 * CDEXP(ztp) / ll

        ztp2 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*eta**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*eta**2

!       dpsi_0/dz in Eq.(29)
        ztpdEdz = ztpdEdz + ztp2 &
        &                  *2.0d0*ztponej*(-1.0d0+rho**2+2.0d0*ztponej*zeta)*psi_0*QQ**2 &
        &                  *CDEXP(ztp) / ll

        ztp2 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*eta**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*eta**2

!       dExp(ztp)/dz in Eq.(29)
        ztpdEdz = ztpdEdz + ztp2 * psi_0 * (-ztponej*exk_EM) * CDEXP(ztp)

!       dexE1y_EM/dn
        dmdexE1ydnn_EM = ztpdEdx * dmnnx + ztpdEdy * dmnny + ztpdEdz * dmnnz
        dmdexE1ydt1_EM = ztpdEdx * dmt1x + ztpdEdy * dmt1y + ztpdEdz * dmt1z
        dmdexE1ydt2_EM = ztpdEdx * dmt2x + ztpdEdy * dmt2y + ztpdEdz * dmt2z


        ztpdEdx = ztpzero
        ztpdEdy = ztpzero
        ztpdEdz = ztpzero

!       d[...]/dxi in Eq.(30)
        ztp3 = (ss**3)*(12.0d0*rho**0*QQ**3*xi-ztponej*8.0d0*rho**2*QQ**4*xi)*eta&
        &     +(ss**5)*(-80.0d0*rho**2*QQ**5*xi+ztponej*60.0d0*rho**4*QQ**6*xi+8.0d0*rho**6*QQ**7*xi)*eta

!       d[...]/dx in Eq.(30)
        ztpdEdx = ztpdEdx + ztp3 * psi_0 * CDEXP(ztp) / omgo

        ztp3 = ss*(-2.0d0*QQ*eta)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*eta&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*eta

!       dpsi_0/dx in Eq.(30)
        ztpdEdx = ztpdEdx + ztp3 * (-2.0d0*ztponej*xi*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/deta in Eq.(30)
        ztp3 = ss*(-2.0d0*QQ)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)&
        &     +(ss**3)*(12.0d0*rho**0*QQ**3*eta-ztponej*8.0d0*rho**2*QQ**4*eta)*eta&
        &     +(ss**5)*(-80.0d0*rho**2*QQ**5*eta+ztponej*60.0d0*rho**4*QQ**6*eta+8.0d0*rho**6*QQ**7*eta)*eta

!       d[...]/dy in Eq.(30)
        ztpdEdy = ztpdEdy + ztp3 * psi_0 * CDEXP(ztp) / omgo

        ztp3 = ss*(-2.0d0*QQ*eta)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*eta&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*eta

!       dpsi_0/dy in Eq.(30)
        ztpdEdy = ztpdEdy + ztp3 * (-2.0d0*ztponej*eta*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/dzeta in Eq.(30)
        ztp3 = ss*(-2.0d0*eta)&
        &     +(ss**3)*(18.0d0*rho**2*QQ**2-ztponej*8.0d0*rho**4*QQ**3)*eta&
        &     +(ss**5)*(-100.0d0*rho**4*QQ**4+ztponej*60.0d0*rho**6*QQ**5+7.0d0*rho**8*QQ**6)*eta

!       d[...]/dz in Eq.(30)
        ztpdEdz = ztpdEdz + ztp3 * (-2.0d0*QQ**2) * psi_0 * CDEXP(ztp) / ll

        ztp3 = ss*(-2.0d0*QQ*eta)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*eta&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*eta

!       dpsi_0/dz in Eq.(30)
        ztpdEdz = ztpdEdz + ztp3 &
        &                  *2.0d0*ztponej*(-1.0d0+rho**2+2.0d0*ztponej*zeta)*psi_0*QQ**2 &
        &                  *CDEXP(ztp) / ll

        ztp3 = ss*(-2.0d0*QQ*eta)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*eta&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*eta

!       dExp(ztp)/dz in Eq.(30)
        ztpdEdz = ztpdEdz + ztp3 * psi_0 * (-ztponej*exk_EM) * CDEXP(ztp)

!       dexE1z_EM/dn
        dmdexE1zdnn_EM = ztpdEdx * dmnnx + ztpdEdy * dmnny + ztpdEdz * dmnnz
        dmdexE1zdt1_EM = ztpdEdx * dmt1x + ztpdEdy * dmt1y + ztpdEdz * dmt1z
        dmdexE1zdt2_EM = ztpdEdx * dmt2x + ztpdEdy * dmt2y + ztpdEdz * dmt2z

        ztp4 = exk_EM/(vcm_mu0*exmiu_EM*AngFrqnc_EM)
        dmdexE1xdnn_EM = dmdexE1xdnn_EM * ztp4
        dmdexE1ydnn_EM = dmdexE1ydnn_EM * ztp4
        dmdexE1zdnn_EM = dmdexE1zdnn_EM * ztp4
        dmdexE1xdt1_EM = dmdexE1xdt1_EM * ztp4
        dmdexE1ydt1_EM = dmdexE1ydt1_EM * ztp4
        dmdexE1zdt1_EM = dmdexE1zdt1_EM * ztp4
        dmdexE1xdt2_EM = dmdexE1xdt2_EM * ztp4
        dmdexE1ydt2_EM = dmdexE1ydt2_EM * ztp4
        dmdexE1zdt2_EM = dmdexE1zdt2_EM * ztp4

!----

        ztp2 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*eta**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*eta**2

        ztp1 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

        ztp3 = ss*(-2.0d0*QQ*eta)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*eta&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*eta

!----

    END SUBROUTINE






    SUBROUTINE GauBM5th_Einc(dmi,dmexE1x_EM,dmexE1y_EM,dmexE1z_EM)

! Fifth-order Gaussian-beam expansion following Barton and Alexander, J. Appl. Phys. 66 (1989) 2800.
! The symbol i in the cited expressions is represented by tponej = (0,-1).
! Sign adjusted for the exp(-i*omega*t) convention used here.

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: dmi
        COMPLEX(KIND=KIND(1.0D0)), INTENT (OUT) :: dmexE1x_EM,dmexE1y_EM,dmexE1z_EM

        DOUBLE PRECISION :: dmx,dmy,dmz

        COMPLEX(KIND=KIND(1.0D0)) :: ztp, ztp1, ztp2, ztp3, ztp4

        DOUBLE PRECISION :: omgo, xi, eta, zeta, ss, ll, rho, zR
        COMPLEX(KIND=KIND(1.0D0)) :: QQ, psi_0

        dmx = xnd(dmi)
        dmy = ynd(dmi)
        dmz = znd(dmi)

        omgo = incFeature_EM
        ll = ((2.0d0*pai)/(vcmwl_EM/exepsn_EM))*omgo**2
        ss = omgo/ll

        xi   = dmx/omgo
        eta  = dmy/omgo
        zeta = dmz/ll

        rho = DSQRT(xi**2 + eta**2)
        QQ = 1.0d0/(ztponej + 2.0d0*zeta)
        psi_0 = ztponej*QQ*CDEXP(-ztponej*rho**2*QQ)

        zR = (pai*(omgo**2)*exepsn_EM)/vcmwl_EM

        ztp1 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*xi**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi**2

        ztp2 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

        ztp3 = ss*(-2.0d0*QQ*xi)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*xi&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*xi

        ztp = -ztponej*(exk_EM*dmz)

        dmexE1x_EM = ztp1 * psi_0 * CDEXP(ztp)
        dmexE1y_EM = ztp2 * psi_0 * CDEXP(ztp)
        dmexE1z_EM = ztp3 * psi_0 * CDEXP(ztp)

    END SUBROUTINE

    SUBROUTINE GauBM5th_Hinc(dmi,dmexH1x_EM,dmexH1y_EM,dmexH1z_EM)

! Fifth-order Gaussian-beam expansion following Barton and Alexander, J. Appl. Phys. 66 (1989) 2800.
! The symbol i in the cited expressions is represented by tponej = (0,-1).
! Sign adjusted for the exp(-i*omega*t) convention used here.

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: dmi
        COMPLEX(KIND=KIND(1.0D0)), INTENT (OUT) :: dmexH1x_EM,dmexH1y_EM,dmexH1z_EM

        DOUBLE PRECISION :: dmx,dmy,dmz

        COMPLEX(KIND=KIND(1.0D0)) :: ztp, ztp1, ztp2, ztp3, ztp4
        DOUBLE PRECISION :: tp1,tp2,tp3,tp4

        DOUBLE PRECISION :: omgo, xi, eta, zeta, ss, ll, rho, zR
        COMPLEX(KIND=KIND(1.0D0)) :: QQ, psi_0

        dmx = xnd(dmi)
        dmy = ynd(dmi)
        dmz = znd(dmi)

        omgo = incFeature_EM
        ll = ((2.0d0*pai)/(vcmwl_EM/exepsn_EM))*omgo**2
        ss = omgo/ll

        xi   = dmx/omgo
        eta  = dmy/omgo
        zeta = dmz/ll

        rho = DSQRT(xi**2 + eta**2)
        QQ = 1.0d0/(ztponej + 2.0d0*zeta)
        psi_0 = ztponej*QQ*CDEXP(-ztponej*rho**2*QQ)

        zR = (pai*(omgo**2)*exepsn_EM)/vcmwl_EM

        ztp2 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*eta**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*eta**2

        ztp1 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

        ztp3 = ss*(-2.0d0*QQ*eta)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*eta&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*eta

        ztp = -ztponej*(exk_EM*dmz)

        ztp4 = (vcm_eps0*exeps_EM)/(vcm_mu0*exmiu_EM)
        tp1 =  REAL(ztp4)
        tp2 = AIMAG(ztp4)
        tp3 = (1.0d0/DSQRT(2.0d0)) * DSQRT( tp1 + DSQRT(tp1**2+tp2**2) )
        tp4 = SIGN(1.0d0,tp2)*(1.0d0/DSQRT(2.0d0)) * DSQRT( -tp1 + DSQRT(tp1**2+tp2**2) )

        ztp4 = exk_EM/(vcm_mu0*exmiu_EM*AngFrqnc_EM)
        dmexH1x_EM = ztp1 * psi_0 * CDEXP(ztp) * ztp4
        dmexH1y_EM = ztp2 * psi_0 * CDEXP(ztp) * ztp4
        dmexH1z_EM = ztp3 * psi_0 * CDEXP(ztp) * ztp4

    END SUBROUTINE





    SUBROUTINE GauBM5th_dEdninc(dmi,dmdexE1xdnn_EM,dmdexE1ydnn_EM,dmdexE1zdnn_EM, &
    &                           dmdexE1xdt1_EM,dmdexE1ydt1_EM,dmdexE1zdt1_EM, &
    &                           dmdexE1xdt2_EM,dmdexE1ydt2_EM,dmdexE1zdt2_EM)

! Fifth-order Gaussian-beam expansion following Barton and Alexander, J. Appl. Phys. 66 (1989) 2800.
! The symbol i in the cited expressions is represented by tponej = (0,-1).
! Sign adjusted for the exp(-i*omega*t) convention used here.

!***
!   Note:
!   drho/dxi = xi/rho;
!   drho/deta = eta/rho
!   dQQ/dzeta = -2*QQ^2
!   dpsi_0/dxi = -2*I*xi*psi_0*QQ;
!   dpsi_0/deta = -2*I*eta*psi_0*QQ
!   dpsi_0/dzeta = 2*I*(-1+rho**2+2*I*zeta)*psi_0*QQ**2
!   dxi/dx = 1/omgo
!   deta/dy = 1/omgo
!   dzeta/dz = 1/ll
!***

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: dmi
        COMPLEX(KIND=KIND(1.0D0)), INTENT (OUT) :: &
        &                           dmdexE1xdnn_EM,dmdexE1ydnn_EM,dmdexE1zdnn_EM, &
        &                           dmdexE1xdt1_EM,dmdexE1ydt1_EM,dmdexE1zdt1_EM, &
        &                           dmdexE1xdt2_EM,dmdexE1ydt2_EM,dmdexE1zdt2_EM

        DOUBLE PRECISION :: dmx,dmy,dmz,dmnnx,dmnny,dmnnz,&
        &                   dmt1x,dmt1y,dmt1z,dmt2x,dmt2y,dmt2z

        COMPLEX(KIND=KIND(1.0D0)) :: ztp, ztp1, ztp2, ztp3, ztp4
        COMPLEX(KIND=KIND(1.0D0)) :: ztpdEdx, ztpdEdy, ztpdEdz

        DOUBLE PRECISION :: omgo, xi, eta, zeta, ss, ll, rho, zR
        COMPLEX(KIND=KIND(1.0D0)) :: QQ, psi_0

        dmx = xnd(dmi)
        dmy = ynd(dmi)
        dmz = znd(dmi)
        dmnnx = nnx(dmi)
        dmnny = nny(dmi)
        dmnnz = nnz(dmi)
        dmt1x = t1x(dmi)
        dmt1y = t1y(dmi)
        dmt1z = t1z(dmi)
        dmt2x = t2x(dmi)
        dmt2y = t2y(dmi)
        dmt2z = t2z(dmi)

        omgo = incFeature_EM
        ll = ((2.0d0*pai)/(vcmwl_EM/exepsn_EM))*omgo**2
        ss = omgo/ll

        xi   = dmx/omgo
        eta  = dmy/omgo
        zeta = dmz/ll

        rho = DSQRT(xi**2 + eta**2)
        QQ = 1.0d0/(ztponej + 2.0d0*zeta)
        psi_0 = ztponej*QQ*CDEXP(-ztponej*rho**2*QQ)

        zR = (pai*(omgo**2)*exepsn_EM)/vcmwl_EM

        ztp = -ztponej*(exk_EM*dmz)

        ztpdEdx = ztpzero
        ztpdEdy = ztpzero
        ztpdEdz = ztpzero

!       d[...]/dxi in Eq.(25)
        ztp1 = (ss**2)*(-2.0d0*rho**0*QQ**2*xi+4.0d0*ztponej*rho**2*QQ**3*xi-4.0d0*QQ**2*xi**1)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4*xi-ztponej*18.0d0*rho**4*QQ**5*xi-4.0d0*rho**6*QQ**6*xi)&
        &     +(ss**4)*(16.0d0*rho**0*QQ**4*xi-ztponej*8.0d0*rho**2*QQ**5*xi)*xi**2 &
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*2.0d0*xi**1

!       d[...]/dx in Eq.(25)
        ztpdEdx = ztpdEdx + ztp1 * psi_0 * CDEXP(ztp) / omgo

        ztp1 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*xi**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi**2

!       dpsi_0/dx in Eq.(25)
        ztpdEdx = ztpdEdx + ztp1 * (-2.0d0*ztponej*xi*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/deta in Eq.(25)
        ztp1 = (ss**2)*(-2.0d0*rho**0*QQ**2*eta+4.0d0*ztponej*rho**2*QQ**3*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4*eta-ztponej*18.0d0*rho**4*QQ**5*eta-4.0d0*rho**6*QQ**6*eta)&
        &     +(ss**4)*(16.0d0*rho**0*QQ**4*eta-ztponej*8.0d0*rho**2*QQ**5*eta)*xi**2

!       d[...]/dy in Eq.(25)
        ztpdEdy = ztpdEdy + ztp1 * psi_0 * CDEXP(ztp) / omgo

        ztp1 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*xi**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi**2

!       dpsi_0/dy in Eq.(25)
        ztpdEdy = ztpdEdy + ztp1 * (-2.0d0*ztponej*eta*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/dzeta in Eq.(25)
        ztp1 = (ss**2)*(-2.0d0*rho**2*QQ**1+3.0d0*ztponej*rho**4*QQ**2-4.0d0*QQ**1*xi**2)&
        &     +(ss**4)*(8.0d0*rho**4*QQ**3-ztponej*15.0d0*rho**6*QQ**4-3.0d0*rho**8*QQ**5)&
        &     +(ss**4)*(32.0d0*rho**2*QQ**3-ztponej*10.0d0*rho**4*QQ**4)*xi**2

!       d[...]/dz in Eq.(25)
        ztpdEdz = ztpdEdz + ztp1 * (-2.0d0*QQ**2) * psi_0 * CDEXP(ztp) / ll

        ztp1 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*xi**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi**2

!       dpsi_0/dz in Eq.(25)
        ztpdEdz = ztpdEdz + ztp1 &
        &                  *2.0d0*ztponej*(-1.0d0+rho**2+2.0d0*ztponej*zeta)*psi_0*QQ**2 &
        &                  *CDEXP(ztp) / ll

        ztp1 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*xi**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi**2

!       dExp(ztp)/dz in Eq.(25)
        ztpdEdz = ztpdEdz + ztp1 * psi_0 * (-ztponej*exk_EM) * CDEXP(ztp)

!       dexE1x_EM/dn
        dmdexE1xdnn_EM = ztpdEdx * dmnnx + ztpdEdy * dmnny + ztpdEdz * dmnnz
        dmdexE1xdt1_EM = ztpdEdx * dmt1x + ztpdEdy * dmt1y + ztpdEdz * dmt1z
        dmdexE1xdt2_EM = ztpdEdx * dmt2x + ztpdEdy * dmt2y + ztpdEdz * dmt2z


        ztpdEdx = ztpzero
        ztpdEdy = ztpzero
        ztpdEdz = ztpzero

!       d[...]/dxi in Eq.(26)
        ztp2 = (ss**2)*(-2.0d0*QQ**2*eta)&
        &     +(ss**4)*(16.0d0*rho**0*QQ**4*xi-ztponej*8.0d0*rho**2*QQ**5*xi)*xi*eta&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*eta

!       d[...]/dx in Eq.(26)
        ztpdEdx = ztpdEdx + ztp2 * psi_0 * CDEXP(ztp) / omgo

        ztp2 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

!       dpsi_0/dx in Eq.(26)
        ztpdEdx = ztpdEdx + ztp2 * (-2.0d0*ztponej*xi*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/deta in Eq.(26)
        ztp2 = (ss**2)*(-2.0d0*QQ**2*xi)&
        &     +(ss**4)*(16.0d0*rho**0*QQ**4*eta-ztponej*8.0d0*rho**2*QQ**5*eta)*xi*eta&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi

!       d[...]/dy in Eq.(26)
        ztpdEdy = ztpdEdy + ztp2 * psi_0 * CDEXP(ztp) / omgo

        ztp2 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

!       dpsi_0/dy in Eq.(26)
        ztpdEdy = ztpdEdy + ztp2 * (-2.0d0*ztponej*eta*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/dzeta in Eq.(26)
        ztp2 = (ss**2)*(-4.0d0*QQ**1*xi*eta)&
        &     +(ss**4)*(32.0d0*rho**2*QQ**3-ztponej*10.0d0*rho**4*QQ**4)*xi*eta

!       d[...]/dz in Eq.(26)
        ztpdEdz = ztpdEdz + ztp2 * (-2.0d0*QQ**2) * psi_0 * CDEXP(ztp) / ll

        ztp2 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

!       dpsi_0/dz in Eq.(26)
        ztpdEdz = ztpdEdz + ztp2 &
        &                  *2.0d0*ztponej*(-1.0d0+rho**2+2.0d0*ztponej*zeta)*psi_0*QQ**2 &
        &                  *CDEXP(ztp) / ll

        ztp2 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

!       dExp(ztp)/dz in Eq.(26)
        ztpdEdz = ztpdEdz + ztp2 * psi_0 * (-ztponej*exk_EM) * CDEXP(ztp)

!       dexE1y_EM/dn
        dmdexE1ydnn_EM = ztpdEdx * dmnnx + ztpdEdy * dmnny + ztpdEdz * dmnnz
        dmdexE1ydt1_EM = ztpdEdx * dmt1x + ztpdEdy * dmt1y + ztpdEdz * dmt1z
        dmdexE1ydt2_EM = ztpdEdx * dmt2x + ztpdEdy * dmt2y + ztpdEdz * dmt2z


        ztpdEdx = ztpzero
        ztpdEdy = ztpzero
        ztpdEdz = ztpzero

!       d[...]/dxi in Eq.(27)
        ztp3 = ss*(-2.0d0*QQ)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)&
        &     +(ss**3)*(12.0d0*rho**0*QQ**3*xi-ztponej*8.0d0*rho**2*QQ**4*xi)*xi&
        &     +(ss**5)*(-80.0d0*rho**2*QQ**5*xi+ztponej*60.0d0*rho**4*QQ**6*xi+8.0d0*rho**6*QQ**7*xi)*xi

!       d[...]/dx in Eq.(27)
        ztpdEdx = ztpdEdx + ztp3 * psi_0 * CDEXP(ztp) / omgo

        ztp3 = ss*(-2.0d0*QQ*xi)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*xi&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*xi

!       dpsi_0/dx in Eq.(27)
        ztpdEdx = ztpdEdx + ztp3 * (-2.0d0*ztponej*xi*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/deta in Eq.(27)
        ztp3 = (ss**3)*(12.0d0*rho**0*QQ**3*eta-ztponej*8.0d0*rho**2*QQ**4*eta)*xi&
        &     +(ss**5)*(-80.0d0*rho**2*QQ**5*eta+ztponej*60.0d0*rho**4*QQ**6*eta+8.0d0*rho**6*QQ**7*eta)*xi

!       d[...]/dy in Eq.(27)
        ztpdEdy = ztpdEdy + ztp3 * psi_0 * CDEXP(ztp) / omgo

        ztp3 = ss*(-2.0d0*QQ*xi)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*xi&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*xi

!       dpsi_0/dy in Eq.(27)
        ztpdEdy = ztpdEdy + ztp3 * (-2.0d0*ztponej*eta*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/dzeta in Eq.(27)
        ztp3 = ss*(-2.0d0*xi)&
        &     +(ss**3)*(18.0d0*rho**2*QQ**2-ztponej*8.0d0*rho**4*QQ**3)*xi&
        &     +(ss**5)*(-100.0d0*rho**4*QQ**4+ztponej*60.0d0*rho**6*QQ**5+7.0d0*rho**8*QQ**6)*xi

!       d[...]/dz in Eq.(27)
        ztpdEdz = ztpdEdz + ztp3 * (-2.0d0*QQ**2) * psi_0 * CDEXP(ztp) / ll

        ztp3 = ss*(-2.0d0*QQ*xi)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*xi&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*xi

!       dpsi_0/dz in Eq.(26)
        ztpdEdz = ztpdEdz + ztp3 &
        &                  *2.0d0*ztponej*(-1.0d0+rho**2+2.0d0*ztponej*zeta)*psi_0*QQ**2 &
        &                  *CDEXP(ztp) / ll

        ztp3 = ss*(-2.0d0*QQ*xi)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*xi&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*xi

!       dExp(ztp)/dz in Eq.(26)
        ztpdEdz = ztpdEdz + ztp3 * psi_0 * (-ztponej*exk_EM) * CDEXP(ztp)

!       dexE1z_EM/dn
        dmdexE1zdnn_EM = ztpdEdx * dmnnx + ztpdEdy * dmnny + ztpdEdz * dmnnz
        dmdexE1zdt1_EM = ztpdEdx * dmt1x + ztpdEdy * dmt1y + ztpdEdz * dmt1z
        dmdexE1zdt2_EM = ztpdEdx * dmt2x + ztpdEdy * dmt2y + ztpdEdz * dmt2z

!----

        ztp1 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*xi**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi**2

        ztp2 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

        ztp3 = ss*(-2.0d0*QQ*xi)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*xi&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*xi

!----

    END SUBROUTINE

    SUBROUTINE GauBM5th_dHdninc(dmi,dmdexE1xdnn_EM,dmdexE1ydnn_EM,dmdexE1zdnn_EM, &
    &                           dmdexE1xdt1_EM,dmdexE1ydt1_EM,dmdexE1zdt1_EM, &
    &                           dmdexE1xdt2_EM,dmdexE1ydt2_EM,dmdexE1zdt2_EM)

! Fifth-order Gaussian-beam expansion following Barton and Alexander, J. Appl. Phys. 66 (1989) 2800.
! The symbol i in the cited expressions is represented by tponej = (0,-1).
! Sign adjusted for the exp(-i*omega*t) convention used here.

!***
!   Note:
!   drho/dxi = xi/rho;
!   drho/deta = eta/rho
!   dQQ/dzeta = -2*QQ^2
!   dpsi_0/dxi = -2*I*xi*psi_0*QQ;
!   dpsi_0/deta = -2*I*eta*psi_0*QQ
!   dpsi_0/dzeta = 2*I*(-1+rho**2+2*I*zeta)*psi_0*QQ**2
!   dxi/dx = 1/omgo
!   deta/dy = 1/omgo
!   dzeta/dz = 1/ll
!***

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: dmi
        COMPLEX(KIND=KIND(1.0D0)), INTENT (OUT) :: &
        &                           dmdexE1xdnn_EM,dmdexE1ydnn_EM,dmdexE1zdnn_EM, &
        &                           dmdexE1xdt1_EM,dmdexE1ydt1_EM,dmdexE1zdt1_EM, &
        &                           dmdexE1xdt2_EM,dmdexE1ydt2_EM,dmdexE1zdt2_EM

        DOUBLE PRECISION :: dmx,dmy,dmz,dmnnx,dmnny,dmnnz,&
        &                   dmt1x,dmt1y,dmt1z,dmt2x,dmt2y,dmt2z

        COMPLEX(KIND=KIND(1.0D0)) :: ztp, ztp1, ztp2, ztp3, ztp4
        COMPLEX(KIND=KIND(1.0D0)) :: ztpdEdx, ztpdEdy, ztpdEdz

        DOUBLE PRECISION :: omgo, xi, eta, zeta, ss, ll, rho, zR
        COMPLEX(KIND=KIND(1.0D0)) :: QQ, psi_0

        dmx = xnd(dmi)
        dmy = ynd(dmi)
        dmz = znd(dmi)
        dmnnx = nnx(dmi)
        dmnny = nny(dmi)
        dmnnz = nnz(dmi)
        dmt1x = t1x(dmi)
        dmt1y = t1y(dmi)
        dmt1z = t1z(dmi)
        dmt2x = t2x(dmi)
        dmt2y = t2y(dmi)
        dmt2z = t2z(dmi)

        omgo = incFeature_EM
        ll = ((2.0d0*pai)/(vcmwl_EM/exepsn_EM))*omgo**2
        ss = omgo/ll

        xi   = dmx/omgo
        eta  = dmy/omgo
        zeta = dmz/ll

        rho = DSQRT(xi**2 + eta**2)
        QQ = 1.0d0/(ztponej + 2.0d0*zeta)
        psi_0 = ztponej*QQ*CDEXP(-ztponej*rho**2*QQ)

        zR = (pai*(omgo**2)*exepsn_EM)/vcmwl_EM

        ztp = -ztponej*(exk_EM*dmz)

        ztpdEdx = ztpzero
        ztpdEdy = ztpzero
        ztpdEdz = ztpzero

!       d[...]/dxi in Eq.(28)
        ztp1 = (ss**2)*(-2.0d0*QQ**2*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*eta&
        &     +(ss**4)*(16.0d0*rho**0*QQ**4*xi-ztponej*8.0d0*rho**2*QQ**5*xi)*xi*eta

!       d[...]/dx in Eq.(28)
        ztpdEdx = ztpdEdx + ztp1 * psi_0 * CDEXP(ztp) / omgo

        ztp1 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

!       dpsi_0/dx in Eq.(28)
        ztpdEdx = ztpdEdx + ztp1 * (-2.0d0*ztponej*xi*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/deta in Eq.(28)
        ztp1 = (ss**2)*(-2.0d0*QQ**2*xi)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi&
        &     +(ss**4)*(16.0d0*rho**0*QQ**4*eta-ztponej*8.0d0*rho**2*QQ**5*eta)*xi*eta

!       d[...]/dy in Eq.(28)
        ztpdEdy = ztpdEdy + ztp1 * psi_0 * CDEXP(ztp) / omgo

        ztp1 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

!       dpsi_0/dy in Eq.(28)
        ztpdEdy = ztpdEdy + ztp1 * (-2.0d0*ztponej*eta*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/dzeta in Eq.(28)
        ztp1 = (ss**2)*(-4.0d0*QQ**1*xi*eta)&
        &     +(ss**4)*(32.0d0*rho**2*QQ**3-ztponej*10.0d0*rho**4*QQ**4)*xi*eta

!       d[...]/dz in Eq.(28)
        ztpdEdz = ztpdEdz + ztp1 * (-2.0d0*QQ**2) * psi_0 * CDEXP(ztp) / ll

        ztp1 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

!       dpsi_0/dz in Eq.(28)
        ztpdEdz = ztpdEdz + ztp1 &
        &                  *2.0d0*ztponej*(-1.0d0+rho**2+2.0d0*ztponej*zeta)*psi_0*QQ**2 &
        &                  *CDEXP(ztp) / ll

        ztp1 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

!       dExp(ztp)/dz in Eq.(28)
        ztpdEdz = ztpdEdz + ztp1 * psi_0 * (-ztponej*exk_EM) * CDEXP(ztp)

!       dexE1x_EM/dn
        dmdexE1xdnn_EM = ztpdEdx * dmnnx + ztpdEdy * dmnny + ztpdEdz * dmnnz
        dmdexE1xdt1_EM = ztpdEdx * dmt1x + ztpdEdy * dmt1y + ztpdEdz * dmt1z
        dmdexE1xdt2_EM = ztpdEdx * dmt2x + ztpdEdy * dmt2y + ztpdEdz * dmt2z


        ztpdEdx = ztpzero
        ztpdEdy = ztpzero
        ztpdEdz = ztpzero

!       d[...]/dxi in Eq.(29)
        ztp2 = (ss**2)*(-2.0d0*rho**0*QQ**2*xi+4.0d0*ztponej*rho**2*QQ**3*xi)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4*xi-ztponej*18.0d0*rho**4*QQ**5*xi-4.0d0*rho**6*QQ**6*xi)&
        &     +(ss**4)*(16.0d0*rho**0*QQ**4*xi-ztponej*8.0d0*rho**2*QQ**5*xi)*eta**2

!       d[...]/dx in Eq.(29)
        ztpdEdx = ztpdEdx + ztp2 * psi_0 * CDEXP(ztp) / omgo

        ztp2 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*eta**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*eta**2

!       dpsi_0/dx in Eq.(29)
        ztpdEdx = ztpdEdx + ztp2 * (-2.0d0*ztponej*xi*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/deta in Eq.(29)
        ztp2 = (ss**2)*(-2.0d0*rho**0*QQ**2*eta+4.0d0*ztponej*rho**2*QQ**3*eta-4.0d0*QQ**2*eta**1)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4*eta-ztponej*18.0d0*rho**4*QQ**5*eta-4.0d0*rho**6*QQ**6*eta)&
        &     +(ss**4)*(16.0d0*rho**0*QQ**4*eta-ztponej*8.0d0*rho**2*QQ**5*eta)*eta**2&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*2.0d0*eta**1

!       d[...]/dy in Eq.(29)
        ztpdEdy = ztpdEdy + ztp2 * psi_0 * CDEXP(ztp) / omgo

        ztp2 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*eta**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*eta**2

!       dpsi_0/dy in Eq.(29)
        ztpdEdy = ztpdEdy + ztp2 * (-2.0d0*ztponej*eta*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/dzeta in Eq.(26)
        ztp2 = (ss**2)*(-2.0d0*rho**2*QQ**1+3.0d0*ztponej*rho**4*QQ**2-4.0d0*QQ**1*eta**2)&
        &     +(ss**4)*(8.0d0*rho**4*QQ**3-ztponej*15.0d0*rho**6*QQ**4-3.0d0*rho**8*QQ**5)&
        &     +(ss**4)*(32.0d0*rho**2*QQ**3-ztponej*10.0d0*rho**4*QQ**4)*eta**2

!       d[...]/dz in Eq.(29)
        ztpdEdz = ztpdEdz + ztp2 * (-2.0d0*QQ**2) * psi_0 * CDEXP(ztp) / ll

        ztp2 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*eta**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*eta**2

!       dpsi_0/dz in Eq.(29)
        ztpdEdz = ztpdEdz + ztp2 &
        &                  *2.0d0*ztponej*(-1.0d0+rho**2+2.0d0*ztponej*zeta)*psi_0*QQ**2 &
        &                  *CDEXP(ztp) / ll

        ztp2 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*eta**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*eta**2

!       dExp(ztp)/dz in Eq.(29)
        ztpdEdz = ztpdEdz + ztp2 * psi_0 * (-ztponej*exk_EM) * CDEXP(ztp)

!       dexE1y_EM/dn
        dmdexE1ydnn_EM = ztpdEdx * dmnnx + ztpdEdy * dmnny + ztpdEdz * dmnnz
        dmdexE1ydt1_EM = ztpdEdx * dmt1x + ztpdEdy * dmt1y + ztpdEdz * dmt1z
        dmdexE1ydt2_EM = ztpdEdx * dmt2x + ztpdEdy * dmt2y + ztpdEdz * dmt2z


        ztpdEdx = ztpzero
        ztpdEdy = ztpzero
        ztpdEdz = ztpzero

!       d[...]/dxi in Eq.(30)
        ztp3 = (ss**3)*(12.0d0*rho**0*QQ**3*xi-ztponej*8.0d0*rho**2*QQ**4*xi)*eta&
        &     +(ss**5)*(-80.0d0*rho**2*QQ**5*xi+ztponej*60.0d0*rho**4*QQ**6*xi+8.0d0*rho**6*QQ**7*xi)*eta

!       d[...]/dx in Eq.(30)
        ztpdEdx = ztpdEdx + ztp3 * psi_0 * CDEXP(ztp) / omgo

        ztp3 = ss*(-2.0d0*QQ*eta)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*eta&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*eta

!       dpsi_0/dx in Eq.(30)
        ztpdEdx = ztpdEdx + ztp3 * (-2.0d0*ztponej*xi*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/deta in Eq.(30)
        ztp3 = ss*(-2.0d0*QQ)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)&
        &     +(ss**3)*(12.0d0*rho**0*QQ**3*eta-ztponej*8.0d0*rho**2*QQ**4*eta)*eta&
        &     +(ss**5)*(-80.0d0*rho**2*QQ**5*eta+ztponej*60.0d0*rho**4*QQ**6*eta+8.0d0*rho**6*QQ**7*eta)*eta

!       d[...]/dy in Eq.(30)
        ztpdEdy = ztpdEdy + ztp3 * psi_0 * CDEXP(ztp) / omgo

        ztp3 = ss*(-2.0d0*QQ*eta)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*eta&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*eta

!       dpsi_0/dy in Eq.(30)
        ztpdEdy = ztpdEdy + ztp3 * (-2.0d0*ztponej*eta*psi_0*QQ) * CDEXP(ztp) / omgo

!       d[...]/dzeta in Eq.(30)
        ztp3 = ss*(-2.0d0*eta)&
        &     +(ss**3)*(18.0d0*rho**2*QQ**2-ztponej*8.0d0*rho**4*QQ**3)*eta&
        &     +(ss**5)*(-100.0d0*rho**4*QQ**4+ztponej*60.0d0*rho**6*QQ**5+7.0d0*rho**8*QQ**6)*eta

!       d[...]/dz in Eq.(30)
        ztpdEdz = ztpdEdz + ztp3 * (-2.0d0*QQ**2) * psi_0 * CDEXP(ztp) / ll

        ztp3 = ss*(-2.0d0*QQ*eta)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*eta&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*eta

!       dpsi_0/dz in Eq.(30)
        ztpdEdz = ztpdEdz + ztp3 &
        &                  *2.0d0*ztponej*(-1.0d0+rho**2+2.0d0*ztponej*zeta)*psi_0*QQ**2 &
        &                  *CDEXP(ztp) / ll

        ztp3 = ss*(-2.0d0*QQ*eta)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*eta&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*eta

!       dExp(ztp)/dz in Eq.(30)
        ztpdEdz = ztpdEdz + ztp3 * psi_0 * (-ztponej*exk_EM) * CDEXP(ztp)

!       dexE1z_EM/dn
        dmdexE1zdnn_EM = ztpdEdx * dmnnx + ztpdEdy * dmnny + ztpdEdz * dmnnz
        dmdexE1zdt1_EM = ztpdEdx * dmt1x + ztpdEdy * dmt1y + ztpdEdz * dmt1z
        dmdexE1zdt2_EM = ztpdEdx * dmt2x + ztpdEdy * dmt2y + ztpdEdz * dmt2z

        ztp4 = exk_EM/(vcm_mu0*exmiu_EM*AngFrqnc_EM)
        dmdexE1xdnn_EM = dmdexE1xdnn_EM * ztp4
        dmdexE1ydnn_EM = dmdexE1ydnn_EM * ztp4
        dmdexE1zdnn_EM = dmdexE1zdnn_EM * ztp4
        dmdexE1xdt1_EM = dmdexE1xdt1_EM * ztp4
        dmdexE1ydt1_EM = dmdexE1ydt1_EM * ztp4
        dmdexE1zdt1_EM = dmdexE1zdt1_EM * ztp4
        dmdexE1xdt2_EM = dmdexE1xdt2_EM * ztp4
        dmdexE1ydt2_EM = dmdexE1ydt2_EM * ztp4
        dmdexE1zdt2_EM = dmdexE1zdt2_EM * ztp4

!----

        ztp2 = 1.0d0 + (ss**2)*(-rho**2*QQ**2+ztponej*rho**4*QQ**3-2.0d0*QQ**2*eta**2)&
        &     +(ss**4)*(2.0d0*rho**4*QQ**4-ztponej*3.0d0*rho**6*QQ**5-0.5d0*rho**8*QQ**6)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*eta**2

        ztp1 = (ss**2)*(-2.0d0*QQ**2*xi*eta)&
        &     +(ss**4)*(8.0d0*rho**2*QQ**4-ztponej*2.0d0*rho**4*QQ**5)*xi*eta

        ztp3 = ss*(-2.0d0*QQ*eta)&
        &     +(ss**3)*(6.0d0*rho**2*QQ**3-ztponej*2.0d0*rho**4*QQ**4)*eta&
        &     +(ss**5)*(-20.0d0*rho**4*QQ**5+ztponej*10.0d0*rho**6*QQ**6+rho**8*QQ**7)*eta

!----

    END SUBROUTINE




    SUBROUTINE BsslBM_Einc(dmi,dmexE1x_EM,dmexE1y_EM,dmexE1z_EM)

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: dmi
        COMPLEX(KIND=KIND(1.0D0)), INTENT (OUT) :: dmexE1x_EM,dmexE1y_EM,dmexE1z_EM

        DOUBLE PRECISION :: dmx,dmy,dmz

        INTEGER :: dmm

        COMPLEX(KIND=KIND(1.0D0)) :: ztp, ztp1, ztp2, ztp3, ztp4
        COMPLEX(KIND=KIND(1.0D0)) :: dmztp1, dmztp2, ztpEr, ztpEf, ztpEz

        DOUBLE PRECISION :: dmq_dmk0, dmk0, dmq, dmbeta, dmr, dmf, dmqr, dmJBS, DdmJBS1
        DOUBLE PRECISION :: tp1, tp2

        dmx = xnd(dmi)
        dmy = ynd(dmi)
        dmz = znd(dmi)

        dmztp1 = poralz_c1
        dmztp2 = poralz_c2

        dmm = incOrder_EM
        dmq_dmk0 = incFeature_EM

        dmk0 = REAL(exk_EM)
        dmq = dmq_dmk0*dmk0

        dmbeta = DSQRT(dmk0**2-dmq**2)

        dmr = DSQRT(dmx**2+dmy**2)
        IF (dmx == 0.0d0 .AND. dmy == 0.0d0) THEN
            dmf = 3.0d0*pai/4.0d0
        ELSE
            dmf = DATAN2(dmy,dmx)
        END IF

        dmqr = dmq*dmr

        dmJBS = BESSEL_JN(dmm,dmqr)
        IF (dmm == 0) THEN
            DdmJBS1 =-BESSEL_JN( 1+dmm,dmqr)
        ELSE
            DdmJBS1 = 0.5d0*(BESSEL_JN(-1+dmm,dmqr) - BESSEL_JN( 1+dmm,dmqr))
        END IF

        IF (dmm == 1) THEN
            IF (dmr == 0.0d0) THEN
                ztp = ztponei*(dmm*dmf + dmbeta*dmz)
                ztpEz = dmztp2*dmJBS * CDEXP(ztp)
                ztpEr =-(dmk0/dmq)*dmztp1*dmm*0.5d0 + (dmbeta/dmq)*dmztp2*ztponei*DdmJBS1
                ztpEr = ztpEr * CDEXP(ztp)
                ztpEf =-(dmk0/dmq)*dmztp1*ztponei*DdmJBS1 - (dmbeta/dmq)*dmztp2*dmm*0.5d0
                ztpEf = ztpEf * CDEXP(ztp)
            ELSE
                ztp = ztponei*(dmm*dmf + dmbeta*dmz)
                ztpEz = dmztp2*dmJBS * CDEXP(ztp)
                ztpEr =-(dmk0/dmq)*dmztp1*dmm/(dmqr)*dmJBS &
                &      +(dmbeta/dmq)*dmztp2*ztponei*DdmJBS1
                ztpEr = ztpEr * CDEXP(ztp)
                ztpEf =-(dmk0/dmq)*dmztp1*ztponei*DdmJBS1 &
                &      -(dmbeta/dmq)*dmztp2*dmm/(dmqr)*dmJBS
                ztpEf = ztpEf * CDEXP(ztp)
            END IF
        ELSE IF (dmm /= 1) THEN
            IF (dmr == 0.0d0) THEN
                ztp = ztponei*(dmm*dmf + dmbeta*dmz)
                ztpEz = dmztp2*dmJBS * CDEXP(ztp)
                ztpEr =-(dmk0/dmq)*dmztp1*dmm*0.0d0 + (dmbeta/dmq)*dmztp2*ztponei*DdmJBS1
                ztpEr = ztpEr * CDEXP(ztp)
                ztpEf =-(dmk0/dmq)*dmztp1*ztponei*DdmJBS1 - (dmbeta/dmq)*dmztp2*dmm*0.0d0
                ztpEf = ztpEf * CDEXP(ztp)
            ELSE
                ztp = ztponei*(dmm*dmf + dmbeta*dmz)
                ztpEz = dmztp2*dmJBS * CDEXP(ztp)
                ztpEr =-(dmk0/dmq)*dmztp1*dmm/(dmqr)*dmJBS &
                &      +(dmbeta/dmq)*dmztp2*ztponei*DdmJBS1
                ztpEr = ztpEr * CDEXP(ztp)
                ztpEf =-(dmk0/dmq)*dmztp1*ztponei*DdmJBS1 &
                &      -(dmbeta/dmq)*dmztp2*dmm/(dmqr)*dmJBS
                ztpEf = ztpEf * CDEXP(ztp)
            END IF
        END IF

        dmexE1x_EM = ztpEr*DCOS(dmf) - ztpEf*DSIN(dmf)
        dmexE1y_EM = ztpEr*DSIN(dmf) + ztpEf*DCOS(dmf)
        dmexE1z_EM = ztpEz

    END SUBROUTINE

    SUBROUTINE BsslBM_Hinc(dmi,dmexH1x_EM,dmexH1y_EM,dmexH1z_EM)

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: dmi
        COMPLEX(KIND=KIND(1.0D0)), INTENT (OUT) :: dmexH1x_EM,dmexH1y_EM,dmexH1z_EM

        DOUBLE PRECISION :: dmx,dmy,dmz

        INTEGER :: dmm

        COMPLEX(KIND=KIND(1.0D0)) :: ztp, ztp1, ztp2, ztp3, ztp4
        COMPLEX(KIND=KIND(1.0D0)) :: dmztp1, dmztp2, ztpEr, ztpEf, ztpEz

        DOUBLE PRECISION :: dmq_dmk0, dmk0, dmq, dmbeta, dmr, dmf, dmqr, dmJBS, DdmJBS1
        DOUBLE PRECISION :: tp1, tp2

        dmx = xnd(dmi)
        dmy = ynd(dmi)
        dmz = znd(dmi)

        dmztp1 = poralz_c1
        dmztp2 = poralz_c2

        dmm = incOrder_EM
        dmq_dmk0 = incFeature_EM

        dmk0 = REAL(exk_EM)
        dmq = dmq_dmk0*dmk0

        dmbeta = DSQRT(dmk0**2-dmq**2)

        dmr = DSQRT(dmx**2+dmy**2)
        IF (dmx == 0.0d0 .AND. dmy == 0.0d0) THEN
            dmf = 3.0d0*pai/4.0d0
        ELSE
            dmf = DATAN2(dmy,dmx)
        END IF

        dmqr = dmq*dmr

        dmJBS = BESSEL_JN(dmm,dmqr)
        IF (dmm == 0) THEN
            DdmJBS1 =-BESSEL_JN( 1+dmm,dmqr)
        ELSE
            DdmJBS1 = 0.5d0*(BESSEL_JN(-1+dmm,dmqr) - BESSEL_JN( 1+dmm,dmqr))
        END IF

        IF (dmm == 1) THEN
            IF (dmr == 0.0d0) THEN
                ztp = ztponei*(dmm*dmf + dmbeta*dmz)
                ztpEz = dmztp1*dmJBS * CDEXP(ztp)
                ztpEr = (dmk0/dmq)*dmztp2*dmm*0.5d0 + (dmbeta/dmq)*dmztp1*ztponei*DdmJBS1
                ztpEr = ztpEr * CDEXP(ztp)
                ztpEf = (dmk0/dmq)*dmztp2*ztponei*DdmJBS1 - (dmbeta/dmq)*dmztp1*dmm*0.5d0
                ztpEf = ztpEf * CDEXP(ztp)
            ELSE
                ztp = ztponei*(dmm*dmf + dmbeta*dmz)
                ztpEz = dmztp1*dmJBS * CDEXP(ztp)
                ztpEr = (dmk0/dmq)*dmztp2*dmm/(dmqr)*dmJBS &
                &      +(dmbeta/dmq)*dmztp1*ztponei*DdmJBS1
                ztpEr = ztpEr * CDEXP(ztp)
                ztpEf = (dmk0/dmq)*dmztp2*ztponei*DdmJBS1 &
                &      -(dmbeta/dmq)*dmztp1*dmm/(dmqr)*dmJBS
                ztpEf = ztpEf * CDEXP(ztp)
            END IF
        ELSE IF (dmm /= 1) THEN
            IF (dmr == 0.0d0) THEN
                ztp = ztponei*(dmm*dmf + dmbeta*dmz)
                ztpEz = dmztp1*dmJBS * CDEXP(ztp)
                ztpEr = (dmk0/dmq)*dmztp2*dmm*0.0d0 + (dmbeta/dmq)*dmztp1*ztponei*DdmJBS1
                ztpEr = ztpEr * CDEXP(ztp)
                ztpEf = (dmk0/dmq)*dmztp2*ztponei*DdmJBS1 - (dmbeta/dmq)*dmztp1*dmm*0.0d0
                ztpEf = ztpEf * CDEXP(ztp)
            ELSE
                ztp = ztponei*(dmm*dmf + dmbeta*dmz)
                ztpEz = dmztp1*dmJBS * CDEXP(ztp)
                ztpEr = (dmk0/dmq)*dmztp2*dmm/(dmqr)*dmJBS &
                &      +(dmbeta/dmq)*dmztp1*ztponei*DdmJBS1
                ztpEr = ztpEr * CDEXP(ztp)
                ztpEf = (dmk0/dmq)*dmztp2*ztponei*DdmJBS1 &
                &      -(dmbeta/dmq)*dmztp1*dmm/(dmqr)*dmJBS
                ztpEf = ztpEf * CDEXP(ztp)
            END IF
        END IF

        ztp = exk_EM/(vcm_mu0*exmiu_EM*AngFrqnc_EM)

        dmexH1x_EM =(ztpEr*DCOS(dmf) - ztpEf*DSIN(dmf))*ztp
        dmexH1y_EM =(ztpEr*DSIN(dmf) + ztpEf*DCOS(dmf))*ztp
        dmexH1z_EM =(ztpEz)*ztp

    END SUBROUTINE

    SUBROUTINE FDMBsslBM_dEdninc(dmi,dmdexE1xdnn_EM,dmdexE1ydnn_EM,dmdexE1zdnn_EM, &
    &                           dmdexE1xdt1_EM,dmdexE1ydt1_EM,dmdexE1zdt1_EM, &
    &                           dmdexE1xdt2_EM,dmdexE1ydt2_EM,dmdexE1zdt2_EM)

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: dmi
        COMPLEX(KIND=KIND(1.0D0)), INTENT (OUT) :: &
        &                           dmdexE1xdnn_EM,dmdexE1ydnn_EM,dmdexE1zdnn_EM, &
        &                           dmdexE1xdt1_EM,dmdexE1ydt1_EM,dmdexE1zdt1_EM, &
        &                           dmdexE1xdt2_EM,dmdexE1ydt2_EM,dmdexE1zdt2_EM

        DOUBLE PRECISION :: dmx,dmy,dmz,dmnnx,dmnny,dmnnz,&
        &                   dmt1x,dmt1y,dmt1z,dmt2x,dmt2y,dmt2z

        INTEGER :: i
        DOUBLE PRECISION :: dlh, dlh_1
        COMPLEX(KIND=KIND(1.0D0)) :: ztp,ztp1(fdm_nmx),ztp2(fdm_nmx),ztp3(fdm_nmx), &
        &             ztpdExdx,ztpdExdy,ztpdExdz,ztpdEydx,ztpdEydy,ztpdEydz, &
        &             ztpdEzdx,ztpdEzdy,ztpdEzdz

        DOUBLE PRECISION fdxgrid(fdm_nmx), fdcff(fdm_nmx), fdstt

        dmx = xnd(dmi)
        dmy = ynd(dmi)
        dmz = znd(dmi)
        dmnnx = nnx(dmi)
        dmnny = nny(dmi)
        dmnnz = nnz(dmi)
        dmt1x = t1x(dmi)
        dmt1y = t1y(dmi)
        dmt1z = t1z(dmi)
        dmt2x = t2x(dmi)
        dmt2y = t2y(dmi)
        dmt2z = t2z(dmi)

        dlh = fdm_dlh
        dlh_1 = 1.0d0/dlh

        fdstt =-dble(((fdm_nmx-1)/2))
        do i = 1, fdm_nmx
            fdxgrid(i) = fdstt + dble(i-1)
        end do
        call fdcoef(2,fdm_nmx,fdxgrid((fdm_nmx+1)/2),fdxgrid(1),fdcff)

        ztpdExdx = ztpzero
        ztpdEydx = ztpzero
        ztpdEzdx = ztpzero
        DO i = 1, fdm_nmx
            CALL BsslBMLoc_Einc(dmx+fdxgrid(i)*dlh,dmy,dmz,ztp1(i),ztp2(i),ztp3(i))
            ztpdExdx = ztpdExdx + ztp1(i)*fdcff(i)
            ztpdEydx = ztpdEydx + ztp2(i)*fdcff(i)
            ztpdEzdx = ztpdEzdx + ztp3(i)*fdcff(i)
        END DO
        ztpdExdx = ztpdExdx * dlh_1
        ztpdEydx = ztpdEydx * dlh_1
        ztpdEzdx = ztpdEzdx * dlh_1

        ztpdExdy = ztpzero
        ztpdEydy = ztpzero
        ztpdEzdy = ztpzero
        DO i = 1, fdm_nmx
            CALL BsslBMLoc_Einc(dmx,dmy+fdxgrid(i)*dlh,dmz,ztp1(i),ztp2(i),ztp3(i))
            ztpdExdy = ztpdExdy + ztp1(i)*fdcff(i)
            ztpdEydy = ztpdEydy + ztp2(i)*fdcff(i)
            ztpdEzdy = ztpdEzdy + ztp3(i)*fdcff(i)
        END DO
        ztpdExdy = ztpdExdy * dlh_1
        ztpdEydy = ztpdEydy * dlh_1
        ztpdEzdy = ztpdEzdy * dlh_1

        ztpdExdz = ztpzero
        ztpdEydz = ztpzero
        ztpdEzdz = ztpzero
        DO i = 1, fdm_nmx
            CALL BsslBMLoc_Einc(dmx,dmy,dmz+fdxgrid(i)*dlh,ztp1(i),ztp2(i),ztp3(i))
            ztpdExdz = ztpdExdz + ztp1(i)*fdcff(i)
            ztpdEydz = ztpdEydz + ztp2(i)*fdcff(i)
            ztpdEzdz = ztpdEzdz + ztp3(i)*fdcff(i)
        END DO
        ztpdExdz = ztpdExdz * dlh_1
        ztpdEydz = ztpdEydz * dlh_1
        ztpdEzdz = ztpdEzdz * dlh_1

        dmdexE1xdnn_EM = ztpdExdx * dmnnx + ztpdExdy * dmnny + ztpdExdz * dmnnz
        dmdexE1ydnn_EM = ztpdEydx * dmnnx + ztpdEydy * dmnny + ztpdEydz * dmnnz
        dmdexE1zdnn_EM = ztpdEzdx * dmnnx + ztpdEzdy * dmnny + ztpdEzdz * dmnnz
        dmdexE1xdt1_EM = ztpdExdx * dmt1x + ztpdExdy * dmt1y + ztpdExdz * dmt1z
        dmdexE1ydt1_EM = ztpdEydx * dmt1x + ztpdEydy * dmt1y + ztpdEydz * dmt1z
        dmdexE1zdt1_EM = ztpdEzdx * dmt1x + ztpdEzdy * dmt1y + ztpdEzdz * dmt1z
        dmdexE1xdt2_EM = ztpdExdx * dmt2x + ztpdExdy * dmt2y + ztpdExdz * dmt2z
        dmdexE1ydt2_EM = ztpdEydx * dmt2x + ztpdEydy * dmt2y + ztpdEydz * dmt2z
        dmdexE1zdt2_EM = ztpdEzdx * dmt2x + ztpdEzdy * dmt2y + ztpdEzdz * dmt2z

    END SUBROUTINE

    SUBROUTINE FDMBsslBM_dHdninc(dmi,dmdexE1xdnn_EM,dmdexE1ydnn_EM,dmdexE1zdnn_EM, &
    &                           dmdexE1xdt1_EM,dmdexE1ydt1_EM,dmdexE1zdt1_EM, &
    &                           dmdexE1xdt2_EM,dmdexE1ydt2_EM,dmdexE1zdt2_EM)

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: dmi
        COMPLEX(KIND=KIND(1.0D0)), INTENT (OUT) :: &
        &                           dmdexE1xdnn_EM,dmdexE1ydnn_EM,dmdexE1zdnn_EM, &
        &                           dmdexE1xdt1_EM,dmdexE1ydt1_EM,dmdexE1zdt1_EM, &
        &                           dmdexE1xdt2_EM,dmdexE1ydt2_EM,dmdexE1zdt2_EM

        DOUBLE PRECISION :: dmx,dmy,dmz,dmnnx,dmnny,dmnnz,&
        &                   dmt1x,dmt1y,dmt1z,dmt2x,dmt2y,dmt2z

        INTEGER :: i
        DOUBLE PRECISION :: dlh, dlh_1
        COMPLEX(KIND=KIND(1.0D0)) :: ztp,ztp1(fdm_nmx),ztp2(fdm_nmx),ztp3(fdm_nmx), &
        &             ztpdExdx,ztpdExdy,ztpdExdz,ztpdEydx,ztpdEydy,ztpdEydz, &
        &             ztpdEzdx,ztpdEzdy,ztpdEzdz

        DOUBLE PRECISION fdxgrid(fdm_nmx), fdcff(fdm_nmx), fdstt

        dmx = xnd(dmi)
        dmy = ynd(dmi)
        dmz = znd(dmi)
        dmnnx = nnx(dmi)
        dmnny = nny(dmi)
        dmnnz = nnz(dmi)
        dmt1x = t1x(dmi)
        dmt1y = t1y(dmi)
        dmt1z = t1z(dmi)
        dmt2x = t2x(dmi)
        dmt2y = t2y(dmi)
        dmt2z = t2z(dmi)

        dlh = fdm_dlh
        dlh_1 = 1.0d0/dlh

        fdstt =-dble(((fdm_nmx-1)/2))
        do i = 1, fdm_nmx
            fdxgrid(i) = fdstt + dble(i-1)
        end do
        call fdcoef(2,fdm_nmx,fdxgrid((fdm_nmx+1)/2),fdxgrid(1),fdcff)

        ztpdExdx = ztpzero
        ztpdEydx = ztpzero
        ztpdEzdx = ztpzero
        DO i = 1, fdm_nmx
            CALL BsslBMLoc_Hinc(dmx+fdxgrid(i)*dlh,dmy,dmz,ztp1(i),ztp2(i),ztp3(i))
            ztpdExdx = ztpdExdx + ztp1(i)*fdcff(i)
            ztpdEydx = ztpdEydx + ztp2(i)*fdcff(i)
            ztpdEzdx = ztpdEzdx + ztp3(i)*fdcff(i)
        END DO
        ztpdExdx = ztpdExdx * dlh_1
        ztpdEydx = ztpdEydx * dlh_1
        ztpdEzdx = ztpdEzdx * dlh_1

        ztpdExdy = ztpzero
        ztpdEydy = ztpzero
        ztpdEzdy = ztpzero
        DO i = 1, fdm_nmx
            CALL BsslBMLoc_Hinc(dmx,dmy+fdxgrid(i)*dlh,dmz,ztp1(i),ztp2(i),ztp3(i))
            ztpdExdy = ztpdExdy + ztp1(i)*fdcff(i)
            ztpdEydy = ztpdEydy + ztp2(i)*fdcff(i)
            ztpdEzdy = ztpdEzdy + ztp3(i)*fdcff(i)
        END DO
        ztpdExdy = ztpdExdy * dlh_1
        ztpdEydy = ztpdEydy * dlh_1
        ztpdEzdy = ztpdEzdy * dlh_1

        ztpdExdz = ztpzero
        ztpdEydz = ztpzero
        ztpdEzdz = ztpzero
        DO i = 1, fdm_nmx
            CALL BsslBMLoc_Hinc(dmx,dmy,dmz+fdxgrid(i)*dlh,ztp1(i),ztp2(i),ztp3(i))
            ztpdExdz = ztpdExdz + ztp1(i)*fdcff(i)
            ztpdEydz = ztpdEydz + ztp2(i)*fdcff(i)
            ztpdEzdz = ztpdEzdz + ztp3(i)*fdcff(i)
        END DO
        ztpdExdz = ztpdExdz * dlh_1
        ztpdEydz = ztpdEydz * dlh_1
        ztpdEzdz = ztpdEzdz * dlh_1

        dmdexE1xdnn_EM = ztpdExdx * dmnnx + ztpdExdy * dmnny + ztpdExdz * dmnnz
        dmdexE1ydnn_EM = ztpdEydx * dmnnx + ztpdEydy * dmnny + ztpdEydz * dmnnz
        dmdexE1zdnn_EM = ztpdEzdx * dmnnx + ztpdEzdy * dmnny + ztpdEzdz * dmnnz
        dmdexE1xdt1_EM = ztpdExdx * dmt1x + ztpdExdy * dmt1y + ztpdExdz * dmt1z
        dmdexE1ydt1_EM = ztpdEydx * dmt1x + ztpdEydy * dmt1y + ztpdEydz * dmt1z
        dmdexE1zdt1_EM = ztpdEzdx * dmt1x + ztpdEzdy * dmt1y + ztpdEzdz * dmt1z
        dmdexE1xdt2_EM = ztpdExdx * dmt2x + ztpdExdy * dmt2y + ztpdExdz * dmt2z
        dmdexE1ydt2_EM = ztpdEydx * dmt2x + ztpdEydy * dmt2y + ztpdEydz * dmt2z
        dmdexE1zdt2_EM = ztpdEzdx * dmt2x + ztpdEzdy * dmt2y + ztpdEzdz * dmt2z

    END SUBROUTINE



    SUBROUTINE monochromaticBM_Einc(dmi,dmexE1x_EM,dmexE1y_EM,dmexE1z_EM)

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: dmi
        COMPLEX(KIND=KIND(1.0D0)), INTENT (OUT) :: dmexE1x_EM,dmexE1y_EM,dmexE1z_EM

        DOUBLE PRECISION :: dmx,dmy,dmz

        INTEGER :: dmm

        COMPLEX(KIND=KIND(1.0D0)) :: ztp, ztp1, ztp2, ztp3, ztp4
        COMPLEX(KIND=KIND(1.0D0)) :: dmztp1, dmztp2, ztpEr, ztpEf, ztpEz

        DOUBLE PRECISION :: dmq_dmk0, dmk0, dmq, dmbeta, dmr, dmf, dmqr, dmJBS, DdmJBS1
        DOUBLE PRECISION :: tp1, tp2

        dmx = xnd(dmi)
        dmy = ynd(dmi)
        dmz = znd(dmi)

        dmztp1 = poralz_c1
        dmztp2 = poralz_c2

        dmq_dmk0 = incFeature_EM

        dmk0 = REAL(exk_EM)
        dmq = dmq_dmk0*dmk0

        dmbeta = DSQRT(1.0d0-dmq_dmk0**2)

        ztp = ztponei*(dmbeta*dmk0*dmx)

        dmexE1x_EM =-ztponei*dmq_dmk0*dmztp2*DSIN(dmq_dmk0*dmk0*dmy) * CDEXP(ztp)
        dmexE1y_EM = dmbeta*dmztp2*DCOS(dmq_dmk0*dmk0*dmy) * CDEXP(ztp)
        dmexE1z_EM = dmztp1*DCOS(dmq_dmk0*dmk0*dmy) * CDEXP(ztp)

    END SUBROUTINE

    SUBROUTINE monochromaticBM_Hinc(dmi,dmexH1x_EM,dmexH1y_EM,dmexH1z_EM)

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: dmi
        COMPLEX(KIND=KIND(1.0D0)), INTENT (OUT) :: dmexH1x_EM,dmexH1y_EM,dmexH1z_EM

        DOUBLE PRECISION :: dmx,dmy,dmz

        INTEGER :: dmm

        COMPLEX(KIND=KIND(1.0D0)) :: ztp, ztp1, ztp2, ztp3, ztp4
        COMPLEX(KIND=KIND(1.0D0)) :: dmztp1, dmztp2, ztpEr, ztpEf, ztpEz

        DOUBLE PRECISION :: dmq_dmk0, dmk0, dmq, dmbeta, dmr, dmf, dmqr, dmJBS, DdmJBS1
        DOUBLE PRECISION :: tp1, tp2

        dmx = xnd(dmi)
        dmy = ynd(dmi)
        dmz = znd(dmi)

        dmztp1 = poralz_c1
        dmztp2 = poralz_c2

        dmq_dmk0 = incFeature_EM

        dmk0 = REAL(exk_EM)
        dmq = dmq_dmk0*dmk0

        dmbeta = DSQRT(1.0d0-dmq_dmk0**2)

        ztp = ztponei*(dmbeta*dmk0*dmx)

        dmexH1x_EM = ztponei*dmq_dmk0*dmztp1*DSIN(dmq_dmk0*dmk0*dmy) * CDEXP(ztp)
        dmexH1y_EM =-dmbeta*dmztp1*DCOS(dmq_dmk0*dmk0*dmy) * CDEXP(ztp)
        dmexH1z_EM = dmztp2*DCOS(dmq_dmk0*dmk0*dmy) * CDEXP(ztp)

        ztp = exk_EM/(vcm_mu0*exmiu_EM*AngFrqnc_EM)

        dmexH1x_EM = dmexH1x_EM * ztp
        dmexH1y_EM = dmexH1y_EM * ztp
        dmexH1z_EM = dmexH1z_EM * ztp

    END SUBROUTINE

    SUBROUTINE monochromaticBM_dEdninc(dmi,dmdexE1xdnn_EM,dmdexE1ydnn_EM,dmdexE1zdnn_EM, &
    &                           dmdexE1xdt1_EM,dmdexE1ydt1_EM,dmdexE1zdt1_EM, &
    &                           dmdexE1xdt2_EM,dmdexE1ydt2_EM,dmdexE1zdt2_EM)

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: dmi
        COMPLEX(KIND=KIND(1.0D0)), INTENT (OUT) :: &
        &                           dmdexE1xdnn_EM,dmdexE1ydnn_EM,dmdexE1zdnn_EM, &
        &                           dmdexE1xdt1_EM,dmdexE1ydt1_EM,dmdexE1zdt1_EM, &
        &                           dmdexE1xdt2_EM,dmdexE1ydt2_EM,dmdexE1zdt2_EM

        DOUBLE PRECISION :: dmx,dmy,dmz,dmnnx,dmnny,dmnnz,&
        &                   dmt1x,dmt1y,dmt1z,dmt2x,dmt2y,dmt2z

        INTEGER :: dmm

        COMPLEX(KIND=KIND(1.0D0)) :: ztp, ztp1, ztp2, ztp3, ztp4
        COMPLEX(KIND=KIND(1.0D0)) :: dmztp1, dmztp2, ztpEr, ztpEf, ztpEz

        DOUBLE PRECISION :: dmq_dmk0, dmk0, dmq, dmbeta, dmr, dmf, dmqr, dmJBS, DdmJBS1
        DOUBLE PRECISION :: tp1, tp2

        dmx = xnd(dmi)
        dmy = ynd(dmi)
        dmz = znd(dmi)
        dmnnx = nnx(dmi)
        dmnny = nny(dmi)
        dmnnz = nnz(dmi)
        dmt1x = t1x(dmi)
        dmt1y = t1y(dmi)
        dmt1z = t1z(dmi)
        dmt2x = t2x(dmi)
        dmt2y = t2y(dmi)
        dmt2z = t2z(dmi)

        dmztp1 = poralz_c1
        dmztp2 = poralz_c2

        dmq_dmk0 = incFeature_EM

        dmk0 = REAL(exk_EM)
        dmq = dmq_dmk0*dmk0

        dmbeta = DSQRT(1.0d0-dmq_dmk0**2)

        ztp = ztponei*(dmbeta*dmk0*dmx)

        dmdexE1xdnn_EM =-ztponei*dmq_dmk0*dmztp2*DSIN(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *ztponei*dmbeta*dmk0*dmnnx &
        &           -ztponei*dmq_dmk0*dmztp2*DCOS(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *dmq_dmk0*dmk0*dmnny
        dmdexE1ydnn_EM = dmbeta*dmztp2*DCOS(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *ztponei*dmbeta*dmk0*dmnnx &
        &           -dmbeta*dmztp2*DSIN(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *dmq_dmk0*dmk0*dmnny
        dmdexE1zdnn_EM = dmztp1*DCOS(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *ztponei*dmbeta*dmk0*dmnnx &
        &           -dmztp1*DSIN(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *dmq_dmk0*dmk0*dmnny
        dmdexE1xdt1_EM=-ztponei*dmq_dmk0*dmztp2*DSIN(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *ztponei*dmbeta*dmk0*dmt1x &
        &           -ztponei*dmq_dmk0*dmztp2*DCOS(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *dmq_dmk0*dmk0*dmt1y
        dmdexE1ydt1_EM= dmbeta*dmztp2*DCOS(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *ztponei*dmbeta*dmk0*dmt1x &
        &           -dmbeta*dmztp2*DSIN(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *dmq_dmk0*dmk0*dmt1y
        dmdexE1zdt1_EM= dmztp1*DCOS(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *ztponei*dmbeta*dmk0*dmt1x &
        &           -dmztp1*DSIN(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *dmq_dmk0*dmk0*dmt1y
        dmdexE1xdt2_EM=-ztponei*dmq_dmk0*dmztp2*DSIN(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *ztponei*dmbeta*dmk0*dmt2x &
        &           -ztponei*dmq_dmk0*dmztp2*DCOS(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *dmq_dmk0*dmk0*dmt2y
        dmdexE1ydt2_EM= dmbeta*dmztp2*DCOS(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *ztponei*dmbeta*dmk0*dmt2x &
        &           -dmbeta*dmztp2*DSIN(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *dmq_dmk0*dmk0*dmt2y
        dmdexE1zdt2_EM= dmztp1*DCOS(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *ztponei*dmbeta*dmk0*dmt2x &
        &           -dmztp1*DSIN(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *dmq_dmk0*dmk0*dmt2y

    END SUBROUTINE

    SUBROUTINE monochromaticBM_dHdninc(dmi,dmdexE1xdnn_EM,dmdexE1ydnn_EM,dmdexE1zdnn_EM, &
    &                           dmdexE1xdt1_EM,dmdexE1ydt1_EM,dmdexE1zdt1_EM, &
    &                           dmdexE1xdt2_EM,dmdexE1ydt2_EM,dmdexE1zdt2_EM)

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: dmi
        COMPLEX(KIND=KIND(1.0D0)), INTENT (OUT) :: &
        &                           dmdexE1xdnn_EM,dmdexE1ydnn_EM,dmdexE1zdnn_EM, &
        &                           dmdexE1xdt1_EM,dmdexE1ydt1_EM,dmdexE1zdt1_EM, &
        &                           dmdexE1xdt2_EM,dmdexE1ydt2_EM,dmdexE1zdt2_EM

        DOUBLE PRECISION :: dmx,dmy,dmz,dmnnx,dmnny,dmnnz,&
        &                   dmt1x,dmt1y,dmt1z,dmt2x,dmt2y,dmt2z

        INTEGER :: dmm

        COMPLEX(KIND=KIND(1.0D0)) :: ztp, ztp1, ztp2, ztp3, ztp4
        COMPLEX(KIND=KIND(1.0D0)) :: dmztp1, dmztp2, ztpEr, ztpEf, ztpEz

        DOUBLE PRECISION :: dmq_dmk0, dmk0, dmq, dmbeta, dmr, dmf, dmqr, dmJBS, DdmJBS1
        DOUBLE PRECISION :: tp1, tp2

        dmx = xnd(dmi)
        dmy = ynd(dmi)
        dmz = znd(dmi)
        dmnnx = nnx(dmi)
        dmnny = nny(dmi)
        dmnnz = nnz(dmi)
        dmt1x = t1x(dmi)
        dmt1y = t1y(dmi)
        dmt1z = t1z(dmi)
        dmt2x = t2x(dmi)
        dmt2y = t2y(dmi)
        dmt2z = t2z(dmi)

        dmztp1 = poralz_c1
        dmztp2 = poralz_c2

        dmq_dmk0 = incFeature_EM

        dmk0 = REAL(exk_EM)
        dmq = dmq_dmk0*dmk0

        dmbeta = DSQRT(1.0d0-dmq_dmk0**2)

        ztp = ztponei*(dmbeta*dmk0*dmx)

        dmdexE1xdnn_EM = ztponei*dmq_dmk0*dmztp1*DSIN(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *ztponei*dmbeta*dmk0*dmnnx &
        &           +ztponei*dmq_dmk0*dmztp1*DCOS(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *dmq_dmk0*dmk0*dmnny
        dmdexE1ydnn_EM =-dmbeta*dmztp1*DCOS(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *ztponei*dmbeta*dmk0*dmnnx &
        &           +dmbeta*dmztp1*DSIN(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *dmq_dmk0*dmk0*dmnny
        dmdexE1zdnn_EM = dmztp2*DCOS(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *ztponei*dmbeta*dmk0*dmnnx &
        &           -dmztp2*DSIN(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *dmq_dmk0*dmk0*dmnny
        dmdexE1xdt1_EM= ztponei*dmq_dmk0*dmztp1*DSIN(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *ztponei*dmbeta*dmk0*dmt1x &
        &           +ztponei*dmq_dmk0*dmztp1*DCOS(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *dmq_dmk0*dmk0*dmt1y
        dmdexE1ydt1_EM=-dmbeta*dmztp1*DCOS(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *ztponei*dmbeta*dmk0*dmt1x &
        &           +dmbeta*dmztp1*DSIN(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *dmq_dmk0*dmk0*dmt1y
        dmdexE1zdt1_EM= dmztp2*DCOS(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *ztponei*dmbeta*dmk0*dmt1x &
        &           -dmztp2*DSIN(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *dmq_dmk0*dmk0*dmt1y
        dmdexE1xdt2_EM= ztponei*dmq_dmk0*dmztp1*DSIN(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *ztponei*dmbeta*dmk0*dmt2x &
        &           +ztponei*dmq_dmk0*dmztp1*DCOS(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *dmq_dmk0*dmk0*dmt2y
        dmdexE1ydt2_EM=-dmbeta*dmztp1*DCOS(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *ztponei*dmbeta*dmk0*dmt2x &
        &           +dmbeta*dmztp1*DSIN(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *dmq_dmk0*dmk0*dmt2y
        dmdexE1zdt2_EM= dmztp2*DCOS(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *ztponei*dmbeta*dmk0*dmt2x &
        &           -dmztp2*DSIN(dmq_dmk0*dmk0*dmy) * CDEXP(ztp) &
        &            *dmq_dmk0*dmk0*dmt2y

        ztp = exk_EM/(vcm_mu0*exmiu_EM*AngFrqnc_EM)

        dmdexE1xdnn_EM = dmdexE1xdnn_EM * ztp
        dmdexE1ydnn_EM = dmdexE1ydnn_EM * ztp
        dmdexE1zdnn_EM = dmdexE1zdnn_EM * ztp
        dmdexE1xdt1_EM = dmdexE1xdt1_EM * ztp
        dmdexE1ydt1_EM = dmdexE1ydt1_EM * ztp
        dmdexE1zdt1_EM = dmdexE1zdt1_EM * ztp
        dmdexE1xdt2_EM = dmdexE1xdt2_EM * ztp
        dmdexE1ydt2_EM = dmdexE1ydt2_EM * ztp
        dmdexE1zdt2_EM = dmdexE1zdt2_EM * ztp

    END SUBROUTINE












    SUBROUTINE BsslBMLoc_Einc(dmx,dmy,dmz,dmexE1x_EM,dmexE1y_EM,dmexE1z_EM)

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT (IN) :: dmx,dmy,dmz
        COMPLEX(KIND=KIND(1.0D0)), INTENT (OUT) :: dmexE1x_EM,dmexE1y_EM,dmexE1z_EM

        INTEGER :: dmm

        COMPLEX(KIND=KIND(1.0D0)) :: ztp, ztp1, ztp2, ztp3, ztp4
        COMPLEX(KIND=KIND(1.0D0)) :: dmztp1, dmztp2, ztpEr, ztpEf, ztpEz

        DOUBLE PRECISION :: dmq_dmk0, dmk0, dmq, dmbeta, dmr, dmf, dmqr, dmJBS, DdmJBS1
        DOUBLE PRECISION :: tp1, tp2

        dmztp1 = poralz_c1
        dmztp2 = poralz_c2

        dmm = incOrder_EM
        dmq_dmk0 = incFeature_EM

        dmk0 = REAL(exk_EM)
        dmq = dmq_dmk0*dmk0

        dmbeta = DSQRT(dmk0**2-dmq**2)

        dmr = DSQRT(dmx**2+dmy**2)
        IF (dmx == 0.0d0 .AND. dmy == 0.0d0) THEN
            dmf = 3.0d0*pai/4.0d0
        ELSE
            dmf = DATAN2(dmy,dmx)
        END IF

        dmqr = dmq*dmr

        dmJBS = BESSEL_JN(dmm,dmqr)
        IF (dmm == 0) THEN
            DdmJBS1 =-BESSEL_JN( 1+dmm,dmqr)
        ELSE
            DdmJBS1 = 0.5d0*(BESSEL_JN(-1+dmm,dmqr) - BESSEL_JN( 1+dmm,dmqr))
        END IF

        IF (dmm == 1) THEN
            IF (dmr == 0.0d0) THEN
                ztp = ztponei*(dmm*dmf + dmbeta*dmz)
                ztpEz = dmztp2*dmJBS * CDEXP(ztp)
                ztpEr =-(dmk0/dmq)*dmztp1*dmm*0.5d0 + (dmbeta/dmq)*dmztp2*ztponei*DdmJBS1
                ztpEr = ztpEr * CDEXP(ztp)
                ztpEf =-(dmk0/dmq)*dmztp1*ztponei*DdmJBS1 - (dmbeta/dmq)*dmztp2*dmm*0.5d0
                ztpEf = ztpEf * CDEXP(ztp)
            ELSE
                ztp = ztponei*(dmm*dmf + dmbeta*dmz)
                ztpEz = dmztp2*dmJBS * CDEXP(ztp)
                ztpEr =-(dmk0/dmq)*dmztp1*dmm/(dmqr)*dmJBS &
                &      +(dmbeta/dmq)*dmztp2*ztponei*DdmJBS1
                ztpEr = ztpEr * CDEXP(ztp)
                ztpEf =-(dmk0/dmq)*dmztp1*ztponei*DdmJBS1 &
                &      -(dmbeta/dmq)*dmztp2*dmm/(dmqr)*dmJBS
                ztpEf = ztpEf * CDEXP(ztp)
            END IF
        ELSE IF (dmm /= 1) THEN
            IF (dmr == 0.0d0) THEN
                ztp = ztponei*(dmm*dmf + dmbeta*dmz)
                ztpEz = dmztp2*dmJBS * CDEXP(ztp)
                ztpEr =-(dmk0/dmq)*dmztp1*dmm*0.0d0 + (dmbeta/dmq)*dmztp2*ztponei*DdmJBS1
                ztpEr = ztpEr * CDEXP(ztp)
                ztpEf =-(dmk0/dmq)*dmztp1*ztponei*DdmJBS1 - (dmbeta/dmq)*dmztp2*dmm*0.0d0
                ztpEf = ztpEf * CDEXP(ztp)
            ELSE
                ztp = ztponei*(dmm*dmf + dmbeta*dmz)
                ztpEz = dmztp2*dmJBS * CDEXP(ztp)
                ztpEr =-(dmk0/dmq)*dmztp1*dmm/(dmqr)*dmJBS &
                &      +(dmbeta/dmq)*dmztp2*ztponei*DdmJBS1
                ztpEr = ztpEr * CDEXP(ztp)
                ztpEf =-(dmk0/dmq)*dmztp1*ztponei*DdmJBS1 &
                &      -(dmbeta/dmq)*dmztp2*dmm/(dmqr)*dmJBS
                ztpEf = ztpEf * CDEXP(ztp)
            END IF
        END IF

        dmexE1x_EM = ztpEr*DCOS(dmf) - ztpEf*DSIN(dmf)
        dmexE1y_EM = ztpEr*DSIN(dmf) + ztpEf*DCOS(dmf)
        dmexE1z_EM = ztpEz

    END SUBROUTINE

    SUBROUTINE BsslBMLoc_Hinc(dmx,dmy,dmz,dmexH1x_EM,dmexH1y_EM,dmexH1z_EM)

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT (IN) :: dmx,dmy,dmz
        COMPLEX(KIND=KIND(1.0D0)), INTENT (OUT) :: dmexH1x_EM,dmexH1y_EM,dmexH1z_EM

        INTEGER :: dmm

        COMPLEX(KIND=KIND(1.0D0)) :: ztp, ztp1, ztp2, ztp3, ztp4
        COMPLEX(KIND=KIND(1.0D0)) :: dmztp1, dmztp2, ztpEr, ztpEf, ztpEz

        DOUBLE PRECISION :: dmq_dmk0, dmk0, dmq, dmbeta, dmr, dmf, dmqr, dmJBS, DdmJBS1
        DOUBLE PRECISION :: tp1, tp2

        dmztp1 = poralz_c1
        dmztp2 = poralz_c2

        dmm = incOrder_EM
        dmq_dmk0 = incFeature_EM

        dmk0 = REAL(exk_EM)
        dmq = dmq_dmk0*dmk0

        dmbeta = DSQRT(dmk0**2-dmq**2)

        dmr = DSQRT(dmx**2+dmy**2)
        IF (dmx == 0.0d0 .AND. dmy == 0.0d0) THEN
            dmf = 3.0d0*pai/4.0d0
        ELSE
            dmf = DATAN2(dmy,dmx)
        END IF

        dmqr = dmq*dmr

        dmJBS = BESSEL_JN(dmm,dmqr)
        IF (dmm == 0) THEN
            DdmJBS1 =-BESSEL_JN( 1+dmm,dmqr)
        ELSE
            DdmJBS1 = 0.5d0*(BESSEL_JN(-1+dmm,dmqr) - BESSEL_JN( 1+dmm,dmqr))
        END IF

        IF (dmm == 1) THEN
            IF (dmr == 0.0d0) THEN
                ztp = ztponei*(dmm*dmf + dmbeta*dmz)
                ztpEz = dmztp1*dmJBS * CDEXP(ztp)
                ztpEr = (dmk0/dmq)*dmztp2*dmm*0.5d0 + (dmbeta/dmq)*dmztp1*ztponei*DdmJBS1
                ztpEr = ztpEr * CDEXP(ztp)
                ztpEf = (dmk0/dmq)*dmztp2*ztponei*DdmJBS1 - (dmbeta/dmq)*dmztp1*dmm*0.5d0
                ztpEf = ztpEf * CDEXP(ztp)
            ELSE
                ztp = ztponei*(dmm*dmf + dmbeta*dmz)
                ztpEz = dmztp1*dmJBS * CDEXP(ztp)
                ztpEr = (dmk0/dmq)*dmztp2*dmm/(dmqr)*dmJBS &
                &      +(dmbeta/dmq)*dmztp1*ztponei*DdmJBS1
                ztpEr = ztpEr * CDEXP(ztp)
                ztpEf = (dmk0/dmq)*dmztp2*ztponei*DdmJBS1 &
                &      -(dmbeta/dmq)*dmztp1*dmm/(dmqr)*dmJBS
                ztpEf = ztpEf * CDEXP(ztp)
            END IF
        ELSE IF (dmm /= 1) THEN
            IF (dmr == 0.0d0) THEN
                ztp = ztponei*(dmm*dmf + dmbeta*dmz)
                ztpEz = dmztp1*dmJBS * CDEXP(ztp)
                ztpEr = (dmk0/dmq)*dmztp2*dmm*0.0d0 + (dmbeta/dmq)*dmztp1*ztponei*DdmJBS1
                ztpEr = ztpEr * CDEXP(ztp)
                ztpEf = (dmk0/dmq)*dmztp2*ztponei*DdmJBS1 - (dmbeta/dmq)*dmztp1*dmm*0.0d0
                ztpEf = ztpEf * CDEXP(ztp)
            ELSE
                ztp = ztponei*(dmm*dmf + dmbeta*dmz)
                ztpEz = dmztp1*dmJBS * CDEXP(ztp)
                ztpEr = (dmk0/dmq)*dmztp2*dmm/(dmqr)*dmJBS &
                &      +(dmbeta/dmq)*dmztp1*ztponei*DdmJBS1
                ztpEr = ztpEr * CDEXP(ztp)
                ztpEf = (dmk0/dmq)*dmztp2*ztponei*DdmJBS1 &
                &      -(dmbeta/dmq)*dmztp1*dmm/(dmqr)*dmJBS
                ztpEf = ztpEf * CDEXP(ztp)
            END IF
        END IF

        ztp = exk_EM/(vcm_mu0*exmiu_EM*AngFrqnc_EM)

        dmexH1x_EM =(ztpEr*DCOS(dmf) - ztpEf*DSIN(dmf))*ztp
        dmexH1y_EM =(ztpEr*DSIN(dmf) + ztpEf*DCOS(dmf))*ztp
        dmexH1z_EM =(ztpEz)*ztp

    END SUBROUTINE



END MODULE
