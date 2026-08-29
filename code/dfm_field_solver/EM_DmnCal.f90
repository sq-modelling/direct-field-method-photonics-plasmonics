
! SPDX-FileCopyrightText: 2026 Qiang Sun
! SPDX-License-Identifier: BSD-3-Clause

! Evaluate complex electric or magnetic fields at off-surface observation
! points from solved surface values and normal derivatives. The same regularised
! Helmholtz representation is used for exterior and bounded material domains.
!
MODULE EM_DmnCal

    USE omp_lib

    USE Pre_Constants
    USE Pre_csvformat

    USE Geom_GlobalData

    USE EM_SurfCal_GlobalData


    IMPLICIT NONE
    CONTAINS

    ! Evaluate one Cartesian field vector at a specified observation point.
    ! dmEorH selects E or H; dmdomain selects the exterior or a material domain.
    ! The three complex output arguments contain the x, y, and z components.

    SUBROUTINE GetDmnE2H2_EM (dmEorH, dmDmnx, dmDmny, dmDmnz, dmdomain, &
    &                         dmEHxdmn, dmEHydmn, dmEHzdmn)

        IMPLICIT NONE

        CHARACTER (LEN=1), INTENT(IN) :: dmEorH
        DOUBLE PRECISION, INTENT(IN) :: dmDmnx, dmDmny, dmDmnz
        INTEGER, INTENT(IN) :: dmdomain
        COMPLEX(KIND=KIND(1.0D0)), INTENT(OUT) :: dmEHxdmn, dmEHydmn, dmEHzdmn

        CHARACTER (LEN=2) :: bdry_type

        INTEGER :: dmprtlid, dmx0id, dmbdryid

        INTEGER :: ithprtl, jthprtl, kthprtl, offstNd, p0id

        INTEGER :: i, j, k, ii, jj, kk

        INTEGER :: x0id, srpj, x0elmnt = 26, elmntA, xpx0_elmnt

        INTEGER :: NdA, NdB, NdC, NdD, NdE, NdF

        DOUBLE PRECISION :: pr0x, pr0y, pr0z, pr0nnx, pr0nny, pr0nnz, tpx, tpy, tpz

        DOUBLE PRECISION :: tp,tp1,tp2,tp3,tp4,tp5,tp6,tp_a,tp_b,tp_g,tpxi,tpet,xircd,etrcd
        DOUBLE PRECISION :: tprx,tpry,tprz,tpEnx,tpEny,tpEnz,JcbDtmn
        DOUBLE PRECISION :: tpendx1,tpendx2,tpendx3,tpendx4,tpendx5,tpendx6
        DOUBLE PRECISION :: tpendy1,tpendy2,tpendy3,tpendy4,tpendy5,tpendy6
        DOUBLE PRECISION :: tpendz1,tpendz2,tpendz3,tpendz4,tpendz5,tpendz6
        DOUBLE PRECISION :: drx_deps,drx_dyet,dry_deps,dry_dyet,drz_deps,drz_dyet

        COMPLEX(KIND=KIND(1.0D0)) :: ztp,ztp1,ztp2,ztp3,ztp4,ztp5,ztp6,ztp7,ztp8,ztp9

        COMPLEX(KIND=KIND(1.0D0)) :: ztpEHx0, ztpEHy0, ztpEHz0, ztpEHxdn0, ztpEHydn0, ztpEHzdn0

        COMPLEX(KIND=KIND(1.0D0)) :: ztpEHx_elmnt, ztpEHy_elmnt, ztpEHz_elmnt

        DOUBLE PRECISION :: mdl_rdmnrd

        DOUBLE PRECISION :: dstntp, mindstntp

        DOUBLE PRECISION :: tpphi_elmnt


        mindstntp = 1E30

        ! Identify the nearest boundary node used by the regularising auxiliary
        ! field and determine whether exterior or interior traces apply.

        DO ithprtl = 1, nmbrprtl

            IF (dmdomain == ithprtl .OR. dmdomain == corelnkshell(ithprtl)) THEN

                DO i = ndstaID(ithprtl), ndendID(ithprtl)
                    dstntp =  (xnd(i)-dmDmnx)**2 &
                            &+(ynd(i)-dmDmny)**2 &
                            &+(znd(i)-dmDmnz)**2
                    dstntp = DSQRT(dstntp)
                    IF (dstntp < mindstntp) THEN
                        mindstntp = dstntp
                        dmx0id = i
                        dmbdryid = ithprtl
                    END IF
                END DO
            END IF
        END DO

        pr0x = xnd(dmx0id)
        pr0y = ynd(dmx0id)
        pr0z = znd(dmx0id)
        pr0nnx = nnx(dmx0id)
        pr0nny = nny(dmx0id)
        pr0nnz = nnz(dmx0id)
        IF (dmdomain == corelnkshell(dmbdryid)) THEN
            IF (dmEorH == 'E') THEN
                ztpEHx0 = exE2x_EM(dmx0id)
                ztpEHy0 = exE2y_EM(dmx0id)
                ztpEHz0 = exE2z_EM(dmx0id)
                ztpEHxdn0 = exE2xdnn_EM(dmx0id)
                ztpEHydn0 = exE2ydnn_EM(dmx0id)
                ztpEHzdn0 = exE2zdnn_EM(dmx0id)
            END IF
            IF (dmEorH == 'H') THEN
                ztpEHx0 = exH2x_EM(dmx0id)
                ztpEHy0 = exH2y_EM(dmx0id)
                ztpEHz0 = exH2z_EM(dmx0id)
                ztpEHxdn0 = exH2xdnn_EM(dmx0id)
                ztpEHydn0 = exH2ydnn_EM(dmx0id)
                ztpEHzdn0 = exH2zdnn_EM(dmx0id)
            END IF
            bdry_type = 'ex'
        END IF
        IF (dmdomain == dmbdryid) THEN
            IF (dmEorH == 'E') THEN
                ztpEHx0 = inE2x_EM(dmx0id)
                ztpEHy0 = inE2y_EM(dmx0id)
                ztpEHz0 = inE2z_EM(dmx0id)
                ztpEHxdn0 = inE2xdnn_EM(dmx0id)
                ztpEHydn0 = inE2ydnn_EM(dmx0id)
                ztpEHzdn0 = inE2zdnn_EM(dmx0id)
            END IF
            IF (dmEorH == 'H') THEN
                ztpEHx0 = inH2x_EM(dmx0id)
                ztpEHy0 = inH2y_EM(dmx0id)
                ztpEHz0 = inH2z_EM(dmx0id)
                ztpEHxdn0 = inH2xdnn_EM(dmx0id)
                ztpEHydn0 = inH2ydnn_EM(dmx0id)
                ztpEHzdn0 = inH2zdnn_EM(dmx0id)
            END IF
            bdry_type = 'in'
        END IF

! Determine the maximum edge length among elements incident on x0.

        tp = 0.0d0
        DO j = 1, mxnmbrndlnkelmnt
            IF (ndlnkelmnt(dmx0id, j) /= 0) THEN
                k = ndlnkelmnt(dmx0id, j)
                IF (MeshType == 'L') THEN
                    DO i = 1, 3
                        NdA = elmntlnknd(k,i)
                        tp1 = (xnd(NdA)-pr0x)**2 &
                        &    +(ynd(NdA)-pr0y)**2 &
                        &    +(znd(NdA)-pr0z)**2
                        tp1 = DSQRT(tp1)
                        tp = MAX(tp,tp1)
                    END DO
                END IF
                IF (MeshType == 'Q') THEN
                    DO i = 1, 6
                        NdA = elmntlnknd(k,i)
                        tp1 = (xnd(NdA)-pr0x)**2 &
                        &    +(ynd(NdA)-pr0y)**2 &
                        &    +(znd(NdA)-pr0z)**2
                        tp1 = DSQRT(tp1)
                        tp = MAX(tp,tp1)
                    END DO
                END IF
            END IF
        END DO


! Detect near-boundary targets that are not aligned with the nearest node normal.
! Such targets require a surface-based regularization point.

        xpx0_elmnt = 0
        IF (mindstntp < 3.0d0*tp) THEN
            tp1 = (pr0x-dmDmnx)*pr0nnx + (pr0y-dmDmny)*pr0nny + (pr0z-dmDmnz)*pr0nnz
            tp1 = DABS(tp1)/mindstntp
            IF (DABS(tp1 - 1.0d0) > 1.d-3) xpx0_elmnt = 1
        END IF



! For a nearby off-normal target, locate the closest surface element.

        IF (xpx0_elmnt == 1) THEN

!-
! Search the candidate object for the closest element.

            i = dmx0id

            mindstntp = 1E30

            DO j = 1, mxnmbrndlnkelmnt

                IF (ndlnkelmnt(i, j) /= 0) THEN

                    k = ndlnkelmnt(i, j)

                    IF (MeshType == 'L') THEN
                        NdA = elmntlnknd(k,1)
                        NdB = elmntlnknd(k,2)
                        NdC = elmntlnknd(k,3)
                        tpx = 1.0d0/3.0d0 * (xnd(NdA)+xnd(NdB)+xnd(NdC))
                        tpy = 1.0d0/3.0d0 * (ynd(NdA)+ynd(NdB)+ynd(NdC))
                        tpz = 1.0d0/3.0d0 * (znd(NdA)+znd(NdB)+znd(NdC))
                        dstntp =  (tpx-dmDmnx)**2 &
                                &+(tpy-dmDmny)**2 &
                                &+(tpz-dmDmnz)**2
                        dstntp = DSQRT(dstntp)
                        IF (dstntp < mindstntp) THEN
                            mindstntp = dstntp
                            dmx0id = k
                            dmbdryid = 1
                        END IF
                    END IF

                    IF (MeshType == 'Q') THEN
                        NdA = elmntlnknd(k,4)
                        NdB = elmntlnknd(k,5)
                        NdC = elmntlnknd(k,6)
                        tpx = 1.0d0/3.0d0 * (xnd(NdA)+xnd(NdB)+xnd(NdC))
                        tpy = 1.0d0/3.0d0 * (ynd(NdA)+ynd(NdB)+ynd(NdC))
                        tpz = 1.0d0/3.0d0 * (znd(NdA)+znd(NdB)+znd(NdC))
                        dstntp =  (tpx-dmDmnx)**2 &
                                &+(tpy-dmDmny)**2 &
                                &+(tpz-dmDmnz)**2
                        dstntp = DSQRT(dstntp)
                        IF (dstntp < mindstntp) THEN
                            mindstntp = dstntp
                            dmx0id = k
                            dmbdryid = 1
                        END IF
                    END IF

                END IF

            END DO

!-

!-
! Refine the regularization point by subdividing the closest element.

            k = dmx0id

            IF (MeshType == 'L') THEN

                NdA = elmntlnknd(k,1)
                NdB = elmntlnknd(k,2)
                NdC = elmntlnknd(k,3)

                tpendx1 = xnd(elmntlnknd(k,1))
                tpendx2 = xnd(elmntlnknd(k,2))
                tpendx3 = xnd(elmntlnknd(k,3))

                tpendy1 = ynd(elmntlnknd(k,1))
                tpendy2 = ynd(elmntlnknd(k,2))
                tpendy3 = ynd(elmntlnknd(k,3))

                tpendz1 = znd(elmntlnknd(k,1))
                tpendz2 = znd(elmntlnknd(k,2))
                tpendz3 = znd(elmntlnknd(k,3))

                mindstntp = 1E30
                do i=1,x0elmnt
                    do j=1,x0elmnt

                        tpxi = dble(i-1)/dble(x0elmnt-1)
                        tpet = dble(j-1)/dble(x0elmnt-1)

                        tp2 = tpxi
                        tp3 = tpet
                        tp1 = 1.0d0 - tp2 - tp3

                        tpx = tp1*tpendx1+tp2*tpendx2+tp3*tpendx3
                        tpy = tp1*tpendy1+tp2*tpendy2+tp3*tpendy3
                        tpz = tp1*tpendz1+tp2*tpendz2+tp3*tpendz3

                        dstntp = (tpx-dmDmnx)**2+(tpy-dmDmny)**2+(tpz-dmDmnz)**2
                        dstntp = DSQRT(dstntp)
                        if (dstntp < mindstntp) then
                            mindstntp = dstntp
                            xircd = tpxi
                            etrcd = tpet
                        end if

                    end do
                end do

                tpxi = xircd
                tpet = etrcd

                tp2 = tpxi
                tp3 = tpet
                tp1 = 1.0d0-tp2-tp3

                pr0x = tp1*xnd(NdA)+tp2*xnd(NdB)+tp3*xnd(NdC)
                pr0y = tp1*ynd(NdA)+tp2*ynd(NdB)+tp3*ynd(NdC)
                pr0z = tp1*znd(NdA)+tp2*znd(NdB)+tp3*znd(NdC)

                pr0nnx = nnxelmnt(k)
                pr0nny = nnyelmnt(k)
                pr0nnz = nnzelmnt(k)

                tp = dsqrt(pr0nnx**2 + pr0nny**2 + pr0nnz**2)
                pr0nnx = pr0nnx/tp
                pr0nny = pr0nny/tp
                pr0nnz = pr0nnz/tp

                IF (bdry_type == 'ex') THEN
                    IF (dmEorH == 'E') THEN
                        ztpEHx0 = tp1*exE2x_EM(NdA)+tp2*exE2x_EM(NdB) &
                        &        +tp3*exE2x_EM(NdC)
                        ztpEHxdn0= tp1*exE2xdnn_EM(NdA)+tp2*exE2xdnn_EM(NdB) &
                        &         +tp3*exE2xdnn_EM(NdC)
                        ztpEHy0 = tp1*exE2y_EM(NdA)+tp2*exE2y_EM(NdB) &
                        &        +tp3*exE2y_EM(NdC)
                        ztpEHydn0= tp1*exE2ydnn_EM(NdA)+tp2*exE2ydnn_EM(NdB) &
                        &         +tp3*exE2ydnn_EM(NdC)
                        ztpEHz0 = tp1*exE2z_EM(NdA)+tp2*exE2z_EM(NdB) &
                        &        +tp3*exE2z_EM(NdC)
                        ztpEHzdn0= tp1*exE2zdnn_EM(NdA)+tp2*exE2zdnn_EM(NdB) &
                        &         +tp3*exE2zdnn_EM(NdC)
                    END IF
                    IF (dmEorH == 'H') THEN
                        ztpEHx0 = tp1*exH2x_EM(NdA)+tp2*exH2x_EM(NdB) &
                        &        +tp3*exH2x_EM(NdC)
                        ztpEHxdn0= tp1*exH2xdnn_EM(NdA)+tp2*exH2xdnn_EM(NdB) &
                        &         +tp3*exH2xdnn_EM(NdC)
                        ztpEHy0 = tp1*exH2y_EM(NdA)+tp2*exH2y_EM(NdB) &
                        &        +tp3*exH2y_EM(NdC)
                        ztpEHydn0= tp1*exH2ydnn_EM(NdA)+tp2*exH2ydnn_EM(NdB) &
                        &         +tp3*exH2ydnn_EM(NdC)
                        ztpEHz0 = tp1*exH2z_EM(NdA)+tp2*exH2z_EM(NdB) &
                        &        +tp3*exH2z_EM(NdC)
                        ztpEHzdn0= tp1*exH2zdnn_EM(NdA)+tp2*exH2zdnn_EM(NdB) &
                        &         +tp3*exH2zdnn_EM(NdC)
                    END IF
                END IF
                IF (bdry_type == 'in') THEN
                    IF (dmEorH == 'E') THEN
                        ztpEHx0 = tp1*inE2x_EM(NdA)+tp2*inE2x_EM(NdB) &
                        &        +tp3*inE2x_EM(NdC)
                        ztpEHxdn0= tp1*inE2xdnn_EM(NdA)+tp2*inE2xdnn_EM(NdB) &
                        &         +tp3*inE2xdnn_EM(NdC)
                        ztpEHy0 = tp1*inE2y_EM(NdA)+tp2*inE2y_EM(NdB) &
                        &        +tp3*inE2y_EM(NdC)
                        ztpEHydn0= tp1*inE2ydnn_EM(NdA)+tp2*inE2ydnn_EM(NdB) &
                        &         +tp3*inE2ydnn_EM(NdC)
                        ztpEHz0 = tp1*inE2z_EM(NdA)+tp2*inE2z_EM(NdB) &
                        &        +tp3*inE2z_EM(NdC)
                        ztpEHzdn0= tp1*inE2zdnn_EM(NdA)+tp2*inE2zdnn_EM(NdB) &
                        &         +tp3*inE2zdnn_EM(NdC)
                    END IF
                    IF (dmEorH == 'H') THEN
                        ztpEHx0 = tp1*inH2x_EM(NdA)+tp2*inH2x_EM(NdB) &
                        &        +tp3*inH2x_EM(NdC)
                        ztpEHxdn0= tp1*inH2xdnn_EM(NdA)+tp2*inH2xdnn_EM(NdB) &
                        &         +tp3*inH2xdnn_EM(NdC)
                        ztpEHy0 = tp1*inH2y_EM(NdA)+tp2*inH2y_EM(NdB) &
                        &        +tp3*inH2y_EM(NdC)
                        ztpEHydn0= tp1*inH2ydnn_EM(NdA)+tp2*inH2ydnn_EM(NdB) &
                        &         +tp3*inH2ydnn_EM(NdC)
                        ztpEHz0 = tp1*inH2z_EM(NdA)+tp2*inH2z_EM(NdB) &
                        &        +tp3*inH2z_EM(NdC)
                        ztpEHzdn0= tp1*inH2zdnn_EM(NdA)+tp2*inH2zdnn_EM(NdB) &
                        &         +tp3*inH2zdnn_EM(NdC)
                    END IF
                END IF

            END IF
            IF (MeshType == 'Q') THEN

                NdA = elmntlnknd(k,1)
                NdB = elmntlnknd(k,2)
                NdC = elmntlnknd(k,3)
                NdD = elmntlnknd(k,4)
                NdE = elmntlnknd(k,5)
                NdF = elmntlnknd(k,6)

                tpendx1 = xnd(elmntlnknd(k,1))
                tpendx2 = xnd(elmntlnknd(k,2))
                tpendx3 = xnd(elmntlnknd(k,3))
                tpendx4 = xnd(elmntlnknd(k,4))
                tpendx5 = xnd(elmntlnknd(k,5))
                tpendx6 = xnd(elmntlnknd(k,6))
                tpendy1 = ynd(elmntlnknd(k,1))
                tpendy2 = ynd(elmntlnknd(k,2))
                tpendy3 = ynd(elmntlnknd(k,3))
                tpendy4 = ynd(elmntlnknd(k,4))
                tpendy5 = ynd(elmntlnknd(k,5))
                tpendy6 = ynd(elmntlnknd(k,6))
                tpendz1 = znd(elmntlnknd(k,1))
                tpendz2 = znd(elmntlnknd(k,2))
                tpendz3 = znd(elmntlnknd(k,3))
                tpendz4 = znd(elmntlnknd(k,4))
                tpendz5 = znd(elmntlnknd(k,5))
                tpendz6 = znd(elmntlnknd(k,6))

                tp1 = DSQRT( (tpendx4-tpendx2)**2&
                &           +(tpendy4-tpendy2)**2&
                &           +(tpendz4-tpendz2)**2 )
                tp2 = DSQRT( (tpendx4-tpendx1)**2&
                &           +(tpendy4-tpendy1)**2&
                &           +(tpendz4-tpendz1)**2 )
                tp_a = 1.0d0/(1.0d0 + tp1/tp2)
                tp1 = DSQRT( (tpendx6-tpendx3)**2&
                &           +(tpendy6-tpendy3)**2&
                &           +(tpendz6-tpendz3)**2 )
                tp2 = DSQRT( (tpendx6-tpendx1)**2&
                &           +(tpendy6-tpendy1)**2&
                &           +(tpendz6-tpendz1)**2 )
                tp_b = 1.0d0/(1.0d0 + tp1/tp2)
                tp1 = DSQRT( (tpendx5-tpendx2)**2&
                &           +(tpendy5-tpendy2)**2&
                &           +(tpendz5-tpendz2)**2 )
                tp2 = DSQRT( (tpendx5-tpendx3)**2&
                &           +(tpendy5-tpendy3)**2&
                &           +(tpendz5-tpendz3)**2 )
                tp_g = 1.0d0/(1.0d0 + tp1/tp2)

                mindstntp = 1E30
                do i=1,x0elmnt
                    do j=1,x0elmnt

                        tpxi = dble(i-1)/dble(x0elmnt-1)
                        tpet = dble(j-1)/dble(x0elmnt-1)

                        tp2 = 1.0d0/(1.0d0-tp_a)*tpxi*(tpxi-tp_a+(tp_a-tp_g)/(1.0d0-tp_g)*tpet)
                        tp3 = 1.0d0/(1.0d0-tp_b)*tpet*(tpet-tp_b+(tp_b+tp_g-1.0d0)/(tp_g)*tpxi)
                        tp4 = 1.0d0/(tp_a*(1.0d0-tp_a))*tpxi*(1.0d0-tpxi-tpet)
                        tp5 = 1.0d0/(tp_g*(1.0d0-tp_g))*tpxi*tpet
                        tp6 = 1.0d0/(tp_b*(1.0d0-tp_b))*tpet*(1.0d0-tpxi-tpet)
                        tp1 = 1.0d0-tp2-tp3-tp4-tp5-tp6

                        tpx = tp1*tpendx1+tp2*tpendx2+tp3*tpendx3 &
                        &    +tp4*tpendx4+tp5*tpendx5+tp6*tpendx6
                        tpy = tp1*tpendy1+tp2*tpendy2+tp3*tpendy3 &
                        &    +tp4*tpendy4+tp5*tpendy5+tp6*tpendy6
                        tpz = tp1*tpendz1+tp2*tpendz2+tp3*tpendz3 &
                        &    +tp4*tpendz4+tp5*tpendz5+tp6*tpendz6

                        dstntp = (tpx-dmDmnx)**2+(tpy-dmDmny)**2+(tpz-dmDmnz)**2
                        dstntp = DSQRT(dstntp)
                        if (dstntp < mindstntp) then
                            mindstntp = dstntp
                            xircd = tpxi
                            etrcd = tpet
                        end if

                    end do
                end do

                tpxi = xircd
                tpet = etrcd

                tp2= 1.0d0/(1.0d0-tp_a)*(tpxi-tp_a+(tp_a-tp_g)/(1.0d0-tp_g)*tpet+tpxi)
                tp3= 1.0d0/(1.0d0-tp_b)*tpet*(tp_b+tp_g-1.0d0)/(tp_g)
                tp4= 1.0d0/(tp_a*(1.0d0-tp_a))*(1.0d0-tpxi-tpet-tpxi)
                tp5= 1.0d0/(tp_g*(1.0d0-tp_g))*tpet
                tp6=-1.0d0/(tp_b*(1.0d0-tp_b))*tpet
                tp1=-tp2-tp3-tp4-tp5-tp6

                drx_deps = tp1*tpendx1+tp2*tpendx2+tp3*tpendx3 &
                &         +tp4*tpendx4+tp5*tpendx5+tp6*tpendx6
                dry_deps = tp1*tpendy1+tp2*tpendy2+tp3*tpendy3 &
                &         +tp4*tpendy4+tp5*tpendy5+tp6*tpendy6
                drz_deps = tp1*tpendz1+tp2*tpendz2+tp3*tpendz3 &
                &         +tp4*tpendz4+tp5*tpendz5+tp6*tpendz6

                tp2= 1.0d0/(1.0d0-tp_a)*tpxi*(tp_a-tp_g)/(1.0d0-tp_g)
                tp3= 1.0d0/(1.0d0-tp_b)*(tpet-tp_b+(tp_b+tp_g-1.0d0)/(tp_g)*tpxi+tpet)
                tp4=-1.0d0/(tp_a*(1.0d0-tp_a))*tpxi
                tp5= 1.0d0/(tp_g*(1.0d0-tp_g))*tpxi
                tp6= 1.0d0/(tp_b*(1.0d0-tp_b))*(1.0d0-tpxi-tpet-tpet)
                tp1=-tp2-tp3-tp4-tp5-tp6

                drx_dyet = tp1*tpendx1+tp2*tpendx2+tp3*tpendx3 &
                &         +tp4*tpendx4+tp5*tpendx5+tp6*tpendx6
                dry_dyet = tp1*tpendy1+tp2*tpendy2+tp3*tpendy3 &
                &         +tp4*tpendy4+tp5*tpendy5+tp6*tpendy6
                drz_dyet = tp1*tpendz1+tp2*tpendz2+tp3*tpendz3 &
                &         +tp4*tpendz4+tp5*tpendz5+tp6*tpendz6

                tpEnx = dry_deps*drz_dyet - drz_deps*dry_dyet
                tpEny = drz_deps*drx_dyet - drx_deps*drz_dyet
                tpEnz = drx_deps*dry_dyet - dry_deps*drx_dyet

                JcbDtmn = DSQRT(tpEnx**2+tpEny**2+tpEnz**2)

                tpEnx = tpEnx/JcbDtmn
                tpEny = tpEny/JcbDtmn
                tpEnz = tpEnz/JcbDtmn

                tp2 = 1.0d0/(1.0d0-tp_a)*tpxi*(tpxi-tp_a+(tp_a-tp_g)/(1.0d0-tp_g)*tpet)
                tp3 = 1.0d0/(1.0d0-tp_b)*tpet*(tpet-tp_b+(tp_b+tp_g-1.0d0)/(tp_g)*tpxi)
                tp4 = 1.0d0/(tp_a*(1.0d0-tp_a))*tpxi*(1.0d0-tpxi-tpet)
                tp5 = 1.0d0/(tp_g*(1.0d0-tp_g))*tpxi*tpet
                tp6 = 1.0d0/(tp_b*(1.0d0-tp_b))*tpet*(1.0d0-tpxi-tpet)
                tp1 = 1.0d0-tp2-tp3-tp4-tp5-tp6

                pr0x = tp1*xnd(NdA)+tp2*xnd(NdB)+tp3*xnd(NdC) &
                &     +tp4*xnd(NdD)+tp5*xnd(NdE)+tp6*xnd(NdF)
                pr0y = tp1*ynd(NdA)+tp2*ynd(NdB)+tp3*ynd(NdC) &
                &     +tp4*ynd(NdD)+tp5*ynd(NdE)+tp6*ynd(NdF)
                pr0z = tp1*znd(NdA)+tp2*znd(NdB)+tp3*znd(NdC) &
                &     +tp4*znd(NdD)+tp5*znd(NdE)+tp6*znd(NdF)

                pr0nnx = tpEnx
                pr0nny = tpEny
                pr0nnz = tpEnz

                tp = dsqrt(pr0nnx**2 + pr0nny**2 + pr0nnz**2)
                pr0nnx = pr0nnx/tp
                pr0nny = pr0nny/tp
                pr0nnz = pr0nnz/tp

                IF (bdry_type == 'ex') THEN
                    IF (dmEorH == 'E') THEN
                        ztpEHx0 = tp1*exE2x_EM(NdA)+tp2*exE2x_EM(NdB) &
                        &        +tp3*exE2x_EM(NdC)+tp4*exE2x_EM(NdD) &
                        &        +tp5*exE2x_EM(NdE)+tp6*exE2x_EM(NdF)
                        ztpEHxdn0= tp1*exE2xdnn_EM(NdA)+tp2*exE2xdnn_EM(NdB) &
                        &         +tp3*exE2xdnn_EM(NdC)+tp4*exE2xdnn_EM(NdD) &
                        &         +tp5*exE2xdnn_EM(NdE)+tp6*exE2xdnn_EM(NdF)
                        ztpEHy0 = tp1*exE2y_EM(NdA)+tp2*exE2y_EM(NdB) &
                        &        +tp3*exE2y_EM(NdC)+tp4*exE2y_EM(NdD) &
                        &        +tp5*exE2y_EM(NdE)+tp6*exE2y_EM(NdF)
                        ztpEHydn0= tp1*exE2ydnn_EM(NdA)+tp2*exE2ydnn_EM(NdB) &
                        &         +tp3*exE2ydnn_EM(NdC)+tp4*exE2ydnn_EM(NdD) &
                        &         +tp5*exE2ydnn_EM(NdE)+tp6*exE2ydnn_EM(NdF)
                        ztpEHz0 = tp1*exE2z_EM(NdA)+tp2*exE2z_EM(NdB) &
                        &        +tp3*exE2z_EM(NdC)+tp4*exE2z_EM(NdD) &
                        &        +tp5*exE2z_EM(NdE)+tp6*exE2z_EM(NdF)
                        ztpEHzdn0= tp1*exE2zdnn_EM(NdA)+tp2*exE2zdnn_EM(NdB) &
                        &         +tp3*exE2zdnn_EM(NdC)+tp4*exE2zdnn_EM(NdD) &
                        &         +tp5*exE2zdnn_EM(NdE)+tp6*exE2zdnn_EM(NdF)
                    END IF
                    IF (dmEorH == 'H') THEN
                        ztpEHx0 = tp1*exH2x_EM(NdA)+tp2*exH2x_EM(NdB) &
                        &        +tp3*exH2x_EM(NdC)+tp4*exH2x_EM(NdD) &
                        &        +tp5*exH2x_EM(NdE)+tp6*exH2x_EM(NdF)
                        ztpEHxdn0= tp1*exH2xdnn_EM(NdA)+tp2*exH2xdnn_EM(NdB) &
                        &         +tp3*exH2xdnn_EM(NdC)+tp4*exH2xdnn_EM(NdD) &
                        &         +tp5*exH2xdnn_EM(NdE)+tp6*exH2xdnn_EM(NdF)
                        ztpEHy0 = tp1*exH2y_EM(NdA)+tp2*exH2y_EM(NdB) &
                        &        +tp3*exH2y_EM(NdC)+tp4*exH2y_EM(NdD) &
                        &        +tp5*exH2y_EM(NdE)+tp6*exH2y_EM(NdF)
                        ztpEHydn0= tp1*exH2ydnn_EM(NdA)+tp2*exH2ydnn_EM(NdB) &
                        &         +tp3*exH2ydnn_EM(NdC)+tp4*exH2ydnn_EM(NdD) &
                        &         +tp5*exH2ydnn_EM(NdE)+tp6*exH2ydnn_EM(NdF)
                        ztpEHz0 = tp1*exH2z_EM(NdA)+tp2*exH2z_EM(NdB) &
                        &        +tp3*exH2z_EM(NdC)+tp4*exH2z_EM(NdD) &
                        &        +tp5*exH2z_EM(NdE)+tp6*exH2z_EM(NdF)
                        ztpEHzdn0= tp1*exH2zdnn_EM(NdA)+tp2*exH2zdnn_EM(NdB) &
                        &         +tp3*exH2zdnn_EM(NdC)+tp4*exH2zdnn_EM(NdD) &
                        &         +tp5*exH2zdnn_EM(NdE)+tp6*exH2zdnn_EM(NdF)
                    END IF
                END IF
                IF (bdry_type == 'in') THEN
                    IF (dmEorH == 'E') THEN
                        ztpEHx0 = tp1*inE2x_EM(NdA)+tp2*inE2x_EM(NdB) &
                        &        +tp3*inE2x_EM(NdC)+tp4*inE2x_EM(NdD) &
                        &        +tp5*inE2x_EM(NdE)+tp6*inE2x_EM(NdF)
                        ztpEHxdn0= tp1*inE2xdnn_EM(NdA)+tp2*inE2xdnn_EM(NdB) &
                        &         +tp3*inE2xdnn_EM(NdC)+tp4*inE2xdnn_EM(NdD) &
                        &         +tp5*inE2xdnn_EM(NdE)+tp6*inE2xdnn_EM(NdF)
                        ztpEHy0 = tp1*inE2y_EM(NdA)+tp2*inE2y_EM(NdB) &
                        &        +tp3*inE2y_EM(NdC)+tp4*inE2y_EM(NdD) &
                        &        +tp5*inE2y_EM(NdE)+tp6*inE2y_EM(NdF)
                        ztpEHydn0= tp1*inE2ydnn_EM(NdA)+tp2*inE2ydnn_EM(NdB) &
                        &         +tp3*inE2ydnn_EM(NdC)+tp4*inE2ydnn_EM(NdD) &
                        &         +tp5*inE2ydnn_EM(NdE)+tp6*inE2ydnn_EM(NdF)
                        ztpEHz0 = tp1*inE2z_EM(NdA)+tp2*inE2z_EM(NdB) &
                        &        +tp3*inE2z_EM(NdC)+tp4*inE2z_EM(NdD) &
                        &        +tp5*inE2z_EM(NdE)+tp6*inE2z_EM(NdF)
                        ztpEHzdn0= tp1*inE2zdnn_EM(NdA)+tp2*inE2zdnn_EM(NdB) &
                        &         +tp3*inE2zdnn_EM(NdC)+tp4*inE2zdnn_EM(NdD) &
                        &         +tp5*inE2zdnn_EM(NdE)+tp6*inE2zdnn_EM(NdF)
                    END IF
                    IF (dmEorH == 'H') THEN
                        ztpEHx0 = tp1*inH2x_EM(NdA)+tp2*inH2x_EM(NdB) &
                        &        +tp3*inH2x_EM(NdC)+tp4*inH2x_EM(NdD) &
                        &        +tp5*inH2x_EM(NdE)+tp6*inH2x_EM(NdF)
                        ztpEHxdn0= tp1*inH2xdnn_EM(NdA)+tp2*inH2xdnn_EM(NdB) &
                        &         +tp3*inH2xdnn_EM(NdC)+tp4*inH2xdnn_EM(NdD) &
                        &         +tp5*inH2xdnn_EM(NdE)+tp6*inH2xdnn_EM(NdF)
                        ztpEHy0 = tp1*inH2y_EM(NdA)+tp2*inH2y_EM(NdB) &
                        &        +tp3*inH2y_EM(NdC)+tp4*inH2y_EM(NdD) &
                        &        +tp5*inH2y_EM(NdE)+tp6*inH2y_EM(NdF)
                        ztpEHydn0= tp1*inH2ydnn_EM(NdA)+tp2*inH2ydnn_EM(NdB) &
                        &         +tp3*inH2ydnn_EM(NdC)+tp4*inH2ydnn_EM(NdD) &
                        &         +tp5*inH2ydnn_EM(NdE)+tp6*inH2ydnn_EM(NdF)
                        ztpEHz0 = tp1*inH2z_EM(NdA)+tp2*inH2z_EM(NdB) &
                        &        +tp3*inH2z_EM(NdC)+tp4*inH2z_EM(NdD) &
                        &        +tp5*inH2z_EM(NdE)+tp6*inH2z_EM(NdF)
                        ztpEHzdn0= tp1*inH2zdnn_EM(NdA)+tp2*inH2zdnn_EM(NdB) &
                        &         +tp3*inH2zdnn_EM(NdC)+tp4*inH2zdnn_EM(NdD) &
                        &         +tp5*inH2zdnn_EM(NdE)+tp6*inH2zdnn_EM(NdF)
                    END IF
                END IF

            END IF

!-

        END IF

        IF (dmdomain == 0) THEN

            dmEHxdmn = ztpEHx0
            dmEHydmn = ztpEHy0
            dmEHzdmn = ztpEHz0

            DO ithprtl = 1, nmbrprtl

                IF (corelnkshell(ithprtl) == dmdomain) THEN

                    ztp = exK_EM

                    DO k = elstaID(ithprtl), elendID(ithprtl)

                        IF (MeshType == "L") THEN

                            CALL CalGHLnrEHDmnCOMP(dmEorH, ztp, k, 'ex', &
                            &       pr0x, pr0y, pr0z, pr0nnx, pr0nny, pr0nnz, &
                            &       dmDmnx, dmDmny, dmDmnz, &
                            &       ztpEHx0, ztpEHy0, ztpEHz0, &
                            &       ztpEHxdn0, ztpEHydn0, ztpEHzdn0, &
                            &       ztpEHx_elmnt, ztpEHy_elmnt, ztpEHz_elmnt)

                        END IF

                        IF (MeshType == "Q") THEN

                            CALL CalGHQdrEHDmnCOMP(dmEorH, ztp, k, 'ex', &
                            &       pr0x, pr0y, pr0z, pr0nnx, pr0nny, pr0nnz, &
                            &       dmDmnx, dmDmny, dmDmnz, &
                            &       ztpEHx0, ztpEHy0, ztpEHz0, &
                            &       ztpEHxdn0, ztpEHydn0, ztpEHzdn0, &
                            &       ztpEHx_elmnt, ztpEHy_elmnt, ztpEHz_elmnt)

                        END IF

                        IF (NrmlInOut(ithprtl) == 1) THEN
                            dmEHxdmn = dmEHxdmn + ztpEHx_elmnt
                            dmEHydmn = dmEHydmn + ztpEHy_elmnt
                            dmEHzdmn = dmEHzdmn + ztpEHz_elmnt
                        END IF
                        IF (NrmlInOut(ithprtl) ==-1) THEN
                            dmEHxdmn = dmEHxdmn - ztpEHx_elmnt
                            dmEHydmn = dmEHydmn - ztpEHy_elmnt
                            dmEHzdmn = dmEHzdmn - ztpEHz_elmnt
                        END IF

                    END DO

                END IF

            END DO

        ELSE

            tp=pr0nnx*(dmDmnx-pr0x)+pr0nny*(dmDmny-pr0y)+pr0nnz*(dmDmnz-pr0z)
            dmEHxdmn = ztpEHx0 + ztpEHxdn0*tp
            dmEHydmn = ztpEHy0 + ztpEHydn0*tp
            dmEHzdmn = ztpEHz0 + ztpEHzdn0*tp

            DO ithprtl = 1, nmbrprtl

                IF (ithprtl == dmdomain) THEN

                    ztp = ink_EM(dmdomain)

                    DO k = elstaID(ithprtl), elendID(ithprtl)

                        IF (MeshType == "L") THEN

                            CALL CalGHLnrEHDmnCOMP(dmEorH, ztp, k, 'in', &
                            &       pr0x, pr0y, pr0z, pr0nnx, pr0nny, pr0nnz, &
                            &       dmDmnx, dmDmny, dmDmnz, &
                            &       ztpEHx0, ztpEHy0, ztpEHz0, &
                            &       ztpEHxdn0, ztpEHydn0, ztpEHzdn0, &
                            &       ztpEHx_elmnt, ztpEHy_elmnt, ztpEHz_elmnt)

                        END IF

                        IF (MeshType == "Q") THEN

                            CALL CalGHQdrEHDmnCOMP(dmEorH, ztp, k, 'in', &
                            &       pr0x, pr0y, pr0z, pr0nnx, pr0nny, pr0nnz, &
                            &       dmDmnx, dmDmny, dmDmnz, &
                            &       ztpEHx0, ztpEHy0, ztpEHz0, &
                            &       ztpEHxdn0, ztpEHydn0, ztpEHzdn0, &
                            &       ztpEHx_elmnt, ztpEHy_elmnt, ztpEHz_elmnt)

                        END IF

                        IF (NrmlInOut(dmdomain) == 1) THEN
                            dmEHxdmn = dmEHxdmn - ztpEHx_elmnt
                            dmEHydmn = dmEHydmn - ztpEHy_elmnt
                            dmEHzdmn = dmEHzdmn - ztpEHz_elmnt
                        END IF
                        IF (NrmlInOut(dmdomain) ==-1) THEN
                            dmEHxdmn = dmEHxdmn + ztpEHx_elmnt
                            dmEHydmn = dmEHydmn + ztpEHy_elmnt
                            dmEHzdmn = dmEHzdmn + ztpEHz_elmnt
                        END IF

                    END DO

                END IF

                IF (corelnkshell(ithprtl) == dmdomain) THEN

                    ztp = ink_EM(dmdomain)

                    DO k = elstaID(ithprtl), elendID(ithprtl)

                        IF (MeshType == "L") THEN

                            CALL CalGHLnrEHDmnCOMP(dmEorH, ztp, k, 'ex', &
                            &       pr0x, pr0y, pr0z, pr0nnx, pr0nny, pr0nnz, &
                            &       dmDmnx, dmDmny, dmDmnz, &
                            &       ztpEHx0, ztpEHy0, ztpEHz0, &
                            &       ztpEHxdn0, ztpEHydn0, ztpEHzdn0, &
                            &       ztpEHx_elmnt, ztpEHy_elmnt, ztpEHz_elmnt)

                        END IF

                        IF (MeshType == "Q") THEN

                            CALL CalGHQdrEHDmnCOMP(dmEorH, ztp, k, 'ex', &
                            &       pr0x, pr0y, pr0z, pr0nnx, pr0nny, pr0nnz, &
                            &       dmDmnx, dmDmny, dmDmnz, &
                            &       ztpEHx0, ztpEHy0, ztpEHz0, &
                            &       ztpEHxdn0, ztpEHydn0, ztpEHzdn0, &
                            &       ztpEHx_elmnt, ztpEHy_elmnt, ztpEHz_elmnt)

                        END IF

                        IF (NrmlInOut(dmdomain) == 1) THEN
                            dmEHxdmn = dmEHxdmn - ztpEHx_elmnt
                            dmEHydmn = dmEHydmn - ztpEHy_elmnt
                            dmEHzdmn = dmEHzdmn - ztpEHz_elmnt
                        END IF
                        IF (NrmlInOut(dmdomain) ==-1) THEN
                            dmEHxdmn = dmEHxdmn + ztpEHx_elmnt
                            dmEHydmn = dmEHydmn + ztpEHy_elmnt
                            dmEHzdmn = dmEHzdmn + ztpEHz_elmnt
                        END IF

                    END DO

                END IF

            END DO

        END IF

    END SUBROUTINE

    ! Return one L3 element's contribution to the regularised field integral.

    SUBROUTINE CalGHLnrEHDmnCOMP(dmdmEorH, dmKwn, dmelmnid, dmbdry_type, &
    &                            dmp0x, dmp0y, dmp0z, dmp0nnx, dmp0nny, dmp0nnz, &
    &                            dmdmDmnx, dmdmDmny, dmdmDmnz, &
    &                            dmEHx0, dmEHy0, dmEHz0, dmEHxdn0, dmEHydn0, dmEHzdn0, &
    &                            dmEHx_elmnt, dmEHy_elmnt, dmEHz_elmnt)

        CHARACTER (LEN=1), INTENT(IN) :: dmdmEorH
        COMPLEX(KIND=KIND(1.0D0)), INTENT(IN) ::  dmKwn
        INTEGER, INTENT(IN) ::  dmelmnid
        CHARACTER (LEN=2), INTENT(IN) :: dmbdry_type
        DOUBLE PRECISION, INTENT(IN) ::  dmp0x, dmp0y, dmp0z, dmp0nnx, dmp0nny, dmp0nnz

        DOUBLE PRECISION, INTENT(IN) ::  dmdmDmnx, dmdmDmny, dmdmDmnz

        COMPLEX(KIND=KIND(1.0D0)), INTENT (IN) :: dmEHx0,dmEHy0,dmEHz0,dmEHxdn0,dmEHydn0,dmEHzdn0

        COMPLEX(KIND=KIND(1.0D0)), INTENT(OUT) :: dmEHx_elmnt, dmEHy_elmnt, dmEHz_elmnt

        DOUBLE PRECISION :: tpEnx, tpEny, tpEnz

        COMPLEX(KIND=KIND(1.0D0)) :: ztpGren,ztpdGdn,ztpGren0,ztpdGdn0, &
        &                            ztpGren_p,ztpdGdn_p,ztpGren0_p,ztpdGdn0_p

        COMPLEX(KIND=KIND(1.0D0)) :: ztpxsum, ztpysum, ztpzsum

        DOUBLE PRECISION :: tprx, tpry, tprz

        DOUBLE PRECISION :: tprr0x, tprr0y, tprr0z

        DOUBLE PRECISION :: mdl_r0r, over_r0r1, over_r0r2, over_r0r3, tpdp

        DOUBLE PRECISION :: tpn0xx0, tpn0En

        COMPLEX(KIND=KIND(1.0D0)) :: ztp11, ztp12, ztp13, ztp14, ztp15, ztp16
        COMPLEX(KIND=KIND(1.0D0)) :: ztp21, ztp22, ztp23, ztp24, ztp25, ztp26
        COMPLEX(KIND=KIND(1.0D0)) :: ztp1, ztp2

        COMPLEX(KIND=KIND(1.0D0)) :: ztpkwni, ztpExpkwnr

        INTEGER :: GLQi, icnt

        ztpkwni = ztponei*dmKwn

        ztpxsum = ztpzero
        ztpysum = ztpzero
        ztpzsum = ztpzero

        DO GLQi = 1, n_glqtr2d

            icnt = n_glqtr2d*(dmelmnid-1)+GLQi

            tprx = srcfmm_vec(1, icnt)
            tpry = srcfmm_vec(2, icnt)
            tprz = srcfmm_vec(3, icnt)

            tpEnx = srcfmm_nrm(1, icnt)
            tpEny = srcfmm_nrm(2, icnt)
            tpEnz = srcfmm_nrm(3, icnt)

            !Green's function and so on : ztpGren, ztpdGdn, ztpGren0, ztpdGdn0
            tprr0x = tprx - dmp0x
            tprr0y = tpry - dmp0y
            tprr0z = tprz - dmp0z
            mdl_r0r = DSQRT(tprr0x**2 + tprr0y**2 + tprr0z**2)
            over_r0r1 = 1.0d0/mdl_r0r
            over_r0r2 = over_r0r1*over_r0r1
            over_r0r3 = over_r0r2*over_r0r1
            tpdp = -(tpEnx*tprr0x + tpEny*tprr0y + tpEnz*tprr0z)
            ztpExpkwnr = CDEXP(ztpkwni*mdl_r0r)
            ztpGren = over_r0r1*ztpExpkwnr
            ztpdGdn = tpdp*(over_r0r1-ztpkwni)*ztpExpkwnr*over_r0r2
            ztpGren0 = DCMPLX(over_r0r1,0.0d0)
            ztpdGdn0 = DCMPLX(tpdp*over_r0r3,0.0d0)

            tpn0xx0 = dmp0nnx*tprr0x+dmp0nny*tprr0y+dmp0nnz*tprr0z
            tpn0En  = dmp0nnx*tpEnx+dmp0nny*tpEny+dmp0nnz*tpEnz

            !Green's function_p and so on : ztpGren_p, ztpdGdn_p, ztpGren0_p, ztpdGdn0_p
            tprr0x = tprx - dmdmDmnx
            tprr0y = tpry - dmdmDmny
            tprr0z = tprz - dmdmDmnz
            mdl_r0r = DSQRT(tprr0x**2 + tprr0y**2 + tprr0z**2)
            over_r0r1 = 1.0d0/mdl_r0r
            over_r0r2 = over_r0r1*over_r0r1
            over_r0r3 = over_r0r2*over_r0r1
            tpdp = -(tpEnx*tprr0x + tpEny*tprr0y + tpEnz*tprr0z)
            ztpExpkwnr = CDEXP(ztpkwni*mdl_r0r)
            ztpGren_p = over_r0r1*ztpExpkwnr
            ztpdGdn_p = tpdp*(over_r0r1-ztpkwni)*ztpExpkwnr*over_r0r2
            ztpGren0_p = DCMPLX(over_r0r1,0.0d0)
            ztpdGdn0_p = DCMPLX(tpdp*over_r0r3,0.0d0)

            ! Apply the non-singular boundary-integral correction.

            IF (dmbdry_type == 'ex') THEN
                IF (dmdmEorH == 'E') THEN
                    ztp11 = exE2xdnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = exE2xdnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = exE2xdnn_EM(elmntlnknd(dmelmnid,3))
                    ztp21 =-exE2x_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-exE2x_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-exE2x_EM(elmntlnknd(dmelmnid,3))
                END IF
                IF (dmdmEorH == 'H') THEN
                    ztp11 = exH2xdnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = exH2xdnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = exH2xdnn_EM(elmntlnknd(dmelmnid,3))
                    ztp21 =-exH2x_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-exH2x_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-exH2x_EM(elmntlnknd(dmelmnid,3))
                END IF
            END IF
            IF (dmbdry_type == 'in') THEN
                IF (dmdmEorH == 'E') THEN
                    ztp11 = inE2xdnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = inE2xdnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = inE2xdnn_EM(elmntlnknd(dmelmnid,3))
                    ztp21 =-inE2x_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-inE2x_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-inE2x_EM(elmntlnknd(dmelmnid,3))
                END IF
                IF (dmdmEorH == 'H') THEN
                    ztp11 = inH2xdnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = inH2xdnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = inH2xdnn_EM(elmntlnknd(dmelmnid,3))
                    ztp21 =-inH2x_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-inH2x_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-inH2x_EM(elmntlnknd(dmelmnid,3))
                END IF
            END IF
            ztp1 = srcfmm_wtnd(1,icnt)*ztp11 &
                & +srcfmm_wtnd(2,icnt)*ztp12 &
                & +srcfmm_wtnd(3,icnt)*ztp13
            ztp2 = srcfmm_wtnd(1,icnt)*ztp21 &
                & +srcfmm_wtnd(2,icnt)*ztp22 &
                & +srcfmm_wtnd(3,icnt)*ztp23
            ztpxsum = ztpxsum &
            &        +ztp2*(ztpdGdn_p-ztpdGdn) &
            &        +ztp1*(ztpGren_p-ztpGren) &
            &        +srcfmm_wght(icnt)*(  &
            &           +dmEHx0*(ztpdGdn0_p-ztpdGdn0) &
            &           +dmEHxdn0*tpn0xx0*(ztpdGdn0_p-ztpdGdn0) &
            &           -dmEHxdn0*tpn0En*(ztpGren0_p-ztpGren0)   )

            IF (dmbdry_type == 'ex') THEN
                IF (dmdmEorH == 'E') THEN
                    ztp11 = exE2ydnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = exE2ydnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = exE2ydnn_EM(elmntlnknd(dmelmnid,3))
                    ztp21 =-exE2y_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-exE2y_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-exE2y_EM(elmntlnknd(dmelmnid,3))
                END IF
                IF (dmdmEorH == 'H') THEN
                    ztp11 = exH2ydnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = exH2ydnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = exH2ydnn_EM(elmntlnknd(dmelmnid,3))
                    ztp21 =-exH2y_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-exH2y_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-exH2y_EM(elmntlnknd(dmelmnid,3))
                END IF
            END IF
            IF (dmbdry_type == 'in') THEN
                IF (dmdmEorH == 'E') THEN
                    ztp11 = inE2ydnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = inE2ydnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = inE2ydnn_EM(elmntlnknd(dmelmnid,3))
                    ztp21 =-inE2y_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-inE2y_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-inE2y_EM(elmntlnknd(dmelmnid,3))
                END IF
                IF (dmdmEorH == 'H') THEN
                    ztp11 = inH2ydnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = inH2ydnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = inH2ydnn_EM(elmntlnknd(dmelmnid,3))
                    ztp21 =-inH2y_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-inH2y_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-inH2y_EM(elmntlnknd(dmelmnid,3))
                END IF
            END IF
            ztp1 = srcfmm_wtnd(1,icnt)*ztp11 &
                & +srcfmm_wtnd(2,icnt)*ztp12 &
                & +srcfmm_wtnd(3,icnt)*ztp13
            ztp2 = srcfmm_wtnd(1,icnt)*ztp21 &
                & +srcfmm_wtnd(2,icnt)*ztp22 &
                & +srcfmm_wtnd(3,icnt)*ztp23
            ztpysum = ztpysum &
            &        +ztp2*(ztpdGdn_p-ztpdGdn) &
            &        +ztp1*(ztpGren_p-ztpGren) &
            &        +srcfmm_wght(icnt)*(  &
            &           +dmEHy0*(ztpdGdn0_p-ztpdGdn0) &
            &           +dmEHydn0*tpn0xx0*(ztpdGdn0_p-ztpdGdn0) &
            &           -dmEHydn0*tpn0En*(ztpGren0_p-ztpGren0)   )

            IF (dmbdry_type == 'ex') THEN
                IF (dmdmEorH == 'E') THEN
                    ztp11 = exE2zdnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = exE2zdnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = exE2zdnn_EM(elmntlnknd(dmelmnid,3))
                    ztp21 =-exE2z_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-exE2z_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-exE2z_EM(elmntlnknd(dmelmnid,3))
                END IF
                IF (dmdmEorH == 'H') THEN
                    ztp11 = exH2zdnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = exH2zdnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = exH2zdnn_EM(elmntlnknd(dmelmnid,3))
                    ztp21 =-exH2z_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-exH2z_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-exH2z_EM(elmntlnknd(dmelmnid,3))
                END IF
            END IF
            IF (dmbdry_type == 'in') THEN
                IF (dmdmEorH == 'E') THEN
                    ztp11 = inE2zdnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = inE2zdnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = inE2zdnn_EM(elmntlnknd(dmelmnid,3))
                    ztp21 =-inE2z_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-inE2z_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-inE2z_EM(elmntlnknd(dmelmnid,3))
                END IF
                IF (dmdmEorH == 'H') THEN
                    ztp11 = inH2zdnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = inH2zdnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = inH2zdnn_EM(elmntlnknd(dmelmnid,3))
                    ztp21 =-inH2z_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-inH2z_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-inH2z_EM(elmntlnknd(dmelmnid,3))
                END IF
            END IF
            ztp1 = srcfmm_wtnd(1,icnt)*ztp11 &
                & +srcfmm_wtnd(2,icnt)*ztp12 &
                & +srcfmm_wtnd(3,icnt)*ztp13
            ztp2 = srcfmm_wtnd(1,icnt)*ztp21 &
                & +srcfmm_wtnd(2,icnt)*ztp22 &
                & +srcfmm_wtnd(3,icnt)*ztp23
            ztpzsum = ztpzsum &
            &        +ztp2*(ztpdGdn_p-ztpdGdn) &
            &        +ztp1*(ztpGren_p-ztpGren) &
            &        +srcfmm_wght(icnt)*(  &
            &           +dmEHz0*(ztpdGdn0_p-ztpdGdn0) &
            &           +dmEHzdn0*tpn0xx0*(ztpdGdn0_p-ztpdGdn0) &
            &           -dmEHzdn0*tpn0En*(ztpGren0_p-ztpGren0)   )

        END DO

        dmEHx_elmnt = ztpxsum / (4.0d0*pai)
        dmEHy_elmnt = ztpysum / (4.0d0*pai)
        dmEHz_elmnt = ztpzsum / (4.0d0*pai)


    END SUBROUTINE

    ! Return one Q6 element's contribution to the regularised field integral.
    ! The accumulated element vector includes the conventional 1/(4*pi) factor.

    SUBROUTINE CalGHQdrEHDmnCOMP(dmdmEorH, dmKwn, dmelmnid, dmbdry_type, &
    &                            dmp0x, dmp0y, dmp0z, dmp0nnx, dmp0nny, dmp0nnz, &
    &                            dmdmDmnx, dmdmDmny, dmdmDmnz, &
    &                            dmEHx0, dmEHy0, dmEHz0, dmEHxdn0, dmEHydn0, dmEHzdn0, &
    &                            dmEHx_elmnt, dmEHy_elmnt, dmEHz_elmnt)

        CHARACTER (LEN=1), INTENT(IN) :: dmdmEorH
        COMPLEX(KIND=KIND(1.0D0)), INTENT(IN) ::  dmKwn
        INTEGER, INTENT(IN) ::  dmelmnid
        CHARACTER (LEN=2), INTENT(IN) :: dmbdry_type
        DOUBLE PRECISION, INTENT(IN) ::  dmp0x, dmp0y, dmp0z, dmp0nnx, dmp0nny, dmp0nnz

        DOUBLE PRECISION, INTENT(IN) ::  dmdmDmnx, dmdmDmny, dmdmDmnz

        COMPLEX(KIND=KIND(1.0D0)), INTENT (IN) :: dmEHx0,dmEHy0,dmEHz0,dmEHxdn0,dmEHydn0,dmEHzdn0

        COMPLEX(KIND=KIND(1.0D0)), INTENT(OUT) :: dmEHx_elmnt, dmEHy_elmnt, dmEHz_elmnt

        COMPLEX(KIND=KIND(1.0D0)) :: ztpGren,ztpdGdn,ztpGren0,ztpdGdn0, &
        &                            ztpGren_p,ztpdGdn_p,ztpGren0_p,ztpdGdn0_p

        COMPLEX(KIND=KIND(1.0D0)) :: ztpxsum, ztpysum, ztpzsum

        DOUBLE PRECISION :: tprx, tpry, tprz

        DOUBLE PRECISION :: tprr0x, tprr0y, tprr0z

        DOUBLE PRECISION :: mdl_r0r, over_r0r1, over_r0r2, over_r0r3, tpdp

        DOUBLE PRECISION :: tpn0xx0, tpn0En

        COMPLEX(KIND=KIND(1.0D0)) :: ztpx, tpExdn, tpEy, tpEydn, tpEz, tpEzdn

        COMPLEX(KIND=KIND(1.0D0)) :: ztp11, ztp12, ztp13, ztp14, ztp15, ztp16
        COMPLEX(KIND=KIND(1.0D0)) :: ztp21, ztp22, ztp23, ztp24, ztp25, ztp26
        COMPLEX(KIND=KIND(1.0D0)) :: ztp1, ztp2

        COMPLEX(KIND=KIND(1.0D0)) :: ztpkwni, ztpExpkwnr

        DOUBLE PRECISION :: tpendx1,tpendx2,tpendx3,tpendx4,tpendx5,tpendx6
        DOUBLE PRECISION :: tpendy1,tpendy2,tpendy3,tpendy4,tpendy5,tpendy6
        DOUBLE PRECISION :: tpendz1,tpendz2,tpendz3,tpendz4,tpendz5,tpendz6
        DOUBLE PRECISION :: drx_deps,drx_dyet,dry_deps,dry_dyet,drz_deps,drz_dyet
        DOUBLE PRECISION :: tpEnx, tpEny, tpEnz, JcbDtmn
        DOUBLE PRECISION :: tpxieta

        INTEGER :: GLQi, icnt

        ztpkwni = ztponei*dmKwn

        ztpxsum = ztpzero
        ztpysum = ztpzero
        ztpzsum = ztpzero

        DO GLQi = 1, n_glqtr2d

            icnt = n_glqtr2d*(dmelmnid-1)+GLQi

            tprx = srcfmm_vec(1, icnt)
            tpry = srcfmm_vec(2, icnt)
            tprz = srcfmm_vec(3, icnt)

            tpEnx = srcfmm_nrm(1, icnt)
            tpEny = srcfmm_nrm(2, icnt)
            tpEnz = srcfmm_nrm(3, icnt)

            !Green's function and so on : ztpGren, ztpdGdn, ztpGren0, ztpdGdn0
            tprr0x = tprx - dmp0x
            tprr0y = tpry - dmp0y
            tprr0z = tprz - dmp0z
            mdl_r0r = DSQRT(tprr0x**2 + tprr0y**2 + tprr0z**2)
            over_r0r1 = 1.0d0/mdl_r0r
            over_r0r2 = over_r0r1*over_r0r1
            over_r0r3 = over_r0r2*over_r0r1
            tpdp = -(tpEnx*tprr0x + tpEny*tprr0y + tpEnz*tprr0z)
            ztpExpkwnr = CDEXP(ztpkwni*mdl_r0r)
            ztpGren = over_r0r1*ztpExpkwnr
            ztpdGdn = tpdp*(over_r0r1-ztpkwni)*ztpExpkwnr*over_r0r2
            ztpGren0 = DCMPLX(over_r0r1,0.0d0)
            ztpdGdn0 = DCMPLX(tpdp*over_r0r3,0.0d0)

            tpn0xx0 = dmp0nnx*tprr0x+dmp0nny*tprr0y+dmp0nnz*tprr0z
            tpn0En  = dmp0nnx*tpEnx+dmp0nny*tpEny+dmp0nnz*tpEnz

            !Green's function_p and so on : ztpGren_p, ztpdGdn_p, ztpGren0_p, ztpdGdn0_p
            tprr0x = tprx - dmdmDmnx
            tprr0y = tpry - dmdmDmny
            tprr0z = tprz - dmdmDmnz
            mdl_r0r = DSQRT(tprr0x**2 + tprr0y**2 + tprr0z**2)
            over_r0r1 = 1.0d0/mdl_r0r
            over_r0r2 = over_r0r1*over_r0r1
            over_r0r3 = over_r0r2*over_r0r1
            tpdp = -(tpEnx*tprr0x + tpEny*tprr0y + tpEnz*tprr0z)
            ztpExpkwnr = CDEXP(ztpkwni*mdl_r0r)
            ztpGren_p = over_r0r1*ztpExpkwnr
            ztpdGdn_p = tpdp*(over_r0r1-ztpkwni)*ztpExpkwnr*over_r0r2
            ztpGren0_p = DCMPLX(over_r0r1,0.0d0)
            ztpdGdn0_p = DCMPLX(tpdp*over_r0r3,0.0d0)

            ! Apply the non-singular boundary-integral correction.

            IF (dmbdry_type == 'ex') THEN
                IF (dmdmEorH == 'E') THEN
                    ztp11 = exE2xdnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = exE2xdnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = exE2xdnn_EM(elmntlnknd(dmelmnid,3))
                    ztp14 = exE2xdnn_EM(elmntlnknd(dmelmnid,4))
                    ztp15 = exE2xdnn_EM(elmntlnknd(dmelmnid,5))
                    ztp16 = exE2xdnn_EM(elmntlnknd(dmelmnid,6))
                    ztp21 =-exE2x_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-exE2x_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-exE2x_EM(elmntlnknd(dmelmnid,3))
                    ztp24 =-exE2x_EM(elmntlnknd(dmelmnid,4))
                    ztp25 =-exE2x_EM(elmntlnknd(dmelmnid,5))
                    ztp26 =-exE2x_EM(elmntlnknd(dmelmnid,6))
                END IF
                IF (dmdmEorH == 'H') THEN
                    ztp11 = exH2xdnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = exH2xdnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = exH2xdnn_EM(elmntlnknd(dmelmnid,3))
                    ztp14 = exH2xdnn_EM(elmntlnknd(dmelmnid,4))
                    ztp15 = exH2xdnn_EM(elmntlnknd(dmelmnid,5))
                    ztp16 = exH2xdnn_EM(elmntlnknd(dmelmnid,6))
                    ztp21 =-exH2x_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-exH2x_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-exH2x_EM(elmntlnknd(dmelmnid,3))
                    ztp24 =-exH2x_EM(elmntlnknd(dmelmnid,4))
                    ztp25 =-exH2x_EM(elmntlnknd(dmelmnid,5))
                    ztp26 =-exH2x_EM(elmntlnknd(dmelmnid,6))
                END IF
            END IF
            IF (dmbdry_type == 'in') THEN
                IF (dmdmEorH == 'E') THEN
                    ztp11 = inE2xdnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = inE2xdnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = inE2xdnn_EM(elmntlnknd(dmelmnid,3))
                    ztp14 = inE2xdnn_EM(elmntlnknd(dmelmnid,4))
                    ztp15 = inE2xdnn_EM(elmntlnknd(dmelmnid,5))
                    ztp16 = inE2xdnn_EM(elmntlnknd(dmelmnid,6))
                    ztp21 =-inE2x_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-inE2x_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-inE2x_EM(elmntlnknd(dmelmnid,3))
                    ztp24 =-inE2x_EM(elmntlnknd(dmelmnid,4))
                    ztp25 =-inE2x_EM(elmntlnknd(dmelmnid,5))
                    ztp26 =-inE2x_EM(elmntlnknd(dmelmnid,6))
                END IF
                IF (dmdmEorH == 'H') THEN
                    ztp11 = inH2xdnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = inH2xdnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = inH2xdnn_EM(elmntlnknd(dmelmnid,3))
                    ztp14 = inH2xdnn_EM(elmntlnknd(dmelmnid,4))
                    ztp15 = inH2xdnn_EM(elmntlnknd(dmelmnid,5))
                    ztp16 = inH2xdnn_EM(elmntlnknd(dmelmnid,6))
                    ztp21 =-inH2x_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-inH2x_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-inH2x_EM(elmntlnknd(dmelmnid,3))
                    ztp24 =-inH2x_EM(elmntlnknd(dmelmnid,4))
                    ztp25 =-inH2x_EM(elmntlnknd(dmelmnid,5))
                    ztp26 =-inH2x_EM(elmntlnknd(dmelmnid,6))
                END IF
            END IF
            ztp1 = srcfmm_wtnd(1,icnt)*ztp11 &
                & +srcfmm_wtnd(2,icnt)*ztp12 &
                & +srcfmm_wtnd(3,icnt)*ztp13 &
                & +srcfmm_wtnd(4,icnt)*ztp14 &
                & +srcfmm_wtnd(5,icnt)*ztp15 &
                & +srcfmm_wtnd(6,icnt)*ztp16
            ztp2 = srcfmm_wtnd(1,icnt)*ztp21 &
                & +srcfmm_wtnd(2,icnt)*ztp22 &
                & +srcfmm_wtnd(3,icnt)*ztp23 &
                & +srcfmm_wtnd(4,icnt)*ztp24 &
                & +srcfmm_wtnd(5,icnt)*ztp25 &
                & +srcfmm_wtnd(6,icnt)*ztp26
            ztpxsum = ztpxsum &
            &        +ztp2*(ztpdGdn_p-ztpdGdn) &
            &        +ztp1*(ztpGren_p-ztpGren) &
            &        +srcfmm_wght(icnt)*(  &
            &           +dmEHx0*(ztpdGdn0_p-ztpdGdn0) &
            &           +dmEHxdn0*tpn0xx0*(ztpdGdn0_p-ztpdGdn0) &
            &           -dmEHxdn0*tpn0En*(ztpGren0_p-ztpGren0)   )

            IF (dmbdry_type == 'ex') THEN
                IF (dmdmEorH == 'E') THEN
                    ztp11 = exE2ydnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = exE2ydnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = exE2ydnn_EM(elmntlnknd(dmelmnid,3))
                    ztp14 = exE2ydnn_EM(elmntlnknd(dmelmnid,4))
                    ztp15 = exE2ydnn_EM(elmntlnknd(dmelmnid,5))
                    ztp16 = exE2ydnn_EM(elmntlnknd(dmelmnid,6))
                    ztp21 =-exE2y_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-exE2y_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-exE2y_EM(elmntlnknd(dmelmnid,3))
                    ztp24 =-exE2y_EM(elmntlnknd(dmelmnid,4))
                    ztp25 =-exE2y_EM(elmntlnknd(dmelmnid,5))
                    ztp26 =-exE2y_EM(elmntlnknd(dmelmnid,6))
                END IF
                IF (dmdmEorH == 'H') THEN
                    ztp11 = exH2ydnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = exH2ydnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = exH2ydnn_EM(elmntlnknd(dmelmnid,3))
                    ztp14 = exH2ydnn_EM(elmntlnknd(dmelmnid,4))
                    ztp15 = exH2ydnn_EM(elmntlnknd(dmelmnid,5))
                    ztp16 = exH2ydnn_EM(elmntlnknd(dmelmnid,6))
                    ztp21 =-exH2y_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-exH2y_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-exH2y_EM(elmntlnknd(dmelmnid,3))
                    ztp24 =-exH2y_EM(elmntlnknd(dmelmnid,4))
                    ztp25 =-exH2y_EM(elmntlnknd(dmelmnid,5))
                    ztp26 =-exH2y_EM(elmntlnknd(dmelmnid,6))
                END IF
            END IF
            IF (dmbdry_type == 'in') THEN
                IF (dmdmEorH == 'E') THEN
                    ztp11 = inE2ydnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = inE2ydnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = inE2ydnn_EM(elmntlnknd(dmelmnid,3))
                    ztp14 = inE2ydnn_EM(elmntlnknd(dmelmnid,4))
                    ztp15 = inE2ydnn_EM(elmntlnknd(dmelmnid,5))
                    ztp16 = inE2ydnn_EM(elmntlnknd(dmelmnid,6))
                    ztp21 =-inE2y_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-inE2y_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-inE2y_EM(elmntlnknd(dmelmnid,3))
                    ztp24 =-inE2y_EM(elmntlnknd(dmelmnid,4))
                    ztp25 =-inE2y_EM(elmntlnknd(dmelmnid,5))
                    ztp26 =-inE2y_EM(elmntlnknd(dmelmnid,6))
                END IF
                IF (dmdmEorH == 'H') THEN
                    ztp11 = inH2ydnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = inH2ydnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = inH2ydnn_EM(elmntlnknd(dmelmnid,3))
                    ztp14 = inH2ydnn_EM(elmntlnknd(dmelmnid,4))
                    ztp15 = inH2ydnn_EM(elmntlnknd(dmelmnid,5))
                    ztp16 = inH2ydnn_EM(elmntlnknd(dmelmnid,6))
                    ztp21 =-inH2y_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-inH2y_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-inH2y_EM(elmntlnknd(dmelmnid,3))
                    ztp24 =-inH2y_EM(elmntlnknd(dmelmnid,4))
                    ztp25 =-inH2y_EM(elmntlnknd(dmelmnid,5))
                    ztp26 =-inH2y_EM(elmntlnknd(dmelmnid,6))
                END IF
            END IF
            ztp1 = srcfmm_wtnd(1,icnt)*ztp11 &
                & +srcfmm_wtnd(2,icnt)*ztp12 &
                & +srcfmm_wtnd(3,icnt)*ztp13 &
                & +srcfmm_wtnd(4,icnt)*ztp14 &
                & +srcfmm_wtnd(5,icnt)*ztp15 &
                & +srcfmm_wtnd(6,icnt)*ztp16
            ztp2 = srcfmm_wtnd(1,icnt)*ztp21 &
                & +srcfmm_wtnd(2,icnt)*ztp22 &
                & +srcfmm_wtnd(3,icnt)*ztp23 &
                & +srcfmm_wtnd(4,icnt)*ztp24 &
                & +srcfmm_wtnd(5,icnt)*ztp25 &
                & +srcfmm_wtnd(6,icnt)*ztp26
            ztpysum = ztpysum &
            &        +ztp2*(ztpdGdn_p-ztpdGdn) &
            &        +ztp1*(ztpGren_p-ztpGren) &
            &        +srcfmm_wght(icnt)*(  &
            &           +dmEHy0*(ztpdGdn0_p-ztpdGdn0) &
            &           +dmEHydn0*tpn0xx0*(ztpdGdn0_p-ztpdGdn0) &
            &           -dmEHydn0*tpn0En*(ztpGren0_p-ztpGren0)   )

            IF (dmbdry_type == 'ex') THEN
                IF (dmdmEorH == 'E') THEN
                    ztp11 = exE2zdnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = exE2zdnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = exE2zdnn_EM(elmntlnknd(dmelmnid,3))
                    ztp14 = exE2zdnn_EM(elmntlnknd(dmelmnid,4))
                    ztp15 = exE2zdnn_EM(elmntlnknd(dmelmnid,5))
                    ztp16 = exE2zdnn_EM(elmntlnknd(dmelmnid,6))
                    ztp21 =-exE2z_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-exE2z_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-exE2z_EM(elmntlnknd(dmelmnid,3))
                    ztp24 =-exE2z_EM(elmntlnknd(dmelmnid,4))
                    ztp25 =-exE2z_EM(elmntlnknd(dmelmnid,5))
                    ztp26 =-exE2z_EM(elmntlnknd(dmelmnid,6))
                END IF
                IF (dmdmEorH == 'H') THEN
                    ztp11 = exH2zdnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = exH2zdnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = exH2zdnn_EM(elmntlnknd(dmelmnid,3))
                    ztp14 = exH2zdnn_EM(elmntlnknd(dmelmnid,4))
                    ztp15 = exH2zdnn_EM(elmntlnknd(dmelmnid,5))
                    ztp16 = exH2zdnn_EM(elmntlnknd(dmelmnid,6))
                    ztp21 =-exH2z_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-exH2z_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-exH2z_EM(elmntlnknd(dmelmnid,3))
                    ztp24 =-exH2z_EM(elmntlnknd(dmelmnid,4))
                    ztp25 =-exH2z_EM(elmntlnknd(dmelmnid,5))
                    ztp26 =-exH2z_EM(elmntlnknd(dmelmnid,6))
                END IF
            END IF
            IF (dmbdry_type == 'in') THEN
                IF (dmdmEorH == 'E') THEN
                    ztp11 = inE2zdnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = inE2zdnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = inE2zdnn_EM(elmntlnknd(dmelmnid,3))
                    ztp14 = inE2zdnn_EM(elmntlnknd(dmelmnid,4))
                    ztp15 = inE2zdnn_EM(elmntlnknd(dmelmnid,5))
                    ztp16 = inE2zdnn_EM(elmntlnknd(dmelmnid,6))
                    ztp21 =-inE2z_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-inE2z_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-inE2z_EM(elmntlnknd(dmelmnid,3))
                    ztp24 =-inE2z_EM(elmntlnknd(dmelmnid,4))
                    ztp25 =-inE2z_EM(elmntlnknd(dmelmnid,5))
                    ztp26 =-inE2z_EM(elmntlnknd(dmelmnid,6))
                END IF
                IF (dmdmEorH == 'H') THEN
                    ztp11 = inH2zdnn_EM(elmntlnknd(dmelmnid,1))
                    ztp12 = inH2zdnn_EM(elmntlnknd(dmelmnid,2))
                    ztp13 = inH2zdnn_EM(elmntlnknd(dmelmnid,3))
                    ztp14 = inH2zdnn_EM(elmntlnknd(dmelmnid,4))
                    ztp15 = inH2zdnn_EM(elmntlnknd(dmelmnid,5))
                    ztp16 = inH2zdnn_EM(elmntlnknd(dmelmnid,6))
                    ztp21 =-inH2z_EM(elmntlnknd(dmelmnid,1))
                    ztp22 =-inH2z_EM(elmntlnknd(dmelmnid,2))
                    ztp23 =-inH2z_EM(elmntlnknd(dmelmnid,3))
                    ztp24 =-inH2z_EM(elmntlnknd(dmelmnid,4))
                    ztp25 =-inH2z_EM(elmntlnknd(dmelmnid,5))
                    ztp26 =-inH2z_EM(elmntlnknd(dmelmnid,6))
                END IF
            END IF
            ztp1 = srcfmm_wtnd(1,icnt)*ztp11 &
                & +srcfmm_wtnd(2,icnt)*ztp12 &
                & +srcfmm_wtnd(3,icnt)*ztp13 &
                & +srcfmm_wtnd(4,icnt)*ztp14 &
                & +srcfmm_wtnd(5,icnt)*ztp15 &
                & +srcfmm_wtnd(6,icnt)*ztp16
            ztp2 = srcfmm_wtnd(1,icnt)*ztp21 &
                & +srcfmm_wtnd(2,icnt)*ztp22 &
                & +srcfmm_wtnd(3,icnt)*ztp23 &
                & +srcfmm_wtnd(4,icnt)*ztp24 &
                & +srcfmm_wtnd(5,icnt)*ztp25 &
                & +srcfmm_wtnd(6,icnt)*ztp26
            ztpzsum = ztpzsum &
            &        +ztp2*(ztpdGdn_p-ztpdGdn) &
            &        +ztp1*(ztpGren_p-ztpGren) &
            &        +srcfmm_wght(icnt)*(  &
            &           +dmEHz0*(ztpdGdn0_p-ztpdGdn0) &
            &           +dmEHzdn0*tpn0xx0*(ztpdGdn0_p-ztpdGdn0) &
            &           -dmEHzdn0*tpn0En*(ztpGren0_p-ztpGren0)   )

        END DO

        dmEHx_elmnt = ztpxsum / (4.0d0*pai)
        dmEHy_elmnt = ztpysum / (4.0d0*pai)
        dmEHz_elmnt = ztpzsum / (4.0d0*pai)

    END SUBROUTINE



END MODULE EM_DmnCal
