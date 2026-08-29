! SPDX-FileCopyrightText: 2026 Qiang Sun
! SPDX-License-Identifier: BSD-3-Clause

! Construct nodal normals, local tangent frames, surface derivatives, and
! curvature data from the oriented triangular mesh. Nodal quantities are
! assembled from incident elements; NrmlInOut fixes the physical normal side.
! Q6 cases use the quadratic geometry and derivative reconstruction required
! by the direct-field boundary equations.
!
MODULE Geom_NormVec

    USE omp_lib

    USE Pre_Constants
    USE Pre_csvformat

    USE Geom_GlobalData

    IMPLICIT NONE

    INTEGER, PRIVATE :: dt_option = 2

    CONTAINS

    ! Compute nodal normals and differential geometry for an L3 mesh.

    SUBROUTINE Getndnrml

        INTEGER :: ithprtl, i, j, k, ii, ij, ik, jj, kk, itmp, elmntA, ndA, ndB, ndC
        INTEGER :: ndst,nded,elst,eled,id_tp

        DOUBLE PRECISION :: tpvctACx, tpvctACy, tpvctACz, &
                         &  tpvctBCx, tpvctBCy, tpvctBCz, &
                         &  tpelmntnrmlx, tpelmntnrmly, tpelmntnrmlz, &
                         &  tpelmntareanrml

        DOUBLE PRECISION :: tpsumndnnx, tpsumndnny, tpsumndnnz, tpsumndnnmdl
        DOUBLE PRECISION :: tpnx, tpny, tpnz, tparea

        DOUBLE PRECISION ::  dmnx, dmny, dmnz, tplmn1, tplmn2, tpll1, tpll2, tpll3, &
                            &tpmm1, tpmm2, tpmm3, tpnn1, tpnn2, tpnn3

        DOUBLE PRECISION :: dtetax,dtetay,dtetaz,dleta,dt_xix,dt_xiy,dt_xiz,dt_xi,&
                           &dteta1,dteta2,dt_xi1,dt_xi2,dfac,dl_xi
        INTEGER :: slfNdA,slfNdB,slfNdC

        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: &
        &   cbcmtrxA,cbcmtrxAA,cbcmtrxAArcd,cbcmtrxAA_inv,cbcmtrxEnd,cbcmtrxA_tr
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: &
        &   cbcmtrxB,cbcmtrxW,cbcmtrxBB,cbcmtrxXX,cbcmtrxS
        INTEGER :: slv_ipiv(6),slv_D,slv_info
        DOUBLE PRECISION :: tp, tpcbcx, tpcbcy, tpcbcz, lcbcx, lcbcy, lcbcz

        INTEGER ::  GLQi, icnt
        DOUBLE PRECISION :: tp1,tp2,tp3,tp4,tp5,tp6,tp_a,tp_b,tp_g,tpxi,tpet
        DOUBLE PRECISION :: tprx,tpry,tprz,tpEnx,tpEny,tpEnz,JcbDtmn
        DOUBLE PRECISION :: tpendx1,tpendx2,tpendx3,tpendx4,tpendx5,tpendx6
        DOUBLE PRECISION :: tpendy1,tpendy2,tpendy3,tpendy4,tpendy5,tpendy6
        DOUBLE PRECISION :: tpendz1,tpendz2,tpendz3,tpendz4,tpendz5,tpendz6
        DOUBLE PRECISION :: drx_deps,drx_dyet,dry_deps,dry_dyet,drz_deps,drz_dyet
        DOUBLE PRECISION :: dphi_deps1,dphi_deps2,dphi_deps3,&
        &                   dphi_deps4,dphi_deps5,dphi_deps6
        DOUBLE PRECISION :: dphi_dyet1,dphi_dyet2,dphi_dyet3,&
        &                   dphi_dyet4,dphi_dyet5,dphi_dyet6
        DOUBLE PRECISION :: h_eps,h_yet,tpkappa,t_teps,t_tyet

        DOUBLE PRECISION :: tpAx, tpAy, tpAz, tpBx, tpBy, tpBz, tpCx, tpCy, tpCz

!==============
!   Element area calculation

!$OMP PARALLEL PRIVATE(k,ndA,ndB,ndC) &
!$OMP & PRIVATE(tpvctACx,tpvctACy,tpvctACz) &
!$OMP & PRIVATE(tpvctBCx,tpvctBCy,tpvctBCz) &
!$OMP & PRIVATE(tpelmntnrmlx,tpelmntnrmly,tpelmntnrmlz) &
!$OMP & PRIVATE(tpelmntareanrml)
!$OMP DO
        DO k = 1, ttlnmbrelmnt

            ndA = elmntlnknd(k,1)
            ndB = elmntlnknd(k,2)
            ndC = elmntlnknd(k,3)

            tpvctACx = xnd(ndA) - xnd(ndC)
            tpvctACy = ynd(ndA) - ynd(ndC)
            tpvctACz = znd(ndA) - znd(ndC)

            tpvctBCx = xnd(ndB) - xnd(ndC)
            tpvctBCy = ynd(ndB) - ynd(ndC)
            tpvctBCz = znd(ndB) - znd(ndC)

            tpelmntnrmlx = tpvctACy*tpvctBCz - tpvctACz*tpvctBCy
            tpelmntnrmly = tpvctACz*tpvctBCx - tpvctACx*tpvctBCz
            tpelmntnrmlz = tpvctACx*tpvctBCy - tpvctACy*tpvctBCx
            tpelmntareanrml = DSQRT(  tpelmntnrmlx**2 &
            &                       + tpelmntnrmly**2 &
            &                       + tpelmntnrmlz**2   )

            elmntarea(k) = 0.5d0 * tpelmntareanrml

            nnxelmnt(k) = tpelmntnrmlx / tpelmntareanrml
            nnyelmnt(k) = tpelmntnrmly / tpelmntareanrml
            nnzelmnt(k) = tpelmntnrmlz / tpelmntareanrml

        END DO
!$OMP END DO
!$OMP END PARALLEL

!==============

!==============
!   Gauss point weight, nodal position, nodal normal vector on each element

!$OMP PARALLEL PRIVATE(k,GLQi,icnt) &
!$OMP & PRIVATE(tp1,tp2,tp3,tpxi,tpet)&
!$OMP & PRIVATE(tprx,tpry,tprz,tpEnx,tpEny,tpEnz,JcbDtmn) &
!$OMP & PRIVATE(tpendx1,tpendx2,tpendx3) &
!$OMP & PRIVATE(tpendy1,tpendy2,tpendy3) &
!$OMP & PRIVATE(tpendz1,tpendz2,tpendz3)
!$OMP DO
        DO k = 1, ttlnmbrelmnt

            tpendx1 = xnd(elmntlnknd(k,1))
            tpendx2 = xnd(elmntlnknd(k,2))
            tpendx3 = xnd(elmntlnknd(k,3))

            tpendy1 = ynd(elmntlnknd(k,1))
            tpendy2 = ynd(elmntlnknd(k,2))
            tpendy3 = ynd(elmntlnknd(k,3))

            tpendz1 = znd(elmntlnknd(k,1))
            tpendz2 = znd(elmntlnknd(k,2))
            tpendz3 = znd(elmntlnknd(k,3))

            DO GLQi = 1, n_glqtr2d

                icnt = n_glqtr2d*(k-1)+GLQi

                tpxi = xg_glqtr2d(GLQi)
                tpet = yg_glqtr2d(GLQi)

                JcbDtmn = 2.0d0*elmntarea(k)

                tpEnx = nnxelmnt(k)
                tpEny = nnyelmnt(k)
                tpEnz = nnzelmnt(k)

                tp2 = tpxi
                tp3 = tpet
                tp1 = 1.0d0 - tp2 - tp3

                tprx = tp1*tpendx1+tp2*tpendx2+tp3*tpendx3
                tpry = tp1*tpendy1+tp2*tpendy2+tp3*tpendy3
                tprz = tp1*tpendz1+tp2*tpendz2+tp3*tpendz3

                srcfmm_vec(1, icnt) = tprx
                srcfmm_vec(2, icnt) = tpry
                srcfmm_vec(3, icnt) = tprz

                srcfmm_nrm(1, icnt) = tpEnx
                srcfmm_nrm(2, icnt) = tpEny
                srcfmm_nrm(3, icnt) = tpEnz

                srcfmm_wght(icnt) = wg_glqtr2d(GLQi)*JcbDtmn

                srcfmm_wtnd(1,icnt) = wg_glqtr2d(GLQi)*JcbDtmn*tp1
                srcfmm_wtnd(2,icnt) = wg_glqtr2d(GLQi)*JcbDtmn*tp2
                srcfmm_wtnd(3,icnt) = wg_glqtr2d(GLQi)*JcbDtmn*tp3

            END DO

        END DO
!$OMP END DO
!$OMP END PARALLEL

!==============

!==============
!   Normal, tangential vector calculation
!   Normal is based on the facet elements connecting to the node of interest
!   Weights are chosen as the combination of the element area
!   and how far the element is relative to the node
!   combine "MAX, N. Journal of Graphics Tools, Vol. 4, No. 2." &
!   "Chen, Wu Computer Aided Geometric Design 21 (2004) 447–458"
!   as: Area_j / (|g_j - v|^2) where g_j is the center of the triangle

!$OMP PARALLEL PRIVATE(i,k,elmntA,ndA,ndB,ndC) &
!$OMP & PRIVATE(tpsumndnnx,tpsumndnny,tpsumndnnz)&
!$OMP & PRIVATE(tpvctACx,tpvctACy,tpvctACz) &
!$OMP & PRIVATE(tpvctBCx,tpvctBCy,tpvctBCz) &
!$OMP & PRIVATE(tpelmntnrmlx,tpelmntnrmly,tpelmntnrmlz) &
!$OMP & PRIVATE(tpelmntareanrml) &
!$OMP & PRIVATE(tpnx,tpny,tpnz,tparea,tpsumndnnmdl) &
!$OMP & PRIVATE(tpll1,tpmm1,tpnn1,tplmn1,tpll2,tpmm2,tpnn2,tplmn2)
!$OMP DO
        DO i = 1, ttlnmbrnd

            tpsumndnnx = 0.0d0
            tpsumndnny = 0.0d0
            tpsumndnnz = 0.0d0

            DO k = 1, mxnmbrndlnkelmnt

                IF (ndlnkelmnt(i,k) /= 0) THEN

                    elmntA = ndlnkelmnt(i,k)

                    IF (i == elmntlnknd(elmntA,1)) THEN
                        ndC = elmntlnknd(elmntA,1)
                        ndA = elmntlnknd(elmntA,2)
                        ndB = elmntlnknd(elmntA,3)
                    END IF

                    IF (i == elmntlnknd(elmntA,2)) THEN
                        ndC = elmntlnknd(elmntA,2)
                        ndA = elmntlnknd(elmntA,3)
                        ndB = elmntlnknd(elmntA,1)
                    END IF

                    IF (i == elmntlnknd(elmntA,3)) THEN
                        ndC = elmntlnknd(elmntA,3)
                        ndA = elmntlnknd(elmntA,1)
                        ndB = elmntlnknd(elmntA,2)
                    END IF

                    tpvctACx = xnd(ndA) - xnd(ndC)
                    tpvctACy = ynd(ndA) - ynd(ndC)
                    tpvctACz = znd(ndA) - znd(ndC)

                    tpvctBCx = xnd(ndB) - xnd(ndC)
                    tpvctBCy = ynd(ndB) - ynd(ndC)
                    tpvctBCz = znd(ndB) - znd(ndC)

                    tpelmntnrmlx = tpvctACy*tpvctBCz - tpvctACz*tpvctBCy
                    tpelmntnrmly = tpvctACz*tpvctBCx - tpvctACx*tpvctBCz
                    tpelmntnrmlz = tpvctACx*tpvctBCy - tpvctACy*tpvctBCx

                    tpelmntareanrml=DSQRT(tpelmntnrmlx**2+tpelmntnrmly**2+tpelmntnrmlz**2)

                    tpnx = tpelmntnrmlx / tpelmntareanrml
                    tpny = tpelmntnrmly / tpelmntareanrml
                    tpnz = tpelmntnrmlz / tpelmntareanrml

                    tpelmntnrmlx = (xnd(ndA) + xnd(ndB) + xnd(ndC))/3.0d0
                    tpelmntnrmly = (ynd(ndA) + ynd(ndB) + ynd(ndC))/3.0d0
                    tpelmntnrmlz = (znd(ndA) + znd(ndB) + znd(ndC))/3.0d0

                    tparea = (tpelmntnrmlx-xnd(ndC))**2 &
                    &       +(tpelmntnrmly-ynd(ndC))**2 &
                    &       +(tpelmntnrmlz-znd(ndC))**2

                    tpsumndnnx = tpsumndnnx + tpnx*tpelmntareanrml/tparea
                    tpsumndnny = tpsumndnny + tpny*tpelmntareanrml/tparea
                    tpsumndnnz = tpsumndnnz + tpnz*tpelmntareanrml/tparea

                END IF

            END DO

            tpsumndnnmdl = DSQRT(tpsumndnnx**2 + tpsumndnny**2 + tpsumndnnz**2)
            nnx(i) = tpsumndnnx/tpsumndnnmdl
            nny(i) = tpsumndnny/tpsumndnnmdl
            nnz(i) = tpsumndnnz/tpsumndnnmdl

            ! Choose t1 perpendicular to n.
            IF (DABS(nnx(i)) > 0.2d0) THEN
                tpll1 = nny(i)
                tpmm1 =-nnx(i)
                tpnn1 = 0.0d0
            ELSE
                IF (DABS(nny(i)) > 0.2d0) THEN
                    tpll1 = 0.0d0
                    tpmm1 =-nnz(i)
                    tpnn1 = nny(i)
                ELSE
                    tpll1 = nnz(i)
                    tpmm1 = 0.0d0
                    tpnn1 =-nnx(i)
                END IF
            END IF
            tplmn1 = 1.0d0/DSQRT(tpll1**2 + tpmm1**2 + tpnn1**2)
            t1x(i) = tpll1*tplmn1
            t1y(i) = tpmm1*tplmn1
            t1z(i) = tpnn1*tplmn1

            ! Complete the right-handed frame with t2 = n cross t1.
            tpll2 = nny(i)*t1z(i)&
            &      -nnz(i)*t1y(i)
            tpmm2 = nnz(i)*t1x(i)&
            &      -nnx(i)*t1z(i)
            tpnn2 = nnx(i)*t1y(i)&
            &      -nny(i)*t1x(i)
            tplmn2 = 1.0d0/DSQRT(tpll2**2 + tpmm2**2 + tpnn2**2)
            t2x(i) = tpll2*tplmn2
            t2y(i) = tpmm2*tplmn2
            t2z(i) = tpnn2*tplmn2

        END DO
!$OMP END DO
!$OMP END PARALLEL

!==============

        IF (dt_option == 1) THEN
!==============
!   Calculate d/dt coefficients
!
!   Weights are chosen as
!   1 / (|g_j - v|) where g_j is the center of the triangle

!$OMP PARALLEL PRIVATE(i,j,k,elmntA,kk,slfNdA,slfNdB,slfNdC) &
!$OMP & PRIVATE(tpendx1,tpendx2,tpendx3) &
!$OMP & PRIVATE(tpendy1,tpendy2,tpendy3) &
!$OMP & PRIVATE(tpendz1,tpendz2,tpendz3) &
!$OMP & PRIVATE(tpvctACx,tpvctACy,tpvctACz) &
!$OMP & PRIVATE(tpvctBCx,tpvctBCy,tpvctBCz) &
!$OMP & PRIVATE(tpelmntnrmlx,tpelmntnrmly,tpelmntnrmlz) &
!$OMP & PRIVATE(tpelmntareanrml,tparea) &
!$OMP & PRIVATE(tp,tp1,tp2,tp3,tp4,tp5,tp6,tpxi,tpet)&
!$OMP & PRIVATE(drx_deps,drx_dyet,dry_deps,dry_dyet,drz_deps,drz_dyet) &
!$OMP & PRIVATE(h_eps,h_yet,tpkappa,t_teps,t_tyet) &
!$OMP & PRIVATE(dphi_deps1,dphi_deps2,dphi_deps3) &
!$OMP & PRIVATE(dphi_dyet1,dphi_dyet2,dphi_dyet3)
!$OMP DO
            DO i = 1, ttlnmbrnd

                DO j = 1, ttlddtnd
                    d_dt1(i,j) = 0.0d0
                    d_dt2(i,j) = 0.0d0
                END DO
                tparea = 0.0d0

                DO k = 1, mxnmbrndlnkelmnt

                    IF (ndlnkelmnt(i,k) /= 0) THEN

                        elmntA = ndlnkelmnt(i,k)

                        DO kk = 1, mxnmbrndlnknd1stslf
                            IF (elmntlnknd(elmntA,1) == ndlnknd1stslf(i,kk)) THEN
                                slfNdA = kk
                            END IF
                            IF (elmntlnknd(elmntA,2) == ndlnknd1stslf(i,kk)) THEN
                                slfNdB = kk
                            END IF
                            IF (elmntlnknd(elmntA,3) == ndlnknd1stslf(i,kk)) THEN
                                slfNdC = kk
                            END IF
                        END DO

                        tpendx1 = xnd(elmntlnknd(elmntA,1))
                        tpendx2 = xnd(elmntlnknd(elmntA,2))
                        tpendx3 = xnd(elmntlnknd(elmntA,3))
                        tpendy1 = ynd(elmntlnknd(elmntA,1))
                        tpendy2 = ynd(elmntlnknd(elmntA,2))
                        tpendy3 = ynd(elmntlnknd(elmntA,3))
                        tpendz1 = znd(elmntlnknd(elmntA,1))
                        tpendz2 = znd(elmntlnknd(elmntA,2))
                        tpendz3 = znd(elmntlnknd(elmntA,3))

                        tpvctACx = tpendx1 - tpendx3
                        tpvctACy = tpendy1 - tpendy3
                        tpvctACz = tpendz1 - tpendz3
                        tpvctBCx = tpendx2 - tpendx3
                        tpvctBCy = tpendy2 - tpendy3
                        tpvctBCz = tpendz2 - tpendz3
                        tpelmntnrmlx = tpvctACy*tpvctBCz - tpvctACz*tpvctBCy
                        tpelmntnrmly = tpvctACz*tpvctBCx - tpvctACx*tpvctBCz
                        tpelmntnrmlz = tpvctACx*tpvctBCy - tpvctACy*tpvctBCx
                        tpelmntareanrml=DSQRT( tpelmntnrmlx**2 &
                        &                     +tpelmntnrmly**2 &
                        &                     +tpelmntnrmlz**2  )
                        tpelmntnrmlx = tpelmntnrmlx/tpelmntareanrml
                        tpelmntnrmly = tpelmntnrmly/tpelmntareanrml
                        tpelmntnrmlz = tpelmntnrmlz/tpelmntareanrml
                        tp = nnx(i)*tpelmntnrmlx &
                        &   +nny(i)*tpelmntnrmly &
                        &   +nnz(i)*tpelmntnrmlz
                        tp = MAX(0.0d0, tp)

                        tpelmntnrmlx = (tpendx1 + tpendx2 + tpendx3)/3.0d0
                        tpelmntnrmly = (tpendy1 + tpendy2 + tpendy3)/3.0d0
                        tpelmntnrmlz = (tpendz1 + tpendz2 + tpendz3)/3.0d0
                        tpelmntareanrml = (tpelmntnrmlx-xnd(i))**2 &
                        &                +(tpelmntnrmly-ynd(i))**2 &
                        &                +(tpelmntnrmlz-znd(i))**2
                        tpelmntareanrml = tp / DSQRT(tpelmntareanrml)


                        IF (i == elmntlnknd(elmntA,1)) THEN
                            tpxi = 0.0d0
                            tpet = 0.0d0
                        END IF
                        IF (i == elmntlnknd(elmntA,2)) THEN
                            tpxi = 1.0d0
                            tpet = 0.0d0
                        END IF
                        IF (i == elmntlnknd(elmntA,3)) THEN
                            tpxi = 0.0d0
                            tpet = 1.0d0
                        END IF

                        tp2= 1.0d0
                        tp3= 0.0d0
                        tp1=-tp2-tp3

                        drx_deps = tp1*tpendx1+tp2*tpendx2+tp3*tpendx3
                        dry_deps = tp1*tpendy1+tp2*tpendy2+tp3*tpendy3
                        drz_deps = tp1*tpendz1+tp2*tpendz2+tp3*tpendz3

                        h_eps = DSQRT(drx_deps**2+dry_deps**2+drz_deps**2)
                        drx_deps = drx_deps/h_eps
                        dry_deps = dry_deps/h_eps
                        drz_deps = drz_deps/h_eps

                        dphi_deps1 = tp1/h_eps
                        dphi_deps2 = tp2/h_eps
                        dphi_deps3 = tp3/h_eps

                        tp2= 0.0d0
                        tp3= 1.0d0
                        tp1=-tp2-tp3

                        drx_dyet = tp1*tpendx1+tp2*tpendx2+tp3*tpendx3
                        dry_dyet = tp1*tpendy1+tp2*tpendy2+tp3*tpendy3
                        drz_dyet = tp1*tpendz1+tp2*tpendz2+tp3*tpendz3

                        h_yet = DSQRT(drx_dyet**2+dry_dyet**2+drz_dyet**2)
                        drx_dyet = drx_dyet/h_yet
                        dry_dyet = dry_dyet/h_yet
                        drz_dyet = drz_dyet/h_yet

                        dphi_dyet1 = tp1/h_yet
                        dphi_dyet2 = tp2/h_yet
                        dphi_dyet3 = tp3/h_yet

                        tpkappa = drx_deps*drx_dyet+dry_deps*dry_dyet+drz_deps*drz_dyet

                        t_teps = t1x(i)*drx_deps&
                        &       +t1y(i)*dry_deps&
                        &       +t1z(i)*drz_deps

                        t_tyet = t1x(i)*drx_dyet&
                        &       +t1y(i)*dry_dyet&
                        &       +t1z(i)*drz_dyet

                        tp1 = 1.0d0/(1.0d0-tpkappa**2)  *t_teps
                        tp2 = tpkappa/(1.0d0-tpkappa**2)*t_tyet
                        tp3 = 1.0d0/(1.0d0-tpkappa**2)  *t_tyet
                        tp4 = tpkappa/(1.0d0-tpkappa**2)*t_teps

                        tp = (tp1 - tp2)*dphi_deps1 + (tp3 - tp4)*dphi_dyet1
                        d_dt1(i,slfndA) = d_dt1(i,slfndA) + tp*tpelmntareanrml
                        tp = (tp1 - tp2)*dphi_deps2 + (tp3 - tp4)*dphi_dyet2
                        d_dt1(i,slfndB) = d_dt1(i,slfndB) + tp*tpelmntareanrml
                        tp = (tp1 - tp2)*dphi_deps3 + (tp3 - tp4)*dphi_dyet3
                        d_dt1(i,slfndC) = d_dt1(i,slfndC) + tp*tpelmntareanrml

                        t_teps = t2x(i)*drx_deps&
                        &       +t2y(i)*dry_deps&
                        &       +t2z(i)*drz_deps

                        t_tyet = t2x(i)*drx_dyet&
                        &       +t2y(i)*dry_dyet&
                        &       +t2z(i)*drz_dyet

                        tp1 = 1.0d0/(1.0d0-tpkappa**2)  *t_teps
                        tp2 = tpkappa/(1.0d0-tpkappa**2)*t_tyet
                        tp3 = 1.0d0/(1.0d0-tpkappa**2)  *t_tyet
                        tp4 = tpkappa/(1.0d0-tpkappa**2)*t_teps

                        tp = (tp1 - tp2)*dphi_deps1 + (tp3 - tp4)*dphi_dyet1
                        d_dt2(i,slfndA) = d_dt2(i,slfndA) + tp*tpelmntareanrml
                        tp = (tp1 - tp2)*dphi_deps2 + (tp3 - tp4)*dphi_dyet2
                        d_dt2(i,slfndB) = d_dt2(i,slfndB) + tp*tpelmntareanrml
                        tp = (tp1 - tp2)*dphi_deps3 + (tp3 - tp4)*dphi_dyet3
                        d_dt2(i,slfndC) = d_dt2(i,slfndC) + tp*tpelmntareanrml

                        tparea = tparea + tpelmntareanrml
                    END IF

                END DO

                DO j = 1, ttlddtnd
                    d_dt1(i,j) = d_dt1(i,j)/tparea
                    d_dt2(i,j) = d_dt2(i,j)/tparea
                END DO

            END DO
!$OMP END DO
!$OMP END PARALLEL
!==============

        END IF


        IF (dt_option == 2) THEN

!==============
!   Calculate d/dt coefficients
!
!   Jiao&Zha_Consistent Computation of First- and Second-Order
!   Differential Quantities for Surface Meshes
!
            DO i = 1, ttlnmbrnd

                ! Count nodes in the two-ring stencil.
                itmp = 0
                DO j = 1, mxnmbrndlnknd2ndslf
                    IF (ndlnknd2ndslf(i,j)/=0) THEN

                        itmp = itmp + 1

                    END IF
                END DO

                ALLOCATE (cbcmtrxA(itmp,6))
                ALLOCATE (cbcmtrxA_tr(6,itmp))
                ALLOCATE (cbcmtrxEnd(6,itmp))
                ALLOCATE (cbcmtrxAA(6,6))
                ALLOCATE (cbcmtrxAArcd(6,6))
                ALLOCATE (cbcmtrxAA_inv(6,6))
                ALLOCATE (cbcmtrxB(itmp))
                ALLOCATE (cbcmtrxBB(6))
                ALLOCATE (cbcmtrxXX(6))
                ALLOCATE (cbcmtrxW(itmp))
                ALLOCATE (cbcmtrxS(6))

                DO ij = 1, itmp
                    cbcmtrxW(ij) = 0.0d0    !weight least square coefficient matrix W
                END DO
                DO ij = 1, 6
                    cbcmtrxS(ij) = 0.0d0    !stable matrix S
                END DO

                itmp = 0
                tp = 0.0d0
                DO j = 1, mxnmbrndlnknd2ndslf
                    IF (ndlnknd2ndslf(i,j)/=0) THEN
                        itmp = itmp + 1

                        k = ndlnknd2ndslf(i,j)

                        tpcbcx = xnd(k) - xnd(i)
                        tpcbcy = ynd(k) - ynd(i)
                        tpcbcz = znd(k) - znd(i)

                        tp2= tpcbcx**2+tpcbcy**2+tpcbcz**2      !|u|^2

                        tp = tp + tp2

                    END IF
                END DO

                tp = tp/(100.0d0*DBLE(itmp))                    !eps

                DO ik = 1, itmp
                    cbcmtrxW(ik) = tp                           !eps
                END DO

                itmp = 0
                DO j = 1, mxnmbrndlnknd2ndslf

                    IF (ndlnknd2ndslf(i,j)/=0) THEN
                        itmp = itmp + 1

                        k = ndlnknd2ndslf(i,j)

                        tpcbcx = xnd(k) - xnd(i)
                        tpcbcy = ynd(k) - ynd(i)
                        tpcbcz = znd(k) - znd(i)

                        lcbcx =  tpcbcx*t1x(i)&
                        &       +tpcbcy*t1y(i)&
                        &       +tpcbcz*t1z(i)
                        lcbcy =  tpcbcx*t2x(i)&
                        &       +tpcbcy*t2y(i)&
                        &       +tpcbcz*t2z(i)
                        lcbcz =  tpcbcx*nnx(i)&
                        &       +tpcbcy*nny(i)&
                        &       +tpcbcz*nnz(i)

                        cbcmtrxA(itmp,1) = lcbcx*lcbcx
                        cbcmtrxA(itmp,2) = lcbcx*lcbcy
                        cbcmtrxA(itmp,3) = lcbcy*lcbcy
                        cbcmtrxA(itmp,4) = lcbcx
                        cbcmtrxA(itmp,5) = lcbcy
                        cbcmtrxA(itmp,6) = 1.0d0
                        cbcmtrxB(itmp) = lcbcz

                        tp1= nnx(i)*nnx(k) &
                        &   +nny(i)*nny(k) &
                        &   +nnz(i)*nnz(k)            !n0 * nk

                        tp1= MAX(0.0d0, tp1)

                        tp2= tpcbcx**2+tpcbcy**2+tpcbcz**2

                        cbcmtrxW(itmp)= tp1/(DSQRT(tp2 + cbcmtrxW(itmp)))   !entries of W

                    END IF

                END DO

                ! Weighted least-squares system: A = W V.
                DO ii = 1, itmp
                    DO ij = 1, 6
                        cbcmtrxA(ii,ij) = cbcmtrxW(ii) * cbcmtrxA(ii,ij)
                    END DO
                END DO

                ! Compute diagonal column scaling for the weighted fit.
                DO ii = 1, 6
                    cbcmtrxS(ii) = 0.0d0
                    DO ik = 1, itmp
                        cbcmtrxS(ii) = cbcmtrxS(ii) + cbcmtrxA(ik,ii)**2
                    END DO
                    cbcmtrxS(ii) = 1.0d0/cbcmtrxS(ii)
                END DO

                ! Apply the column scaling.
                DO ii = 1, 6
                    DO ij = 1, itmp
                        cbcmtrxA(ij,ii) = cbcmtrxA(ij,ii) * cbcmtrxS(ii)
                    END DO
                END DO

                !Tr(A)
                DO ii = 1, 6
                    DO ik = 1, itmp
                        cbcmtrxA_tr(ii,ik) = cbcmtrxA(ik,ii)
                    END DO
                END DO

                !Least square
                ![Tr(A) * A]^(-1)
                DO ii = 1, 6
                    DO ij = 1, 6
                        tp = 0.0d0
                        DO ik = 1, itmp
                            tp = tp + cbcmtrxA_tr(ii,ik)*cbcmtrxA(ik,ij)
                        END DO
                        cbcmtrxAArcd(ii,ij) = tp
                    END DO
                END DO

                DO ii = 1, 6
                    DO ij = 1, 6
                        DO ik = 1, 6
                            cbcmtrxAA(ij,ik) = cbcmtrxAArcd(ij,ik)
                        END DO
                        cbcmtrxBB(ij) = 0.0d0
                    END DO
                    cbcmtrxBB(ii) = 1.0d0
                    CALL dgesv (6, 1, cbcmtrxAA, 6, slv_ipiv, &
                    &           cbcmtrxBB, 6, slv_info)
                    DO ij = 1, 6
                        cbcmtrxAA_inv(ij,ii) = cbcmtrxBB(ij)
                    END DO
                END DO

                !S * [Tr(A) * A]^(-1)
                DO ii = 1, 6
                    DO ij = 1, 6
                        cbcmtrxAA_inv(ii,ij) = cbcmtrxS(ii) * cbcmtrxAA_inv(ii,ij)
                    END DO
                END DO

                !(S * [Tr(A) * A]^(-1)) * Tr(A)
                DO ii = 1, 6
                    DO ij = 1, itmp
                        tp = 0.0d0
                        DO ik = 1, 6
                            tp = tp + cbcmtrxAA_inv(ii,ik)*cbcmtrxA_tr(ik,ij)
                        END DO
                        cbcmtrxEnd(ii,ij) = tp
                    END DO
                END DO

                ![(S * [Tr(A) * A]^(-1)) * Tr(A)] * W
                DO ii = 1, itmp
                    DO ij = 1, 6
                        cbcmtrxEnd(ij,ii) = cbcmtrxEnd(ij,ii) * cbcmtrxW(ii)
                    END DO
                END DO

                DO ii = 1, ttlddtnd
                    d_dt1(i,ii) = 0.0d0
                    d_dt2(i,ii) = 0.0d0
                    d2_dt1(i,ii) = 0.0d0
                    d2_dt2(i,ii) = 0.0d0
                END DO

                DO ii = 1, itmp
                    d_dt1(i,ii) = cbcmtrxEnd(4,ii)
                    d_dt2(i,ii) = cbcmtrxEnd(5,ii)
                    d2_dt1(i,ii) = 2.0d0*cbcmtrxEnd(1,ii)
                    d2_dt2(i,ii) = 2.0d0*cbcmtrxEnd(3,ii)
                END DO

                DEALLOCATE(cbcmtrxA,cbcmtrxAA,cbcmtrxAArcd,&
                &          cbcmtrxAA_inv,cbcmtrxEnd,cbcmtrxA_tr,cbcmtrxB,&
                &          cbcmtrxW,cbcmtrxBB,cbcmtrxXX,cbcmtrxS)

            END DO

        END IF

        ! Local-frame directional curvatures:
        ! kappa_t1 = -(dn/dt1 dot t1), kappa_t2 = -(dn/dt2 dot t2).
        ! These are not generally principal curvatures; curvmn is their signed sum.

        DO i = 1, ttlnmbrnd

            ithprtl = 1
            ndst = 1
            nded = nmbrnd(1)
            elst = 1
            eled = nmbrelmnt(1)
            IF (nmbrprtl > 1 .AND. i > nded) THEN
                DO id_tp = 2, nmbrprtl
                    ndst = ndst + nmbrnd(id_tp-1)
                    nded = nded + nmbrnd(id_tp)
                    elst = elst + nmbrelmnt(id_tp-1)
                    eled = eled + nmbrelmnt(id_tp)
                    IF (i >= ndst .AND. i <= nded) THEN
                        ithprtl = id_tp
                        EXIT
                    END IF
                END DO
            END IF

            tp1 = 0.0d0
            tp2 = 0.0d0
            tp3 = 0.0d0
            tp4 = 0.0d0
            tp5 = 0.0d0
            tp6 = 0.0d0
            DO j = 1, mxnmbrndlnknd2ndslf
                IF (ndlnknd2ndslf(i,j) /= 0) THEN

                    NdA = ndlnknd2ndslf(i,j)

                    tp1 = tp1 + nnx(NdA) * d_dt1(i,j)
                    tp2 = tp2 + nny(NdA) * d_dt1(i,j)
                    tp3 = tp3 + nnz(NdA) * d_dt1(i,j)

                    tp4 = tp4 + nnx(NdA) * d_dt2(i,j)
                    tp5 = tp5 + nny(NdA) * d_dt2(i,j)
                    tp6 = tp6 + nnz(NdA) * d_dt2(i,j)

                END IF
            END DO

            curvt1(i) =-t1x(i)*tp1-t1y(i)*tp2-t1z(i)*tp3

            curvt2(i) =-t2x(i)*tp4-t2y(i)*tp5-t2z(i)*tp6

            curvmn(i) = curvt1(i) + curvt2(i)

            CYCLE

            IF (PrtlType(ithprtl) == 'Sphr' .OR. PrtlType(ithprtl) == 'PrSp' .OR. &
            &   PrtlType(ithprtl) == 'ObSp' ) THEN

                CALL dF_curv('z',ithprtl,i,tp2,tp4,tp6)

                CALL dF_curv('x',ithprtl,i,tp1,tp3,tp5)
                curvt1th(i) =-tp2**2*tp3+2.0d0*tp1*tp2*tp5-tp1**2*tp4
                curvt1th(i) = curvt1th(i)/(tp1**2+tp2**2)**(1.5d0)
                curvt1th(i) = curvt1th(i)/sizezoom(ithprtl)

                CALL dF_curv('y',ithprtl,i,tp1,tp3,tp5)
                curvt2th(i) =-tp2**2*tp3+2.0d0*tp1*tp2*tp5-tp1**2*tp4
                curvt2th(i) = curvt2th(i)/(tp1**2+tp2**2)**(1.5d0)
                curvt2th(i) = curvt2th(i)/sizezoom(ithprtl)

                curvmnth(i) = curvt1th(i) + curvt2th(i)

                IF (NrmlInOut(ithprtl) ==-1) THEN
                    curvt1th(i) =-curvt1th(i)
                    curvt2th(i) =-curvt2th(i)
                    curvmnth(i) =-curvmnth(i)
                END IF

            END IF

        END DO


!==============

!==============

        DO ithprtl = 1, nmbrprtl
            tp = 0.0d0
!$OMP PARALLEL PRIVATE(tp2)
            tp2 = 0.0d0
!$OMP DO SCHEDULE(GUIDED,4) PRIVATE(k,GLQi,icnt)
            DO k = elstaID(ithprtl), elendID(ithprtl)
                DO GLQi = 1, n_glqtr2d
                    icnt = n_glqtr2d*(k-1)+GLQi
                    tp2 = tp2 + srcfmm_wght(icnt)
                END DO
            END DO
!$OMP END DO
!$OMP ATOMIC
            tp = tp + tp2
!$OMP END PARALLEL
            surfarea(ithprtl) = tp
        END DO


        DO ithprtl = 1, nmbrprtl
            tp = 0.0d0
!$OMP PARALLEL PRIVATE(tp2)
            tp2 = 0.0d0
!$OMP DO SCHEDULE(GUIDED,4) PRIVATE(tpAx,tpAy,tpAz,tpBx,tpBy,tpBz,tpCx,tpCy,tpCz,tp1)
            DO k = elstaID(ithprtl), elendID(ithprtl)
                tpAx = xnd(elmntlnknd(k,1))
                tpAy = ynd(elmntlnknd(k,1))
                tpAz = znd(elmntlnknd(k,1))
                tpBx = xnd(elmntlnknd(k,2))
                tpBy = ynd(elmntlnknd(k,2))
                tpBz = znd(elmntlnknd(k,2))
                tpCx = xnd(elmntlnknd(k,3))
                tpCy = ynd(elmntlnknd(k,3))
                tpCz = znd(elmntlnknd(k,3))
                tp1 = - tpAx*tpBy*tpCz + tpAx*tpCy*tpBz + tpBx*tpAy*tpCz &
                    & - tpBx*tpCy*tpAz - tpCx*tpAy*tpBz + tpCx*tpBy*tpAz
                tp2 = tp2 + tp1
            END DO
!$OMP END DO
!$OMP ATOMIC
            tp = tp + tp2
!$OMP END PARALLEL
            IF (NrmlInOut(ithprtl) == 1) volume(ithprtl) = tp/6.0d0
            IF (NrmlInOut(ithprtl) ==-1) volume(ithprtl) =-tp/6.0d0
        END DO

    END SUBROUTINE



    ! Compute nodal normals, tangent frames, curvatures, and surface-derivative
    ! weights for the Q6 mesh used by all released manuscript cases.

    SUBROUTINE GetndnrmlQdrtcLnr

        INTEGER :: ithprtl, i, j, k, ii, ij, ik, jj, kk, itmp, elmntA, ndA, ndB, ndC
        INTEGER :: ndst,nded,elst,eled,id_tp
        INTEGER, DIMENSION (mxnmbrndlnknd2nd,2) :: Nd2ndPairRcd

        DOUBLE PRECISION :: tpvctACx, tpvctACy, tpvctACz, &
                         &  tpvctBCx, tpvctBCy, tpvctBCz, &
                         &  tpelmntnrmlx, tpelmntnrmly, tpelmntnrmlz, &
                         &  tpelmntareanrml

        DOUBLE PRECISION :: tpsumndnnx, tpsumndnny, tpsumndnnz, tpsumndnnmdl
        DOUBLE PRECISION :: tpnx, tpny, tpnz, tparea

        DOUBLE PRECISION ::  dmnx, dmny, dmnz, tplmn1, tplmn2, tpll1, tpll2, tpll3, &
                            &tpmm1, tpmm2, tpmm3, tpnn1, tpnn2, tpnn3

        DOUBLE PRECISION :: dtetax,dtetay,dtetaz,dleta,dt_xix,dt_xiy,dt_xiz,dt_xi,&
                           &dteta1,dteta2,dt_xi1,dt_xi2,dfac,dl_xi
        INTEGER :: slfNdA,slfNdB,slfNdC,slfNdD,slfNdE,slfNdF

        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: &
        &   cbcmtrxA,cbcmtrxAA,cbcmtrxAArcd,cbcmtrxAA_inv,cbcmtrxEnd,cbcmtrxA_tr
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: &
        &   cbcmtrxB,cbcmtrxW,cbcmtrxBB,cbcmtrxXX,cbcmtrxS
        INTEGER :: slv_ipiv(6),slv_D,slv_info
        DOUBLE PRECISION :: tp, tpcbcx, tpcbcy, tpcbcz, lcbcx, lcbcy, lcbcz

        INTEGER ::  GLQi, icnt
        DOUBLE PRECISION :: tp1,tp2,tp3,tp4,tp5,tp6,tp_a,tp_b,tp_g,tpxi,tpet
        DOUBLE PRECISION :: tprx,tpry,tprz,tpEnx,tpEny,tpEnz,JcbDtmn
        DOUBLE PRECISION :: tpendx1,tpendx2,tpendx3,tpendx4,tpendx5,tpendx6
        DOUBLE PRECISION :: tpendy1,tpendy2,tpendy3,tpendy4,tpendy5,tpendy6
        DOUBLE PRECISION :: tpendz1,tpendz2,tpendz3,tpendz4,tpendz5,tpendz6
        DOUBLE PRECISION :: drx_deps,drx_dyet,dry_deps,dry_dyet,drz_deps,drz_dyet
        DOUBLE PRECISION :: dphi_deps1,dphi_deps2,dphi_deps3,&
        &                   dphi_deps4,dphi_deps5,dphi_deps6
        DOUBLE PRECISION :: dphi_dyet1,dphi_dyet2,dphi_dyet3,&
        &                   dphi_dyet4,dphi_dyet5,dphi_dyet6
        DOUBLE PRECISION :: h_eps,h_yet,tpkappa,t_teps,t_tyet

        DOUBLE PRECISION :: tpAx, tpAy, tpAz, tpBx, tpBy, tpBz, tpCx, tpCy, tpCz

        CHARACTER (LEN=1) :: ndnrml_mth = 'L'

!==============
!   Gauss point weight, nodal position, nodal normal vector on each element

!$OMP PARALLEL PRIVATE(k,GLQi,icnt) &
!$OMP & PRIVATE(tp1,tp2,tp3,tp4,tp5,tp6,tp_a,tp_b,tp_g,tpxi,tpet)&
!$OMP & PRIVATE(tprx,tpry,tprz,tpEnx,tpEny,tpEnz,JcbDtmn) &
!$OMP & PRIVATE(tpendx1,tpendx2,tpendx3,tpendx4,tpendx5,tpendx6) &
!$OMP & PRIVATE(tpendy1,tpendy2,tpendy3,tpendy4,tpendy5,tpendy6) &
!$OMP & PRIVATE(tpendz1,tpendz2,tpendz3,tpendz4,tpendz5,tpendz6) &
!$OMP & PRIVATE(drx_deps,drx_dyet,dry_deps,dry_dyet,drz_deps,drz_dyet)
!$OMP DO
        DO k = 1, ttlnmbrelmnt

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

            DO GLQi = 1, n_glqtr2d

                icnt = n_glqtr2d*(k-1)+GLQi

                tpxi = xg_glqtr2d(GLQi)
                tpet = yg_glqtr2d(GLQi)

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

                tprx = tp1*tpendx1+tp2*tpendx2+tp3*tpendx3 &
                &     +tp4*tpendx4+tp5*tpendx5+tp6*tpendx6
                tpry = tp1*tpendy1+tp2*tpendy2+tp3*tpendy3 &
                &     +tp4*tpendy4+tp5*tpendy5+tp6*tpendy6
                tprz = tp1*tpendz1+tp2*tpendz2+tp3*tpendz3 &
                &     +tp4*tpendz4+tp5*tpendz5+tp6*tpendz6

                srcfmm_vec(1, icnt) = tprx
                srcfmm_vec(2, icnt) = tpry
                srcfmm_vec(3, icnt) = tprz

                srcfmm_nrm(1, icnt) = tpEnx
                srcfmm_nrm(2, icnt) = tpEny
                srcfmm_nrm(3, icnt) = tpEnz

                srcfmm_wght(icnt) = wg_glqtr2d(GLQi)*JcbDtmn

                srcfmm_wtnd(1,icnt) = wg_glqtr2d(GLQi)*JcbDtmn*tp1
                srcfmm_wtnd(2,icnt) = wg_glqtr2d(GLQi)*JcbDtmn*tp2
                srcfmm_wtnd(3,icnt) = wg_glqtr2d(GLQi)*JcbDtmn*tp3
                srcfmm_wtnd(4,icnt) = wg_glqtr2d(GLQi)*JcbDtmn*tp4
                srcfmm_wtnd(5,icnt) = wg_glqtr2d(GLQi)*JcbDtmn*tp5
                srcfmm_wtnd(6,icnt) = wg_glqtr2d(GLQi)*JcbDtmn*tp6

            END DO

        END DO
!$OMP END DO
!$OMP END PARALLEL
!==============

!==============
!   Normal, tangential vector calculation
!   Normal, tangential vector calculation
!   Normal is based on the facet elements connecting to the node of interest
!   Weights are chosen as the combination of the element area
!   and how far the element is relative to the node
!   combine "MAX, N. Journal of Graphics Tools, Vol. 4, No. 2." &
!   "Chen, Wu Computer Aided Geometric Design 21 (2004) 447–458"
!   as: Area_j / (|g_j - v|^2) where g_j is the center of the triangle

        IF (ndnrml_mth == 'L') THEN

!$OMP PARALLEL PRIVATE(i,k,elmntA,ndA,ndB,ndC) &
!$OMP & PRIVATE(tpsumndnnx,tpsumndnny,tpsumndnnz)&
!$OMP & PRIVATE(tpvctACx,tpvctACy,tpvctACz) &
!$OMP & PRIVATE(tpvctBCx,tpvctBCy,tpvctBCz) &
!$OMP & PRIVATE(tpelmntnrmlx,tpelmntnrmly,tpelmntnrmlz) &
!$OMP & PRIVATE(tpelmntareanrml) &
!$OMP & PRIVATE(tpnx,tpny,tpnz,tparea,tpsumndnnmdl) &
!$OMP & PRIVATE(tpll1,tpmm1,tpnn1,tplmn1,tpll2,tpmm2,tpnn2,tplmn2)
!$OMP DO
            DO i = 1, ttlnmbrnd

                tpsumndnnx = 0.0d0
                tpsumndnny = 0.0d0
                tpsumndnnz = 0.0d0

                DO k = 1, mxnmbrndlnkelmntlnr

                    IF (ndlnkelmntlnr(i,k) /= 0) THEN

                        elmntA = ndlnkelmntlnr(i,k)

                        IF (i == elmntlnkndlnr(elmntA,1)) THEN
                            ndC = elmntlnkndlnr(elmntA,1)
                            ndA = elmntlnkndlnr(elmntA,2)
                            ndB = elmntlnkndlnr(elmntA,3)
                        END IF

                        IF (i == elmntlnkndlnr(elmntA,2)) THEN
                            ndC = elmntlnkndlnr(elmntA,2)
                            ndA = elmntlnkndlnr(elmntA,3)
                            ndB = elmntlnkndlnr(elmntA,1)
                        END IF

                        IF (i == elmntlnkndlnr(elmntA,3)) THEN
                            ndC = elmntlnkndlnr(elmntA,3)
                            ndA = elmntlnkndlnr(elmntA,1)
                            ndB = elmntlnkndlnr(elmntA,2)
                        END IF

                        tpvctACx = xnd(ndA) - xnd(ndC)
                        tpvctACy = ynd(ndA) - ynd(ndC)
                        tpvctACz = znd(ndA) - znd(ndC)

                        tpvctBCx = xnd(ndB) - xnd(ndC)
                        tpvctBCy = ynd(ndB) - ynd(ndC)
                        tpvctBCz = znd(ndB) - znd(ndC)

                        tpelmntnrmlx = tpvctACy*tpvctBCz - tpvctACz*tpvctBCy
                        tpelmntnrmly = tpvctACz*tpvctBCx - tpvctACx*tpvctBCz
                        tpelmntnrmlz = tpvctACx*tpvctBCy - tpvctACy*tpvctBCx

                        tpelmntareanrml=DSQRT( tpelmntnrmlx**2 &
                        &                     +tpelmntnrmly**2 &
                        &                     +tpelmntnrmlz**2  )

                        tpnx = tpelmntnrmlx / tpelmntareanrml
                        tpny = tpelmntnrmly / tpelmntareanrml
                        tpnz = tpelmntnrmlz / tpelmntareanrml

                        tpelmntnrmlx = (xnd(ndA) + xnd(ndB) + xnd(ndC))/3.0d0
                        tpelmntnrmly = (ynd(ndA) + ynd(ndB) + ynd(ndC))/3.0d0
                        tpelmntnrmlz = (znd(ndA) + znd(ndB) + znd(ndC))/3.0d0

                        tparea = (tpelmntnrmlx-xnd(ndC))**2 &
                        &       +(tpelmntnrmly-ynd(ndC))**2 &
                        &       +(tpelmntnrmlz-znd(ndC))**2

                        tpsumndnnx = tpsumndnnx + tpnx*tpelmntareanrml/tparea
                        tpsumndnny = tpsumndnny + tpny*tpelmntareanrml/tparea
                        tpsumndnnz = tpsumndnnz + tpnz*tpelmntareanrml/tparea

                    END IF

                END DO

                tpsumndnnmdl = DSQRT(tpsumndnnx**2 + tpsumndnny**2 + tpsumndnnz**2)
                nnx(i) = tpsumndnnx/tpsumndnnmdl
                nny(i) = tpsumndnny/tpsumndnnmdl
                nnz(i) = tpsumndnnz/tpsumndnnmdl

                ! Choose t1 perpendicular to n.
                IF (DABS(nnx(i)).gt.0.2d0) THEN
                    tpll1 = nny(i)
                    tpmm1 =-nnx(i)
                    tpnn1 = 0.0d0
                ELSE
                    IF (DABS(nny(i)).gt.0.2d0) THEN
                        tpll1 = 0.0d0
                        tpmm1 =-nnz(i)
                        tpnn1 = nny(i)
                    ELSE
                        tpll1 = nnz(i)
                        tpmm1 = 0.0d0
                        tpnn1 =-nnx(i)
                    END IF
                END IF
                tplmn1 = 1.0d0/DSQRT(tpll1**2 + tpmm1**2 + tpnn1**2)
                t1x(i) = tpll1*tplmn1
                t1y(i) = tpmm1*tplmn1
                t1z(i) = tpnn1*tplmn1

            ! Complete the right-handed frame with t2 = n cross t1.
            tpll2 = nny(i)*t1z(i)&
                    &  -nnz(i)*t1y(i)
                tpmm2 = nnz(i)*t1x(i)&
                    &  -nnx(i)*t1z(i)
                tpnn2 = nnx(i)*t1y(i)&
                    &  -nny(i)*t1x(i)
                tplmn2 = 1.0d0/DSQRT(tpll2**2 + tpmm2**2 + tpnn2**2)
                t2x(i) = tpll2*tplmn2
                t2y(i) = tpmm2*tplmn2
                t2z(i) = tpnn2*tplmn2

            END DO
!$OMP END DO
!$OMP END PARALLEL

            IF (dt_option == 1) THEN
!==============
!   Calculate d/dt coefficients
!
!   Weights are chosen as
!   1 / (|g_j - v|) where g_j is the center of the triangle

!$OMP PARALLEL PRIVATE(i,j,k,elmntA,kk,slfNdA,slfNdB,slfNdC) &
!$OMP & PRIVATE(tpendx1,tpendx2,tpendx3) &
!$OMP & PRIVATE(tpendy1,tpendy2,tpendy3) &
!$OMP & PRIVATE(tpendz1,tpendz2,tpendz3) &
!$OMP & PRIVATE(tpvctACx,tpvctACy,tpvctACz) &
!$OMP & PRIVATE(tpvctBCx,tpvctBCy,tpvctBCz) &
!$OMP & PRIVATE(tpelmntnrmlx,tpelmntnrmly,tpelmntnrmlz) &
!$OMP & PRIVATE(tpelmntareanrml,tparea) &
!$OMP & PRIVATE(tp,tp1,tp2,tp3,tp4,tp5,tp6,tpxi,tpet)&
!$OMP & PRIVATE(drx_deps,drx_dyet,dry_deps,dry_dyet,drz_deps,drz_dyet) &
!$OMP & PRIVATE(h_eps,h_yet,tpkappa,t_teps,t_tyet) &
!$OMP & PRIVATE(dphi_deps1,dphi_deps2,dphi_deps3) &
!$OMP & PRIVATE(dphi_dyet1,dphi_dyet2,dphi_dyet3)
!$OMP DO
            DO i = 1, ttlnmbrnd

                DO j = 1, ttlddtnd
                    d_dt1(i,j) = 0.0d0
                    d_dt2(i,j) = 0.0d0
                END DO
                tparea = 0.0d0

                DO k = 1, mxnmbrndlnkelmntlnr

                    IF (ndlnkelmntlnr(i,k) /= 0) THEN

                        elmntA = ndlnkelmntlnr(i,k)

                        DO kk = 1, mxnmbrndlnknd1stslf
                            IF (elmntlnkndlnr(elmntA,1) == ndlnknd1stslf(i,kk)) THEN
                                slfNdA = kk
                            END IF
                            IF (elmntlnkndlnr(elmntA,2) == ndlnknd1stslf(i,kk)) THEN
                                slfNdB = kk
                            END IF
                            IF (elmntlnkndlnr(elmntA,3) == ndlnknd1stslf(i,kk)) THEN
                                slfNdC = kk
                            END IF
                        END DO

                        tpendx1 = xnd(elmntlnkndlnr(elmntA,1))
                        tpendx2 = xnd(elmntlnkndlnr(elmntA,2))
                        tpendx3 = xnd(elmntlnkndlnr(elmntA,3))
                        tpendy1 = ynd(elmntlnkndlnr(elmntA,1))
                        tpendy2 = ynd(elmntlnkndlnr(elmntA,2))
                        tpendy3 = ynd(elmntlnkndlnr(elmntA,3))
                        tpendz1 = znd(elmntlnkndlnr(elmntA,1))
                        tpendz2 = znd(elmntlnkndlnr(elmntA,2))
                        tpendz3 = znd(elmntlnkndlnr(elmntA,3))

                        tpvctACx = tpendx1 - tpendx3
                        tpvctACy = tpendy1 - tpendy3
                        tpvctACz = tpendz1 - tpendz3
                        tpvctBCx = tpendx2 - tpendx3
                        tpvctBCy = tpendy2 - tpendy3
                        tpvctBCz = tpendz2 - tpendz3
                        tpelmntnrmlx = tpvctACy*tpvctBCz - tpvctACz*tpvctBCy
                        tpelmntnrmly = tpvctACz*tpvctBCx - tpvctACx*tpvctBCz
                        tpelmntnrmlz = tpvctACx*tpvctBCy - tpvctACy*tpvctBCx
                        tpelmntareanrml=DSQRT( tpelmntnrmlx**2 &
                        &                     +tpelmntnrmly**2 &
                        &                     +tpelmntnrmlz**2  )
                        tpelmntnrmlx = tpelmntnrmlx/tpelmntareanrml
                        tpelmntnrmly = tpelmntnrmly/tpelmntareanrml
                        tpelmntnrmlz = tpelmntnrmlz/tpelmntareanrml
                        tp = nnx(i)*tpelmntnrmlx &
                        &   +nny(i)*tpelmntnrmly &
                        &   +nnz(i)*tpelmntnrmlz
                        tp = MAX(0.0d0, tp)

                        tpelmntnrmlx = (tpendx1 + tpendx2 + tpendx3)/3.0d0
                        tpelmntnrmly = (tpendy1 + tpendy2 + tpendy3)/3.0d0
                        tpelmntnrmlz = (tpendz1 + tpendz2 + tpendz3)/3.0d0
                        tpelmntareanrml = (tpelmntnrmlx-xnd(i))**2 &
                        &                +(tpelmntnrmly-ynd(i))**2 &
                        &                +(tpelmntnrmlz-znd(i))**2
                        tpelmntareanrml = tp / DSQRT(tpelmntareanrml)


                        IF (i == elmntlnkndlnr(elmntA,1)) THEN
                            tpxi = 0.0d0
                            tpet = 0.0d0
                        END IF
                        IF (i == elmntlnkndlnr(elmntA,2)) THEN
                            tpxi = 1.0d0
                            tpet = 0.0d0
                        END IF
                        IF (i == elmntlnkndlnr(elmntA,3)) THEN
                            tpxi = 0.0d0
                            tpet = 1.0d0
                        END IF

                        tp2= 1.0d0
                        tp3= 0.0d0
                        tp1=-tp2-tp3

                        drx_deps = tp1*tpendx1+tp2*tpendx2+tp3*tpendx3
                        dry_deps = tp1*tpendy1+tp2*tpendy2+tp3*tpendy3
                        drz_deps = tp1*tpendz1+tp2*tpendz2+tp3*tpendz3

                        h_eps = DSQRT(drx_deps**2+dry_deps**2+drz_deps**2)
                        drx_deps = drx_deps/h_eps
                        dry_deps = dry_deps/h_eps
                        drz_deps = drz_deps/h_eps

                        dphi_deps1 = tp1/h_eps
                        dphi_deps2 = tp2/h_eps
                        dphi_deps3 = tp3/h_eps

                        tp2= 0.0d0
                        tp3= 1.0d0
                        tp1=-tp2-tp3

                        drx_dyet = tp1*tpendx1+tp2*tpendx2+tp3*tpendx3
                        dry_dyet = tp1*tpendy1+tp2*tpendy2+tp3*tpendy3
                        drz_dyet = tp1*tpendz1+tp2*tpendz2+tp3*tpendz3

                        h_yet = DSQRT(drx_dyet**2+dry_dyet**2+drz_dyet**2)
                        drx_dyet = drx_dyet/h_yet
                        dry_dyet = dry_dyet/h_yet
                        drz_dyet = drz_dyet/h_yet

                        dphi_dyet1 = tp1/h_yet
                        dphi_dyet2 = tp2/h_yet
                        dphi_dyet3 = tp3/h_yet

                        tpkappa = drx_deps*drx_dyet+dry_deps*dry_dyet+drz_deps*drz_dyet

                        t_teps = t1x(i)*drx_deps&
                        &       +t1y(i)*dry_deps&
                        &       +t1z(i)*drz_deps

                        t_tyet = t1x(i)*drx_dyet&
                        &       +t1y(i)*dry_dyet&
                        &       +t1z(i)*drz_dyet

                        tp1 = 1.0d0/(1.0d0-tpkappa**2)  *t_teps
                        tp2 = tpkappa/(1.0d0-tpkappa**2)*t_tyet
                        tp3 = 1.0d0/(1.0d0-tpkappa**2)  *t_tyet
                        tp4 = tpkappa/(1.0d0-tpkappa**2)*t_teps

                        tp = (tp1 - tp2)*dphi_deps1 + (tp3 - tp4)*dphi_dyet1
                        d_dt1(i,slfndA) = d_dt1(i,slfndA) + tp*tpelmntareanrml
                        tp = (tp1 - tp2)*dphi_deps2 + (tp3 - tp4)*dphi_dyet2
                        d_dt1(i,slfndB) = d_dt1(i,slfndB) + tp*tpelmntareanrml
                        tp = (tp1 - tp2)*dphi_deps3 + (tp3 - tp4)*dphi_dyet3
                        d_dt1(i,slfndC) = d_dt1(i,slfndC) + tp*tpelmntareanrml

                        t_teps = t2x(i)*drx_deps&
                        &       +t2y(i)*dry_deps&
                        &       +t2z(i)*drz_deps

                        t_tyet = t2x(i)*drx_dyet&
                        &       +t2y(i)*dry_dyet&
                        &       +t2z(i)*drz_dyet

                        tp1 = 1.0d0/(1.0d0-tpkappa**2)  *t_teps
                        tp2 = tpkappa/(1.0d0-tpkappa**2)*t_tyet
                        tp3 = 1.0d0/(1.0d0-tpkappa**2)  *t_tyet
                        tp4 = tpkappa/(1.0d0-tpkappa**2)*t_teps

                        tp = (tp1 - tp2)*dphi_deps1 + (tp3 - tp4)*dphi_dyet1
                        d_dt2(i,slfndA) = d_dt2(i,slfndA) + tp*tpelmntareanrml
                        tp = (tp1 - tp2)*dphi_deps2 + (tp3 - tp4)*dphi_dyet2
                        d_dt2(i,slfndB) = d_dt2(i,slfndB) + tp*tpelmntareanrml
                        tp = (tp1 - tp2)*dphi_deps3 + (tp3 - tp4)*dphi_dyet3
                        d_dt2(i,slfndC) = d_dt2(i,slfndC) + tp*tpelmntareanrml

                        tparea = tparea + tpelmntareanrml
                    END IF

                END DO

                DO j = 1, ttlddtnd
                    d_dt1(i,j) = d_dt1(i,j)/tparea
                    d_dt2(i,j) = d_dt2(i,j)/tparea
                END DO

            END DO
!$OMP END DO
!$OMP END PARALLEL

!==============

            END IF

        END IF


        IF (dt_option == 2) THEN

!==============
!   Calculate d/dt coefficients
!
!   Jiao&Zha_Consistent Computation of First- and Second-Order
!   Differential Quantities for Surface Meshes
!
            DO i = 1, ttlnmbrnd

                ! Count nodes in the two-ring stencil.
                itmp = 0
                DO j = 1, mxnmbrndlnknd2ndslf
                    IF (ndlnknd2ndslf(i,j)/=0) THEN

                        itmp = itmp + 1

                    END IF
                END DO

                ALLOCATE (cbcmtrxA(itmp,6))
                ALLOCATE (cbcmtrxA_tr(6,itmp))
                ALLOCATE (cbcmtrxEnd(6,itmp))
                ALLOCATE (cbcmtrxAA(6,6))
                ALLOCATE (cbcmtrxAArcd(6,6))
                ALLOCATE (cbcmtrxAA_inv(6,6))
                ALLOCATE (cbcmtrxB(itmp))
                ALLOCATE (cbcmtrxBB(6))
                ALLOCATE (cbcmtrxXX(6))
                ALLOCATE (cbcmtrxW(itmp))
                ALLOCATE (cbcmtrxS(6))

                DO ij = 1, itmp
                    cbcmtrxW(ij) = 0.0d0    !weight least square coefficient matrix W
                END DO
                DO ij = 1, 6
                    cbcmtrxS(ij) = 0.0d0    !stable matrix S
                END DO

                itmp = 0
                tp = 0.0d0
                DO j = 1, mxnmbrndlnknd2ndslf
                    IF (ndlnknd2ndslf(i,j)/=0) THEN
                        itmp = itmp + 1

                        k = ndlnknd2ndslf(i,j)

                        tpcbcx = xnd(k) - xnd(i)
                        tpcbcy = ynd(k) - ynd(i)
                        tpcbcz = znd(k) - znd(i)

                        tp2= tpcbcx**2+tpcbcy**2+tpcbcz**2      !|u|^2

                        tp = tp + tp2

                    END IF
                END DO

                tp = tp/(100.0d0*DBLE(itmp))                    !eps

                DO ik = 1, itmp
                    cbcmtrxW(ik) = tp                           !eps
                END DO

                itmp = 0
                DO j = 1, mxnmbrndlnknd2ndslf

                    IF (ndlnknd2ndslf(i,j)/=0) THEN
                        itmp = itmp + 1

                        k = ndlnknd2ndslf(i,j)

                        tpcbcx = xnd(k) - xnd(i)
                        tpcbcy = ynd(k) - ynd(i)
                        tpcbcz = znd(k) - znd(i)

                        lcbcx =  tpcbcx*t1x(i)&
                        &       +tpcbcy*t1y(i)&
                        &       +tpcbcz*t1z(i)
                        lcbcy =  tpcbcx*t2x(i)&
                        &       +tpcbcy*t2y(i)&
                        &       +tpcbcz*t2z(i)
                        lcbcz =  tpcbcx*nnx(i)&
                        &       +tpcbcy*nny(i)&
                        &       +tpcbcz*nnz(i)

                        cbcmtrxA(itmp,1) = lcbcx*lcbcx
                        cbcmtrxA(itmp,2) = lcbcx*lcbcy
                        cbcmtrxA(itmp,3) = lcbcy*lcbcy
                        cbcmtrxA(itmp,4) = lcbcx
                        cbcmtrxA(itmp,5) = lcbcy
                        cbcmtrxA(itmp,6) = 1.0d0
                        cbcmtrxB(itmp) = lcbcz

                        tp1= nnx(i)*nnx(k) &
                        &   +nny(i)*nny(k) &
                        &   +nnz(i)*nnz(k)            !n0 * nk

                        tp1= MAX(0.0d0, tp1)

                        tp2= tpcbcx**2+tpcbcy**2+tpcbcz**2

                        cbcmtrxW(itmp)= tp1/(DSQRT(tp2 + cbcmtrxW(itmp)))   !entries of W

                    END IF

                END DO

                ! Weighted least-squares system: A = W V.
                DO ii = 1, itmp
                    DO ij = 1, 6
                        cbcmtrxA(ii,ij) = cbcmtrxW(ii) * cbcmtrxA(ii,ij)
                    END DO
                END DO

                ! Compute diagonal column scaling for the weighted fit.
                DO ii = 1, 6
                    cbcmtrxS(ii) = 0.0d0
                    DO ik = 1, itmp
                        cbcmtrxS(ii) = cbcmtrxS(ii) + cbcmtrxA(ik,ii)**2
                    END DO
                    cbcmtrxS(ii) = 1.0d0/cbcmtrxS(ii)
                END DO

                ! Apply the column scaling.
                DO ii = 1, 6
                    DO ij = 1, itmp
                        cbcmtrxA(ij,ii) = cbcmtrxA(ij,ii) * cbcmtrxS(ii)
                    END DO
                END DO

                !Tr(A)
                DO ii = 1, 6
                    DO ik = 1, itmp
                        cbcmtrxA_tr(ii,ik) = cbcmtrxA(ik,ii)
                    END DO
                END DO

                !Least square
                ![Tr(A) * A]^(-1)
                DO ii = 1, 6
                    DO ij = 1, 6
                        tp = 0.0d0
                        DO ik = 1, itmp
                            tp = tp + cbcmtrxA_tr(ii,ik)*cbcmtrxA(ik,ij)
                        END DO
                        cbcmtrxAArcd(ii,ij) = tp
                    END DO
                END DO

                DO ii = 1, 6
                    DO ij = 1, 6
                        DO ik = 1, 6
                            cbcmtrxAA(ij,ik) = cbcmtrxAArcd(ij,ik)
                        END DO
                        cbcmtrxBB(ij) = 0.0d0
                    END DO
                    cbcmtrxBB(ii) = 1.0d0
                    CALL dgesv (6, 1, cbcmtrxAA, 6, slv_ipiv, &
                    &           cbcmtrxBB, 6, slv_info)
                    DO ij = 1, 6
                        cbcmtrxAA_inv(ij,ii) = cbcmtrxBB(ij)
                    END DO
                END DO

                !S * [Tr(A) * A]^(-1)
                DO ii = 1, 6
                    DO ij = 1, 6
                        cbcmtrxAA_inv(ii,ij) = cbcmtrxS(ii) * cbcmtrxAA_inv(ii,ij)
                    END DO
                END DO

                !(S * [Tr(A) * A]^(-1)) * Tr(A)
                DO ii = 1, 6
                    DO ij = 1, itmp
                        tp = 0.0d0
                        DO ik = 1, 6
                            tp = tp + cbcmtrxAA_inv(ii,ik)*cbcmtrxA_tr(ik,ij)
                        END DO
                        cbcmtrxEnd(ii,ij) = tp
                    END DO
                END DO

                ![(S * [Tr(A) * A]^(-1)) * Tr(A)] * W
                DO ii = 1, itmp
                    DO ij = 1, 6
                        cbcmtrxEnd(ij,ii) = cbcmtrxEnd(ij,ii) * cbcmtrxW(ii)
                    END DO
                END DO

                DO ii = 1, ttlddtnd
                    d_dt1(i,ii) = 0.0d0
                    d_dt2(i,ii) = 0.0d0
                    d2_dt1(i,ii) = 0.0d0
                    d2_dt2(i,ii) = 0.0d0
                END DO

                DO ii = 1, itmp
                    d_dt1(i,ii) = cbcmtrxEnd(4,ii)
                    d_dt2(i,ii) = cbcmtrxEnd(5,ii)
                    d2_dt1(i,ii) = 2.0d0*cbcmtrxEnd(1,ii)
                    d2_dt2(i,ii) = 2.0d0*cbcmtrxEnd(3,ii)
                END DO

                DEALLOCATE(cbcmtrxA,cbcmtrxAA,cbcmtrxAArcd,&
                &          cbcmtrxAA_inv,cbcmtrxEnd,cbcmtrxA_tr,cbcmtrxB,&
                &          cbcmtrxW,cbcmtrxBB,cbcmtrxXX,cbcmtrxS)

            END DO

        END IF

        ! Local-frame directional curvatures:
        ! kappa_t1 = -(dn/dt1 dot t1), kappa_t2 = -(dn/dt2 dot t2).
        ! These are not generally principal curvatures; curvmn is their signed sum.

        DO i = 1, ttlnmbrnd

            ithprtl = 1
            ndst = 1
            nded = nmbrnd(1)
            elst = 1
            eled = nmbrelmnt(1)
            IF (nmbrprtl > 1 .AND. i > nded) THEN
                DO id_tp = 2, nmbrprtl
                    ndst = ndst + nmbrnd(id_tp-1)
                    nded = nded + nmbrnd(id_tp)
                    elst = elst + nmbrelmnt(id_tp-1)
                    eled = eled + nmbrelmnt(id_tp)
                    IF (i >= ndst .AND. i <= nded) THEN
                        ithprtl = id_tp
                        EXIT
                    END IF
                END DO
            END IF

            tp1 = 0.0d0
            tp2 = 0.0d0
            tp3 = 0.0d0
            tp4 = 0.0d0
            tp5 = 0.0d0
            tp6 = 0.0d0
            DO j = 1, mxnmbrndlnknd2ndslf
                IF (ndlnknd2ndslf(i,j) /= 0) THEN

                    NdA = ndlnknd2ndslf(i,j)

                    tp1 = tp1 + nnx(NdA) * d_dt1(i,j)
                    tp2 = tp2 + nny(NdA) * d_dt1(i,j)
                    tp3 = tp3 + nnz(NdA) * d_dt1(i,j)

                    tp4 = tp4 + nnx(NdA) * d_dt2(i,j)
                    tp5 = tp5 + nny(NdA) * d_dt2(i,j)
                    tp6 = tp6 + nnz(NdA) * d_dt2(i,j)

                END IF
            END DO

            curvt1(i) =-t1x(i)*tp1-t1y(i)*tp2-t1z(i)*tp3

            curvt2(i) =-t2x(i)*tp4-t2y(i)*tp5-t2z(i)*tp6

            curvmn(i) = curvt1(i) + curvt2(i)

            CYCLE

            IF (PrtlType(ithprtl) == 'Sphr' .OR. PrtlType(ithprtl) == 'PrSp' .OR. &
            &   PrtlType(ithprtl) == 'ObSp' ) THEN

                CALL dF_curv('z',ithprtl,i,tp2,tp4,tp6)

                CALL dF_curv('x',ithprtl,i,tp1,tp3,tp5)
                curvt1th(i) =-tp2**2*tp3+2.0d0*tp1*tp2*tp5-tp1**2*tp4
                curvt1th(i) = curvt1th(i)/(tp1**2+tp2**2)**(1.5d0)
                curvt1th(i) = curvt1th(i)/sizezoom(ithprtl)

                CALL dF_curv('y',ithprtl,i,tp1,tp3,tp5)
                curvt2th(i) =-tp2**2*tp3+2.0d0*tp1*tp2*tp5-tp1**2*tp4
                curvt2th(i) = curvt2th(i)/(tp1**2+tp2**2)**(1.5d0)
                curvt2th(i) = curvt2th(i)/sizezoom(ithprtl)

                curvmnth(i) = curvt1th(i) + curvt2th(i)

                IF (NrmlInOut(ithprtl) ==-1) THEN
                    curvt1th(i) =-curvt1th(i)
                    curvt2th(i) =-curvt2th(i)
                    curvmnth(i) =-curvmnth(i)
                END IF

            END IF

        END DO

!==============

!==============

        DO ithprtl = 1, nmbrprtl
            tp = 0.0d0
!$OMP PARALLEL PRIVATE(tp2)
            tp2 = 0.0d0
!$OMP DO SCHEDULE(GUIDED,4) PRIVATE(k,GLQi,icnt)
            DO k = elstaID(ithprtl), elendID(ithprtl)
                DO GLQi = 1, n_glqtr2d
                    icnt = n_glqtr2d*(k-1)+GLQi
                    tp2 = tp2 + srcfmm_wght(icnt)
                END DO
            END DO
!$OMP END DO
!$OMP ATOMIC
            tp = tp + tp2
!$OMP END PARALLEL
            surfarea(ithprtl) = tp
        END DO

        DO ithprtl = 1, nmbrprtl
            tp = 0.0d0
!$OMP PARALLEL PRIVATE(tp2)
            tp2 = 0.0d0
!$OMP DO SCHEDULE(GUIDED,4) PRIVATE(tpAx,tpAy,tpAz,tpBx,tpBy,tpBz,tpCx,tpCy,tpCz,tp1)
            DO k = elstaID(ithprtl), elendID(ithprtl)
                tpAx = xnd(elmntlnknd(k,1))
                tpAy = ynd(elmntlnknd(k,1))
                tpAz = znd(elmntlnknd(k,1))
                tpBx = xnd(elmntlnknd(k,4))
                tpBy = ynd(elmntlnknd(k,4))
                tpBz = znd(elmntlnknd(k,4))
                tpCx = xnd(elmntlnknd(k,6))
                tpCy = ynd(elmntlnknd(k,6))
                tpCz = znd(elmntlnknd(k,6))
                tp1 = - tpAx*tpBy*tpCz + tpAx*tpCy*tpBz + tpBx*tpAy*tpCz &
                    & - tpBx*tpCy*tpAz - tpCx*tpAy*tpBz + tpCx*tpBy*tpAz
                tp2 = tp2 + tp1
                tpAx = xnd(elmntlnknd(k,4))
                tpAy = ynd(elmntlnknd(k,4))
                tpAz = znd(elmntlnknd(k,4))
                tpBx = xnd(elmntlnknd(k,2))
                tpBy = ynd(elmntlnknd(k,2))
                tpBz = znd(elmntlnknd(k,2))
                tpCx = xnd(elmntlnknd(k,5))
                tpCy = ynd(elmntlnknd(k,5))
                tpCz = znd(elmntlnknd(k,5))
                tp1 = - tpAx*tpBy*tpCz + tpAx*tpCy*tpBz + tpBx*tpAy*tpCz &
                    & - tpBx*tpCy*tpAz - tpCx*tpAy*tpBz + tpCx*tpBy*tpAz
                tp2 = tp2 + tp1
                tpAx = xnd(elmntlnknd(k,6))
                tpAy = ynd(elmntlnknd(k,6))
                tpAz = znd(elmntlnknd(k,6))
                tpBx = xnd(elmntlnknd(k,5))
                tpBy = ynd(elmntlnknd(k,5))
                tpBz = znd(elmntlnknd(k,5))
                tpCx = xnd(elmntlnknd(k,3))
                tpCy = ynd(elmntlnknd(k,3))
                tpCz = znd(elmntlnknd(k,3))
                tp1 = - tpAx*tpBy*tpCz + tpAx*tpCy*tpBz + tpBx*tpAy*tpCz &
                    & - tpBx*tpCy*tpAz - tpCx*tpAy*tpBz + tpCx*tpBy*tpAz
                tp2 = tp2 + tp1
                tpAx = xnd(elmntlnknd(k,4))
                tpAy = ynd(elmntlnknd(k,4))
                tpAz = znd(elmntlnknd(k,4))
                tpBx = xnd(elmntlnknd(k,5))
                tpBy = ynd(elmntlnknd(k,5))
                tpBz = znd(elmntlnknd(k,5))
                tpCx = xnd(elmntlnknd(k,6))
                tpCy = ynd(elmntlnknd(k,6))
                tpCz = znd(elmntlnknd(k,6))
                tp1 = - tpAx*tpBy*tpCz + tpAx*tpCy*tpBz + tpBx*tpAy*tpCz &
                    & - tpBx*tpCy*tpAz - tpCx*tpAy*tpBz + tpCx*tpBy*tpAz
                tp2 = tp2 + tp1
            END DO
!$OMP END DO
!$OMP ATOMIC
            tp = tp + tp2
!$OMP END PARALLEL
            IF (NrmlInOut(ithprtl) == 1) volume(ithprtl) = tp/6.0d0
            IF (NrmlInOut(ithprtl) ==-1) volume(ithprtl) =-tp/6.0d0
        END DO

    END SUBROUTINE


    ! Evaluate first and second finite-difference derivatives of the local
    ! curvature auxiliary function along one Cartesian direction.

    SUBROUTINE dF_curv(dmdirect,dmithprtl,dmi,dmdF1,dmdF11,dmdF12)

        CHARACTER (LEN=1), INTENT (IN) :: dmdirect
        INTEGER, INTENT (IN) :: dmithprtl,dmi
        DOUBLE PRECISION, INTENT (OUT) :: dmdF1,dmdF11,dmdF12

        INTEGER :: fdm_nmx
        DOUBLE PRECISION :: dlh = 1.0d-6, dlh_1
        DOUBLE PRECISION, ALLOCATABLE :: fdxgrid(:), fdcff(:)
        DOUBLE PRECISION :: fdstt

        INTEGER :: i, j
        DOUBLE PRECISION :: tp,tp1,tp2,tp3,tp4,tp5,tp6

        dlh_1 = 1.0d0/dlh
        fdm_nmx = 5
        ALLOCATE (fdxgrid(fdm_nmx))
        ALLOCATE (fdcff(fdm_nmx))

        fdstt =-dble(((fdm_nmx-1)/2))
        DO i = 1, fdm_nmx
            fdxgrid(i) = fdstt + dble(i-1)
        END DO
        CALL fdcoef(2,fdm_nmx,fdxgrid((fdm_nmx+1)/2),fdxgrid(1),fdcff)



        dmdF1 = 0.0d0
        dmdF11= 0.0d0
        dmdF12= 0.0d0
        DO i = 1, fdm_nmx

            IF (dmdirect == 'x') THEN

                CALL F_curv(dmithprtl,dmi,0.0d0+fdxgrid(i)*dlh,0.0d0,0.0d0,tp1)
                dmdF1 = dmdF1 + tp1*fdcff(i)

                tp2 = 0.0d0
                DO j = 1, fdm_nmx
                    tp3 = fdxgrid(i)*dlh + fdxgrid(j)*dlh
                    CALL F_curv(dmithprtl,dmi,0.0d0+tp3,0.0d0,0.0d0,tp1)
                    tp2 = tp2 + tp1*fdcff(j)
                END DO
                dmdF11 = dmdF11 + (tp2 * dlh_1)*fdcff(i)

                tp2 = 0.0d0
                DO j = 1, fdm_nmx
                    tp3 = fdxgrid(j)*dlh
                    CALL F_curv(dmithprtl,dmi,0.0d0+fdxgrid(i)*dlh,0.0d0,0.0d0+tp3,tp1)
                    tp2 = tp2 + tp1*fdcff(j)
                END DO
                dmdF12 = dmdF12 + (tp2 * dlh_1)*fdcff(i)

            END IF

            IF (dmdirect == 'y') THEN

                CALL F_curv(dmithprtl,dmi,0.0d0,0.0d0+fdxgrid(i)*dlh,0.0d0,tp1)
                dmdF1 = dmdF1 + tp1*fdcff(i)

                tp2 = 0.0d0
                DO j = 1, fdm_nmx
                    tp3 = fdxgrid(i)*dlh + fdxgrid(j)*dlh
                    CALL F_curv(dmithprtl,dmi,0.0d0,0.0d0+tp3,0.0d0,tp1)
                    tp2 = tp2 + tp1*fdcff(j)
                END DO
                dmdF11 = dmdF11 + (tp2 * dlh_1)*fdcff(i)

                tp2 = 0.0d0
                DO j = 1, fdm_nmx
                    tp3 = fdxgrid(j)*dlh
                    CALL F_curv(dmithprtl,dmi,0.0d0,0.0d0+fdxgrid(i)*dlh,0.0d0+tp3,tp1)
                    tp2 = tp2 + tp1*fdcff(j)
                END DO
                dmdF12 = dmdF12 + (tp2 * dlh_1)*fdcff(i)

            END IF

            IF (dmdirect == 'z') THEN

                CALL F_curv(dmithprtl,dmi,0.0d0,0.0d0,0.0d0+fdxgrid(i)*dlh,tp1)
                dmdF1 = dmdF1 + tp1*fdcff(i)

                tp2 = 0.0d0
                DO j = 1, fdm_nmx
                    tp3 = fdxgrid(i)*dlh + fdxgrid(j)*dlh
                    CALL F_curv(dmithprtl,dmi,0.0d0,0.0d0,0.0d0+tp3,tp1)
                    tp2 = tp2 + tp1*fdcff(j)
                END DO
                dmdF11 = dmdF11 + (tp2 * dlh_1)*fdcff(i)

            END IF


        END DO
        dmdF1 = dmdF1 * dlh_1
        dmdF11= dmdF11* dlh_1
        dmdF12= dmdF12* dlh_1


        DEALLOCATE (fdxgrid,fdcff)

    END SUBROUTINE

    ! Evaluate the curvature auxiliary function at a small displacement from
    ! node dmi of object dmithprtl; used only by dF_curv.

    SUBROUTINE F_curv(dmithprtl,dmi,dmu,dmv,dmw,dmFcurv)

        INTEGER, INTENT (IN) :: dmithprtl,dmi
        DOUBLE PRECISION, INTENT (IN) :: dmu,dmv,dmw
        DOUBLE PRECISION, INTENT (OUT) :: dmFcurv
        DOUBLE PRECISION :: tpx,tpy,tpz,tpt1x,tpt1y,tpt1z,tpt2x,tpt2y,tpt2z,&
        &                   tpnnx,tpnny,tpnnz
        DOUBLE PRECISION :: tp,tp1,tp2,tp3,tp4,tp5,tp6
        DOUBLE PRECISION :: tpMtr11,tpMtr12,tpMtr13,tpMtr21,tpMtr22,tpMtr23,&
        &                   tpMtr31,tpMtr32,tpMtr33,tpMtrdet

        tpMtr11 = DCOS(anglecal_y(dmithprtl))*DCOS(anglecal_z(dmithprtl))
        tpMtr12 =-DCOS(anglecal_y(dmithprtl))*DSIN(anglecal_z(dmithprtl))
        tpMtr13 = DSIN(anglecal_y(dmithprtl))
        tpMtr21 = DSIN(anglecal_x(dmithprtl))*DSIN(anglecal_y(dmithprtl))&
        &        *DCOS(anglecal_z(dmithprtl)) &
        &        +DCOS(anglecal_x(dmithprtl))*DSIN(anglecal_z(dmithprtl))
        tpMtr22 =-DSIN(anglecal_x(dmithprtl))*DSIN(anglecal_y(dmithprtl))&
        &        *DSIN(anglecal_z(dmithprtl)) &
        &        +DCOS(anglecal_x(dmithprtl))*DCOS(anglecal_z(dmithprtl))
        tpMtr23 =-DSIN(anglecal_x(dmithprtl))*DCOS(anglecal_y(dmithprtl))
        tpMtr31 =-DCOS(anglecal_x(dmithprtl))*DSIN(anglecal_y(dmithprtl))&
        &        *DCOS(anglecal_z(dmithprtl)) &
        &        +DSIN(anglecal_x(dmithprtl))*DSIN(anglecal_z(dmithprtl))
        tpMtr32 = DCOS(anglecal_x(dmithprtl))*DSIN(anglecal_y(dmithprtl))&
        &        *DSIN(anglecal_z(dmithprtl)) &
        &        +DSIN(anglecal_x(dmithprtl))*DCOS(anglecal_z(dmithprtl))
        tpMtr33 = DCOS(anglecal_x(dmithprtl))*DCOS(anglecal_y(dmithprtl))
        tpMtrdet = tpMtr11*tpMtr22*tpMtr33&
                & +tpMtr21*tpMtr32*tpMtr13&
                & +tpMtr31*tpMtr12*tpMtr23&
                & -tpMtr11*tpMtr32*tpMtr23&
                & -tpMtr31*tpMtr22*tpMtr13&
                & -tpMtr21*tpMtr12*tpMtr33
        tpMtrdet = 1.0d0/tpMtrdet

        tpx = xnd(dmi) - xloctn(dmithprtl)
        tpy = ynd(dmi) - yloctn(dmithprtl)
        tpz = znd(dmi) - zloctn(dmithprtl)

        tp1 = tpMtrdet*(tpMtr22*tpMtr33-tpMtr23*tpMtr32)
        tp2 = tpMtrdet*(tpMtr13*tpMtr32-tpMtr12*tpMtr33)
        tp3 = tpMtrdet*(tpMtr12*tpMtr23-tpMtr13*tpMtr22)
        tp4 = tpx*tp1+tpy*tp2+tpz*tp3
        tpt1x = t1x(dmi)*tp1+t1y(dmi)*tp2+t1z(dmi)*tp3
        tpt2x = t2x(dmi)*tp1+t2y(dmi)*tp2+t2z(dmi)*tp3
        tpnnx = nnx(dmi)*tp1+nny(dmi)*tp2+nnz(dmi)*tp3

        tp1 = tpMtrdet*(tpMtr23*tpMtr31-tpMtr21*tpMtr33)
        tp2 = tpMtrdet*(tpMtr11*tpMtr33-tpMtr13*tpMtr31)
        tp3 = tpMtrdet*(tpMtr13*tpMtr21-tpMtr11*tpMtr23)
        tp5 = tpx*tp1+tpy*tp2+tpz*tp3
        tpt1y = t1x(dmi)*tp1+t1y(dmi)*tp2+t1z(dmi)*tp3
        tpt2y = t2x(dmi)*tp1+t2y(dmi)*tp2+t2z(dmi)*tp3
        tpnny = nnx(dmi)*tp1+nny(dmi)*tp2+nnz(dmi)*tp3

        tp1 = tpMtrdet*(tpMtr21*tpMtr32-tpMtr22*tpMtr31)
        tp2 = tpMtrdet*(tpMtr12*tpMtr31-tpMtr11*tpMtr32)
        tp3 = tpMtrdet*(tpMtr11*tpMtr22-tpMtr12*tpMtr21)
        tp6 = tpx*tp1+tpy*tp2+tpz*tp3
        tpt1z = t1x(dmi)*tp1+t1y(dmi)*tp2+t1z(dmi)*tp3
        tpt2z = t2x(dmi)*tp1+t2y(dmi)*tp2+t2z(dmi)*tp3
        tpnnz = nnx(dmi)*tp1+nny(dmi)*tp2+nnz(dmi)*tp3

        tpx = tp4/sizezoom(dmithprtl)
        tpy = tp5/sizezoom(dmithprtl)
        tpz = tp6/sizezoom(dmithprtl)

        tp4 = 1.0d0/DSQRT(tpt1x**2+tpt1y**2+tpt1z**2)
        tpt1x = tpt1x*tp4
        tpt1y = tpt1y*tp4
        tpt1z = tpt1z*tp4

        tp4 = 1.0d0/DSQRT(tpt2x**2+tpt2y**2+tpt2z**2)
        tpt2x = tpt2x*tp4
        tpt2y = tpt2y*tp4
        tpt2z = tpt2z*tp4

        tp4 = 1.0d0/DSQRT(tpnnx**2+tpnny**2+tpnnz**2)
        tpnnx = tpnnx*tp4
        tpnny = tpnny*tp4
        tpnnz = tpnnz*tp4

        tpx = tpt1x * dmu &
        &    +tpt2x * dmv &
        &    +tpnnx * dmw + tpx
        tpy = tpt1y * dmu &
        &    +tpt2y * dmv &
        &    +tpnny * dmw + tpy
        tpz = tpt1z * dmu &
        &    +tpt2z * dmv &
        &    +tpnnz * dmw + tpz

        IF (PrtlType(dmithprtl) == 'Sphr') THEN
            dmFcurv =-(tpx)**2-(tpy)**2-(tpz)**2
        END IF
        IF (PrtlType(dmithprtl) == 'PrSp') THEN
            dmFcurv =-(tpx/bvsa(dmithprtl))**2-(tpy/bvsa(dmithprtl))**2-(tpz)**2
        END IF
        IF (PrtlType(dmithprtl) == 'ObSp') THEN
            dmFcurv =-(tpx)**2-(tpy)**2-(tpz/bvsa(dmithprtl))**2
        END IF
        IF (PrtlType(dmithprtl) == 'DfSp') THEN
            dmFcurv =-(tpx/dfsp_a(dmithprtl))**dfsp_l(dmithprtl) &
            &        -(tpy/dfsp_b(dmithprtl))**dfsp_m(dmithprtl) &
            &        -(tpz/dfsp_c(dmithprtl))**dfsp_n(dmithprtl)
        END IF
        IF (PrtlType(dmithprtl) == 'DbBl') THEN
            tp = DSQRT((tpx)**2 + (tpy)**2)
            tp1= 1.0d0/((tpx)**2 + (tpy)**2)
            tp2= DSQRT((tpz+cvsa(dmithprtl))**2 + (tpx)**2 + (tpy)**2)
            tp3= DSQRT((tpz-cvsa(dmithprtl))**2 + (tpx)**2 + (tpy)**2)
            dmFcurv = (tp2)**(-3.0d0) &
            &        +(tp3)**(-3.0d0)
        END IF
        IF (PrtlType(dmithprtl) == 'Coin') THEN
            tp = DSQRT((tpx)**2 + (tpy)**2)
            tp2= DSQRT((tp+cvsa(dmithprtl))**2 + (tpz)**2)
            tp3= DSQRT((tp-cvsa(dmithprtl))**2 + (tpz)**2)
            dmFcurv = (tp2)**(-3.0d0) &
            &        +(tp3)**(-3.0d0)
        END IF
        IF (PrtlType(dmithprtl) == 'LgRd') THEN
            tp = DSQRT((tpx)**2 + (tpy)**2)
            tp1= 1.0d0/((tpx)**2 + (tpy)**2)
            tp2= DSQRT((tpz+cvsa(dmithprtl))**2 + (tpx)**2 + (tpy)**2)
            tp3= DSQRT((tpz-cvsa(dmithprtl))**2 + (tpx)**2 + (tpy)**2)
            dmFcurv = (1.0d0-cvsa(dmithprtl)**2)**2 &
            &           *((tpz+cvsa(dmithprtl))/tp2 - (tpz-cvsa(dmithprtl))/tp3) &
            &        -2.0d0*cvsa(dmithprtl) * ((tpx)**2 + (tpy)**2)
        END IF
        IF (PrtlType(dmithprtl) == 'Disk') THEN
            tp = DSQRT((tpx)**2 + (tpy)**2)
            tp1= 1.0d0/((tpz)**2)
            tp2= DSQRT((tp+cvsa(dmithprtl))**2 + (tpz)**2)
            tp3= DSQRT((tp-cvsa(dmithprtl))**2 + (tpz)**2)
            dmFcurv = (1.0d0-cvsa(dmithprtl)**2)**2 &
            &           *((tp+cvsa(dmithprtl))/tp2 - (tp-cvsa(dmithprtl))/tp3) &
            &        -2.0d0*cvsa(dmithprtl) * (tpz)**2
        END IF


    END SUBROUTINE



END MODULE
