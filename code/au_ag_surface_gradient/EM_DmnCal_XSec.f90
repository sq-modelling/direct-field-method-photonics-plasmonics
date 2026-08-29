
! SPDX-FileCopyrightText: 2026 Qiang Sun
! SPDX-License-Identifier: BSD-3-Clause

MODULE EM_DmnCal_XSec

    USE omp_lib

    USE Pre_Constants

    USE Geom_GlobalData
    USE Geom_MeshSphereCircle

    USE EM_SurfCal_GlobalData
    USE EM_SurfCal_PhysBC
    USE EM_DmnCal


    IMPLICIT NONE
    CONTAINS



    SUBROUTINE GetCrossSection_EM (R_crss,xcntr_crss,ycntr_crss,zcntr_crss,&
    &                              inc_crss,sct_crss,exc_crss,abs_crss)

        DOUBLE PRECISION, INTENT (IN) :: R_crss,xcntr_crss,ycntr_crss,zcntr_crss
        DOUBLE PRECISION, INTENT (OUT) :: inc_crss,sct_crss,exc_crss,abs_crss
        DOUBLE PRECISION :: ttl_crss

        CHARACTER (LEN=1) :: Mtype_crss = 'Q'

        DOUBLE PRECISION :: tpendx1,tpendx2,tpendx3,tpendx4,tpendx5,tpendx6
        DOUBLE PRECISION :: tpendy1,tpendy2,tpendy3,tpendy4,tpendy5,tpendy6
        DOUBLE PRECISION :: tpendz1,tpendz2,tpendz3,tpendz4,tpendz5,tpendz6
        DOUBLE PRECISION :: tpvctACx, tpvctACy, tpvctACz, &
        &                   tpvctBCx, tpvctBCy, tpvctBCz, &
        &                   tpvctrefCx, tpvctrefCy, tpvctrefCz, tpnrmlchck, &
        &                   tpsumndnnx,tpsumndnny,tpsumndnnz, tpsumndnnmdl, &
        &                   tpelmntnrmlx, tpelmntnrmly, tpelmntnrmlz, tpelmntareanrml, &
        &                   tpnx, tpny, tpnz

        DOUBLE PRECISION :: drx_deps,drx_dyet,dry_deps,dry_dyet,drz_deps,drz_dyet
        DOUBLE PRECISION :: tpEnx, tpEny, tpEnz, JcbDtmn
        DOUBLE PRECISION :: tpxieta

        DOUBLE PRECISION :: dphi_deps1,dphi_deps2,dphi_deps3,&
        &                   dphi_deps4,dphi_deps5,dphi_deps6
        DOUBLE PRECISION :: dphi_dyet1,dphi_dyet2,dphi_dyet3,&
        &                   dphi_dyet4,dphi_dyet5,dphi_dyet6
        DOUBLE PRECISION :: h_eps,h_yet,tpkappa,t_teps,t_tyet

        DOUBLE PRECISION :: pr0x, pr0y, pr0z, pr0nnx, pr0nny, pr0nnz, Dmnx, Dmny, Dmnz

        DOUBLE PRECISION :: tp,tp1,tp2,tp3,tp4,tp5,tp6,tp_a,tp_b,tp_g,tpxi,tpet

        INTEGER :: IOS, ithprtl, i, j, k, icmpnt, GLQi, icnt, id_tp
        INTEGER :: NdA, NdB, NdC, NdD, NdE, NdF, elmntA
        INTEGER :: rec_kmx, rec_k, mxnmbrndlnkelmntlnr_crss

        COMPLEX(KIND=KIND(1.0D0)) :: ztp, ztp1, ztp2, ztp3, ztp4, ztp5, ztpsum
        COMPLEX(KIND=KIND(1.0D0)) :: Exdmn,Eydmn,Ezdmn,Hxdmn,Hydmn,Hzdmn, &
        &                            ztpEx,ztpEy,ztpEz,ztpHx,ztpHy,ztpHz

        INTEGER :: nmbrnd_crss, nmbrelmnt_crss
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: x_crss, y_crss, z_crss
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: nnx_crss, nny_crss, nnz_crss
        INTEGER, ALLOCATABLE, DIMENSION (:,:) :: elmntlnknd_crss
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: incnd_crss,sctnd_crss, &
        &                                               ttlnd_crss,excnd_crss
        DOUBLE PRECISION :: dmIntI_capomp, tpareaomp, tparea, dmIntI_cap
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: dmE2_cap

        OPEN (9001, FILE = "Rslt_PrcdSmm.dat", STATUS = "OLD", &
            & POSITION="APPEND", ACTION="WRITE")
        WRITE (9001,*)
        WRITE (9001,*) '============================='
        WRITE (9001,*)
        WRITE (9001,*) 'cross section method - NRM'
        WRITE (9001,*)
        WRITE (9001,*) '============================='
        WRITE (9001,*)
        CLOSE (9001)

        tp1 = CDABS(exk_EM)
        tp = tp1 * R_crss

        tp = 5.0d0 + 1.0d0/1.950634514d0*(tp-1.950634514d0)
        icmpnt = NINT(tp)
        IF (MOD(icmpnt, 2) == 1) icmpnt = icmpnt + 1
        IF (icmpnt < 6) icmpnt = 6


!        CALL GetMeshNdElmntIcshdrlPaper(icmpnt,1.0d0,Mtype_crss)
        CALL GetMeshNdElmntIcshdrlSimple(icmpnt,1.0d0,Mtype_crss)

        OPEN (101, FILE = "Prtl_Orgnl.inp", STATUS = "OLD", IOSTAT = IOS)
        IF (IOS /= 0) THEN
            PRINT*, '"Prtl_Orgnl.inp" does not exist! Please check!'
            STOP
        END IF

        READ (101, *)
        READ (101, *)
        READ (101, *) nmbrnd_crss
        READ (101, *) nmbrelmnt_crss
        READ (101, *) tp1, tp2, tp3
        CLOSE (101, STATUS = 'DELETE')

        ALLOCATE (x_crss(nmbrnd_crss))
        ALLOCATE (y_crss(nmbrnd_crss))
        ALLOCATE (z_crss(nmbrnd_crss))
        ALLOCATE (nnx_crss(nmbrnd_crss))
        ALLOCATE (nny_crss(nmbrnd_crss))
        ALLOCATE (nnz_crss(nmbrnd_crss))

        IF (Mtype_crss == 'L') ALLOCATE (elmntlnknd_crss(nmbrelmnt_crss,3))
        IF (Mtype_crss == 'Q') ALLOCATE (elmntlnknd_crss(nmbrelmnt_crss,6))

        OPEN (111, FILE = "Prtl_Orgnl.vrt", STATUS = "OLD", IOSTAT = IOS)
        IF (IOS /= 0) THEN
            PRINT*, "'Prtl_Orgnl.vrt' does not exist! Please check!"
            STOP
        END IF

        DO i = 1, nmbrnd_crss
            READ (111, *) j, x_crss(i), y_crss(i), z_crss(i)
            nnx_crss(i) = x_crss(i)
            nny_crss(i) = y_crss(i)
            nnz_crss(i) = z_crss(i)
            x_crss(i) = x_crss(i) * R_crss
            y_crss(i) = y_crss(i) * R_crss
            z_crss(i) = z_crss(i) * R_crss
            x_crss(i) = x_crss(i) + xcntr_crss
            y_crss(i) = y_crss(i) + ycntr_crss
            z_crss(i) = z_crss(i) + zcntr_crss
        END DO

        CLOSE (111, STATUS = "DELETE")

        OPEN (121, FILE = "Prtl_Orgnl.cel", STATUS = "OLD", IOSTAT = IOS)
        IF (IOS /= 0) THEN
            PRINT*, '"Prtl_Orgnl.cel" does not exist! Please check!'
            STOP
        END IF

        DO k = 1, nmbrelmnt_crss
            IF (Mtype_crss == 'L') THEN
                READ (121, *) j, ndA, ndB, ndC
                elmntlnknd_crss(k,1) = ndA
                elmntlnknd_crss(k,2) = ndB
                elmntlnknd_crss(k,3) = ndC
            END IF
            IF (Mtype_crss == 'Q') THEN
                READ (121, *) j, ndA, ndB, ndC, ndD, ndE, ndF
                elmntlnknd_crss(k,1) = ndA
                elmntlnknd_crss(k,2) = ndB
                elmntlnknd_crss(k,3) = ndC
                elmntlnknd_crss(k,4) = ndD
                elmntlnknd_crss(k,5) = ndE
                elmntlnknd_crss(k,6) = ndF
            END IF
        END DO

        CLOSE (121, STATUS = "DELETE")

        ALLOCATE (incnd_crss(nmbrnd_crss))
        ALLOCATE (sctnd_crss(nmbrnd_crss))
        ALLOCATE (ttlnd_crss(nmbrnd_crss))
        ALLOCATE (excnd_crss(nmbrnd_crss))

!$OMP PARALLEL PRIVATE (i,ithprtl,k) &
!$OMP & PRIVATE (Dmnx,Dmny,Dmnz) &
!$OMP & PRIVATE (Exdmn,Eydmn,Ezdmn,Hxdmn,Hydmn,Hzdmn) &
!$OMP & PRIVATE (ztpEx,ztpEy,ztpEz,ztpHx,ztpHy,ztpHz) &
!$OMP & PRIVATE (ztp,ztp1,ztp2,ztp3)

!$OMP DO

        DO i = 1, nmbrnd_crss
            Dmnx = x_crss(i)
            Dmny = y_crss(i)
            Dmnz = z_crss(i)

            CALL GetDmnE1H1_EM (Dmnx,Dmny,Dmnz,0,ztpEx,ztpEy,ztpEz,ztpHx,ztpHy,ztpHz)
            ztp1 = ztpEy*dconjg(ztpHz) - ztpEz*dconjg(ztpHy)
            ztp2 = ztpEz*dconjg(ztpHx) - ztpEx*dconjg(ztpHz)
            ztp3 = ztpEx*dconjg(ztpHy) - ztpEy*dconjg(ztpHx)
            ztp = ztp1*nnx_crss(i)+ztp2*nny_crss(i)+ztp3*nnz_crss(i)
            incnd_crss(i) = 0.5d0 * REAL(ztp)

            CALL GetDmnE2H2_EM ('E',Dmnx,Dmny,Dmnz,0,Exdmn,Eydmn,Ezdmn)
            CALL GetDmnE2H2_EM ('H',Dmnx,Dmny,Dmnz,0,Hxdmn,Hydmn,Hzdmn)
            ztp1 = Eydmn*dconjg(Hzdmn) - Ezdmn*dconjg(Hydmn)
            ztp2 = Ezdmn*dconjg(Hxdmn) - Exdmn*dconjg(Hzdmn)
            ztp3 = Exdmn*dconjg(Hydmn) - Eydmn*dconjg(Hxdmn)
            ztp = ztp1*nnx_crss(i)+ztp2*nny_crss(i)+ztp3*nnz_crss(i)
            sctnd_crss(i) = 0.5d0 * REAL(ztp)

            ztp1 = Eydmn*dconjg(ztpHz) - Ezdmn*dconjg(ztpHy)
            ztp2 = Ezdmn*dconjg(ztpHx) - Exdmn*dconjg(ztpHz)
            ztp3 = Exdmn*dconjg(ztpHy) - Eydmn*dconjg(ztpHx)
            ztp = ztp1*nnx_crss(i)+ztp2*nny_crss(i)+ztp3*nnz_crss(i)
            excnd_crss(i) = 0.5d0 * REAL(ztp)
            ztp1 = ztpEy*dconjg(Hzdmn) - ztpEz*dconjg(Hydmn)
            ztp2 = ztpEz*dconjg(Hxdmn) - ztpEx*dconjg(Hzdmn)
            ztp3 = ztpEx*dconjg(Hydmn) - ztpEy*dconjg(Hxdmn)
            ztp = ztp1*nnx_crss(i)+ztp2*nny_crss(i)+ztp3*nnz_crss(i)
            excnd_crss(i) = excnd_crss(i) + 0.5d0 * REAL(ztp)

            Exdmn = Exdmn + ztpEx
            Eydmn = Eydmn + ztpEy
            Ezdmn = Ezdmn + ztpEz
            Hxdmn = Hxdmn + ztpHx
            Hydmn = Hydmn + ztpHy
            Hzdmn = Hzdmn + ztpHz
            ztp1 = Eydmn*dconjg(Hzdmn) - Ezdmn*dconjg(Hydmn)
            ztp2 = Ezdmn*dconjg(Hxdmn) - Exdmn*dconjg(Hzdmn)
            ztp3 = Exdmn*dconjg(Hydmn) - Eydmn*dconjg(Hxdmn)
            ztp = ztp1*nnx_crss(i)+ztp2*nny_crss(i)+ztp3*nnz_crss(i)
            ttlnd_crss(i) = 0.5d0 * REAL(ztp)

        END DO

!$OMP END DO
!$OMP END PARALLEL

        ALLOCATE(dmE2_cap(nmbrnd_crss))

        DO icmpnt = 1, 4

!$OMP PARALLEL PRIVATE (i)
!$OMP DO
            DO i = 1, nmbrnd_crss
                dmE2_cap(i) = 0.0d0
                IF (icmpnt == 1) dmE2_cap(i) = incnd_crss(i)
                IF (icmpnt == 2) dmE2_cap(i) = sctnd_crss(i)
                IF (icmpnt == 3) dmE2_cap(i) = excnd_crss(i)
                IF (icmpnt == 4) dmE2_cap(i) = ttlnd_crss(i)
            END DO
!$OMP END DO
!$OMP END PARALLEL

            IF (Mtype_crss == 'L') THEN

                dmIntI_cap = 0.0d0
                tparea = 0.0d0

!$OMP PARALLEL PRIVATE(k,GLQi,icnt,NdA,NdB,NdC) &
!$OMP & PRIVATE(tp,tp1,tp2,tp3,tpxi,tpet)&
!$OMP & PRIVATE(tpEnx,tpEny,tpEnz,JcbDtmn) &
!$OMP & PRIVATE(tpendx1,tpendx2,tpendx3) &
!$OMP & PRIVATE(tpendy1,tpendy2,tpendy3) &
!$OMP & PRIVATE(tpendz1,tpendz2,tpendz3) &
!$OMP & PRIVATE (dmIntI_capomp,tpareaomp)

                dmIntI_capomp = 0.0d0
                tpareaomp = 0.0d0

!$OMP DO
                DO k = 1, nmbrelmnt_crss

                    ndA = elmntlnknd_crss(k,1)
                    ndB = elmntlnknd_crss(k,2)
                    ndC = elmntlnknd_crss(k,3)

                    tpendx1 = x_crss(elmntlnknd_crss(k,1))
                    tpendx2 = x_crss(elmntlnknd_crss(k,2))
                    tpendx3 = x_crss(elmntlnknd_crss(k,3))

                    tpendy1 = y_crss(elmntlnknd_crss(k,1))
                    tpendy2 = y_crss(elmntlnknd_crss(k,2))
                    tpendy3 = y_crss(elmntlnknd_crss(k,3))

                    tpendz1 = z_crss(elmntlnknd_crss(k,1))
                    tpendz2 = z_crss(elmntlnknd_crss(k,2))
                    tpendz3 = z_crss(elmntlnknd_crss(k,3))

                    tpvctACx = x_crss(ndA) - x_crss(ndC)
                    tpvctACy = y_crss(ndA) - y_crss(ndC)
                    tpvctACz = z_crss(ndA) - z_crss(ndC)

                    tpvctBCx = x_crss(ndB) - x_crss(ndC)
                    tpvctBCy = y_crss(ndB) - y_crss(ndC)
                    tpvctBCz = z_crss(ndB) - z_crss(ndC)

                    tpEnx = tpvctACy*tpvctBCz - tpvctACz*tpvctBCy
                    tpEny = tpvctACz*tpvctBCx - tpvctACx*tpvctBCz
                    tpEnz = tpvctACx*tpvctBCy - tpvctACy*tpvctBCx
                    JcbDtmn = DSQRT( tpEnx**2+tpEny**2+tpEnz**2  )

                    DO GLQi = 1, n_glqtr2d

                        icnt = n_glqtr2d*(k-1)+GLQi

                        tpxi = xg_glqtr2d(GLQi)
                        tpet = yg_glqtr2d(GLQi)

                        tp2 = tpxi
                        tp3 = tpet
                        tp1 = 1.0d0 - tp2 - tp3

                        tp =   tp1*dmE2_cap(NdA) &
                        &     +tp2*dmE2_cap(NdB) &
                        &     +tp3*dmE2_cap(NdC)

                        dmIntI_capomp = dmIntI_capomp + wg_glqtr2d(GLQi) * tp * JcbDtmn

                        tp =   tp1+tp2+tp3

                        tpareaomp = tpareaomp + wg_glqtr2d(GLQi) * tp * JcbDtmn

                    END DO

                END DO
!$OMP END DO

!$OMP CRITICAL
                dmIntI_cap = dmIntI_cap + dmIntI_capomp
                tparea = tparea + tpareaomp
!$OMP END CRITICAL

!$OMP END PARALLEL

            END IF

            IF (Mtype_crss == 'Q') THEN

                dmIntI_cap = 0.0d0
                tparea = 0.0d0

!$OMP PARALLEL PRIVATE(k,GLQi,icnt,NdA,NdB,NdC,NdD,NdE,NdF) &
!$OMP & PRIVATE(tp,tp1,tp2,tp3,tp4,tp5,tp6,tp_a,tp_b,tp_g,tpxi,tpet)&
!$OMP & PRIVATE(tpEnx,tpEny,tpEnz,JcbDtmn) &
!$OMP & PRIVATE(tpendx1,tpendx2,tpendx3,tpendx4,tpendx5,tpendx6) &
!$OMP & PRIVATE(tpendy1,tpendy2,tpendy3,tpendy4,tpendy5,tpendy6) &
!$OMP & PRIVATE(tpendz1,tpendz2,tpendz3,tpendz4,tpendz5,tpendz6) &
!$OMP & PRIVATE(drx_deps,drx_dyet,dry_deps,dry_dyet,drz_deps,drz_dyet) &
!$OMP & PRIVATE (dmIntI_capomp,tpareaomp)

                dmIntI_capomp = 0.0d0
                tpareaomp = 0.0d0

!$OMP DO
                DO k = 1, nmbrelmnt_crss

                    NdA = elmntlnknd_crss(k,1)
                    NdB = elmntlnknd_crss(k,2)
                    NdC = elmntlnknd_crss(k,3)
                    NdD = elmntlnknd_crss(k,4)
                    NdE = elmntlnknd_crss(k,5)
                    NdF = elmntlnknd_crss(k,6)

                    tpendx1 = x_crss(elmntlnknd_crss(k,1))
                    tpendx2 = x_crss(elmntlnknd_crss(k,2))
                    tpendx3 = x_crss(elmntlnknd_crss(k,3))
                    tpendx4 = x_crss(elmntlnknd_crss(k,4))
                    tpendx5 = x_crss(elmntlnknd_crss(k,5))
                    tpendx6 = x_crss(elmntlnknd_crss(k,6))
                    tpendy1 = y_crss(elmntlnknd_crss(k,1))
                    tpendy2 = y_crss(elmntlnknd_crss(k,2))
                    tpendy3 = y_crss(elmntlnknd_crss(k,3))
                    tpendy4 = y_crss(elmntlnknd_crss(k,4))
                    tpendy5 = y_crss(elmntlnknd_crss(k,5))
                    tpendy6 = y_crss(elmntlnknd_crss(k,6))
                    tpendz1 = z_crss(elmntlnknd_crss(k,1))
                    tpendz2 = z_crss(elmntlnknd_crss(k,2))
                    tpendz3 = z_crss(elmntlnknd_crss(k,3))
                    tpendz4 = z_crss(elmntlnknd_crss(k,4))
                    tpendz5 = z_crss(elmntlnknd_crss(k,5))
                    tpendz6 = z_crss(elmntlnknd_crss(k,6))

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

                        tp =   tp1*dmE2_cap(NdA) &
                        &     +tp2*dmE2_cap(NdB) &
                        &     +tp3*dmE2_cap(NdC) &
                        &     +tp4*dmE2_cap(NdD) &
                        &     +tp5*dmE2_cap(NdE) &
                        &     +tp6*dmE2_cap(NdF)

                        dmIntI_capomp = dmIntI_capomp + wg_glqtr2d(GLQi) * tp * JcbDtmn

                        tp =   tp1+tp2+tp3+tp4+tp5+tp6

                        tpareaomp = tpareaomp + wg_glqtr2d(GLQi) * tp * JcbDtmn

                    END DO

                END DO

!$OMP END DO

!$OMP CRITICAL
                dmIntI_cap = dmIntI_cap + dmIntI_capomp
                tparea = tparea + tpareaomp
!$OMP END CRITICAL

!$OMP END PARALLEL

            END IF

            IF (icmpnt == 1) inc_crss =-dmIntI_cap
            IF (icmpnt == 2) sct_crss = dmIntI_cap
            IF (icmpnt == 3) exc_crss =-dmIntI_cap
            IF (icmpnt == 4) ttl_crss =-dmIntI_cap

        END DO

        tp = 1.0d0
        IF (     excitetype_EM == 'pwe' .OR. excitetype_EM == 'pwh' &
        &   .OR. excitetype_EM == 'swe' .OR. excitetype_EM == 'swh' ) THEN
            ztp = exeps_EM*exmiu_EM
            ztp = CDSQRT(ztp)
            tp  = 0.5d0*CDABS(ztp)*incFieldmdl_EM**2
            tp  = 1.0d0/tp
        END IF

        inc_crss = inc_crss*tp
        sct_crss = sct_crss*tp
        exc_crss = exc_crss*tp
        abs_crss = ttl_crss*tp

        DEALLOCATE (x_crss, y_crss, z_crss,nnx_crss,nny_crss,nnz_crss)
        DEALLOCATE (elmntlnknd_crss)
        DEALLOCATE (incnd_crss,sctnd_crss,excnd_crss,ttlnd_crss,dmE2_cap)

    END SUBROUTINE
END MODULE
