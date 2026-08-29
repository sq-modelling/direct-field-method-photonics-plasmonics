! SPDX-FileCopyrightText: 2026 Qiang Sun
! SPDX-License-Identifier: BSD-3-Clause

! Direct-field boundary-integral assembly and dense complex linear solves.
!
! Surface differential identities are evaluated in the local (t1,t2,n) frame.
! The divergence constraints use the signed curvature trace
! curvmn = curvt1 + curvt2 for the normal component and differentiate the
! complete tangential vectors rather than treating scalar components
! independently.  This form remains valid
! when the nodal tangent basis rotates or reverses between neighbouring nodes.
! The same construction is used for the electric and magnetic fields.
!
! Normal derivatives are solved as boundary unknowns.  They are not recovered
! from nearly singular post-processing integrals, which would be sensitive to
! arbitrary choices of local tangent-frame orientation on a high-order mesh.
! Legacy PEMC-specific branches remain isolated and are not exercised by the
! released PEC, dielectric, Au/Ag, or optical-force drivers.

MODULE EM_SurfCal_Solver

    USE omp_lib

    USE Pre_Constants
    USE Pre_csvformat

    USE Geom_GlobalData

    USE BRIEFGHComp

    USE EM_SurfCal_GlobalData

    IMPLICIT NONE

    INTEGER,PRIVATE :: PostProcdn_on = 0

    INTEGER,PRIVATE :: DivFree_option

    INTEGER,PRIVATE ::  cal_dn_sct
    INTEGER,PRIVATE ::  cal_dn_tra

    INTEGER :: rowid,colid,rowscl,colscl,bcscl=6

    COMPLEX(KIND=KIND(1.0D0)),PRIVATE, ALLOCATABLE, DIMENSION (:,:) :: DivFrAA
    COMPLEX(KIND=KIND(1.0D0)),PRIVATE, ALLOCATABLE, DIMENSION (:)   :: DivFrBB
    COMPLEX(KIND=KIND(1.0D0)),PRIVATE, ALLOCATABLE, DIMENSION (:,:) :: TanAAt1, TanAAt2
    COMPLEX(KIND=KIND(1.0D0)),PRIVATE, ALLOCATABLE, DIMENSION (:)   :: TanBBt1, TanBBt2

    COMPLEX(KIND=KIND(1.0D0)),PRIVATE, ALLOCATABLE, DIMENSION (:,:)   :: &
    &   bimGGEX, bimHHEX, bimGGIN, bimHHIN
    COMPLEX(KIND=KIND(1.0D0)),PRIVATE, ALLOCATABLE, DIMENSION (:,:)   :: &
    &   bimBMGGEX,bimBMHHEX,bimBMaxGGEX,bimax0GGEX,bimax0HHEX

    COMPLEX(KIND=KIND(1.0D0)),PRIVATE, ALLOCATABLE, DIMENSION (:,:)   :: bimAA,ztpPreCD,ztpPrercd
    COMPLEX(KIND=KIND(1.0D0)),PRIVATE, ALLOCATABLE, DIMENSION (:)     :: bimX,bimBB,ztpPreBB,ztpPreX
    COMPLEX(KIND=KIND(1.0D0)),PRIVATE, ALLOCATABLE, DIMENSION (:,:,:) :: bimPreCD
    INTEGER,PRIVATE :: bimPreCDsz, detCheck_On

    CONTAINS

    ! Assemble the surface-divergence equation for one PEC object.
    ! dmithprtl selects the object and dmcalfield selects the E or H trace.
    ! Coefficients are accumulated in the module-level dense system.

    SUBROUTINE Get_PEC_DivFr(dmithprtl, dmcalfield)

        INTEGER, INTENT(IN) :: dmithprtl
        CHARACTER(LEN=*), INTENT(IN) :: dmcalfield
        INTEGER :: i,j,ii,jj,kk,NdA,ithprtl,id_tp
        COMPLEX(KIND=KIND(1.0D0)) :: ztp,ztp1,ztp2,ztp3,mdm_io,mdm_oi

!*********
! Curvature form:
!div(exE2) = 0 <=> n•∂(exE2)/∂n = curvmn * exE2nn + curvmn * exE1nn - n•∂(exE1)/∂n
! Since exE3t1 = exE3t2 = 0 on PEC surface
! Unknowns are ordered as follows:
![exE2nn, t1•∂(exE2)/∂n, t2•∂(exE2)/∂n]

        IF (dmcalfield == 'H') RETURN

        ithprtl = dmithprtl
        IF (BCType_EM(ithprtl) == 'PEC') THEN
!$OMP PARALLEL PRIVATE (i,ii,jj,kk,NdA,id_tp,mdm_oi,mdm_io,ztp,ztp1,ztp2,ztp3)
!$OMP DO
            DO i = ndstaID(ithprtl), ndendID(ithprtl)

                !curvmn * E2n
                DivFrAA(i,3) = curvmn(i)

                !curvmn * E1n - n•∂E1/∂n
                ztp = nnx(i)*exE1x_EM(i) &
                &    +nny(i)*exE1y_EM(i) &
                &    +nnz(i)*exE1z_EM(i)
                DivFrBB(i) = DivFrBB(i) + ztp * curvmn(i)

                ztp = nnx(i)*exE1xdnn_EM(i) &
                &    +nny(i)*exE1ydnn_EM(i) &
                &    +nnz(i)*exE1zdnn_EM(i)
                DivFrBB(i) = DivFrBB(i) - ztp

            END DO
!$OMP END DO
!$OMP END PARALLEL
        END IF

    END SUBROUTINE

    ! Assemble the two tangential compatibility equations on one PEC surface.
    ! The full tangent vectors are differentiated so that the equations do
    ! not depend on an arbitrary orientation of the local tangent basis.

    SUBROUTINE Get_PEC_TanCn(dmithprtl, dmcalfield)

        INTEGER, INTENT(IN) :: dmithprtl
        CHARACTER(LEN=*), INTENT(IN) :: dmcalfield
        INTEGER :: i,j,ii,jj,kk,NdA,ithprtl,id_tp
        COMPLEX(KIND=KIND(1.0D0)) :: ztp,ztp1,ztp2,ztp3,mdm_io,mdm_oi

!*********
!d/dt1, d/dt2:
!t1•∂(exH2)/∂n = n*∂(exH2t1*t1)/∂t1 + n*∂(exH2t2*t2)/∂t1 &
!&              -∂(exH1n)/∂t1 + n*∂(exH1)/∂t1 - t1*∂(exH1)/∂n
!t2•∂(exH2)/∂n = n*∂(exH2t1*t1)/∂t2 + n*∂(exH2t2*t2)/∂t2 &
!&              -∂(exH1n)/∂t2 + n*∂(exH1)/∂t2 - t2*∂(exH1)/∂n
!
! Unknowns are ordered as follows:
![exH2t1, exH2t2, n•∂(exH2)/∂n]

        IF (dmcalfield == 'E') RETURN

        ithprtl = dmithprtl
        IF (BCType_EM(ithprtl) == 'PEC') THEN
!$OMP PARALLEL PRIVATE (i,ii,jj,kk,NdA,id_tp,mdm_oi,mdm_io,ztp,ztp1,ztp2,ztp3)
!$OMP DO
            DO i = ndstaID(ithprtl), ndendID(ithprtl)

                !n*∂(exH2t1*t1)/∂t1 + n*∂(exH2t2*t2)/∂t1
                !n*∂(exH2t1*t1)/∂t2 + n*∂(exH2t2*t2)/∂t2
                ztp1 = ztpzero
                ztp2 = ztpzero
                DO kk = 1, mxnmbrndlnknd2ndslf
                    IF (ndlnknd2ndslf(i,kk) /= 0) THEN

                        NdA = ndlnknd2ndslf(i,kk)
                        ii = bcscl*(kk-1)

                        TanAAt1(i,ii+1) = TanAAt1(i,ii+1) &
                        &                +                ( nnx(i)*t1x(NdA) &
                        &                                  +nny(i)*t1y(NdA) &
                        &                                  +nnz(i)*t1z(NdA) )*d_dt1(i,kk)
                        TanAAt1(i,ii+2) = TanAAt1(i,ii+2) &
                        &                +                ( nnx(i)*t2x(NdA) &
                        &                                  +nny(i)*t2y(NdA) &
                        &                                  +nnz(i)*t2z(NdA) )*d_dt1(i,kk)

                        TanAAt2(i,ii+1) = TanAAt2(i,ii+1) &
                        &                +                ( nnx(i)*t1x(NdA) &
                        &                                  +nny(i)*t1y(NdA) &
                        &                                  +nnz(i)*t1z(NdA) )*d_dt2(i,kk)
                        TanAAt2(i,ii+2) = TanAAt2(i,ii+2) &
                        &                +                ( nnx(i)*t2x(NdA) &
                        &                                  +nny(i)*t2y(NdA) &
                        &                                  +nnz(i)*t2z(NdA) )*d_dt2(i,kk)

                        !-∂(exH1n)/∂t1
                        !-∂(exH1n)/∂t2
                        IF (dmcalfield == 'H') THEN
                            ztp = nnx(NdA)*exH1x_EM(NdA)+nny(NdA)*exH1y_EM(NdA)+nnz(NdA)*exH1z_EM(NdA)
                            ztp1 = ztp1 - ztp*d_dt1(i,kk)
                            ztp2 = ztp2 - ztp*d_dt2(i,kk)
                        END IF



                    END IF
                END DO

                TanBBt1(i) = TanBBt1(i) + ztp1
                TanBBt2(i) = TanBBt2(i) + ztp2

!&              n*∂(exH1)/∂t1 - t1*∂(exH1)/∂n
                ztp1 = nnx(i)*exH1xdt1_EM(i) &
                &     +nny(i)*exH1ydt1_EM(i) &
                &     +nnz(i)*exH1zdt1_EM(i)
                ztp2 = t1x(i)*exH1xdnn_EM(i) &
                &     +t1y(i)*exH1ydnn_EM(i) &
                &     +t1z(i)*exH1zdnn_EM(i)

                TanBBt1(i) = TanBBt1(i) + ztp1 - ztp2

!&              n*∂(exH1)/∂t2 - t2*∂(exH1)/∂n
                ztp1 = nnx(i)*exH1xdt2_EM(i) &
                &     +nny(i)*exH1ydt2_EM(i) &
                &     +nnz(i)*exH1zdt2_EM(i)
                ztp2 = t2x(i)*exH1xdnn_EM(i) &
                &     +t2y(i)*exH1ydnn_EM(i) &
                &     +t2z(i)*exH1zdnn_EM(i)

                TanBBt2(i) = TanBBt2(i) + ztp1 - ztp2

            END DO
!$OMP END DO
!$OMP END PARALLEL
        END IF

    END SUBROUTINE





    ! Assemble the divergence equation across one dielectric interface.
    ! Exterior scattered and interior transmitted traces are coupled through
    ! the material and orientation conventions stored in the global arrays.

    SUBROUTINE Get_DEL_DivFrsc(dmithprtl, dmcalfield)

        INTEGER, INTENT(IN) :: dmithprtl
        CHARACTER(LEN=*), INTENT(IN) :: dmcalfield
        INTEGER :: i,j,ii,jj,kk,NdA,ithprtl,id_tp
        COMPLEX(KIND=KIND(1.0D0)) :: ztp,ztp1,ztp2,ztp3,mdm_io,mdm_oi

!*********
!curvature
!divE = 0 on both sides <=> n•∂(inE2)/∂n = n•∂(exE2)/∂n + (esp_oi-1) * curvmn * exE2n &
!&                              + n•∂(exE1)/∂n + (esp_oi-1) * curvmn * exE1n - n•∂(inE1)/∂n
!divH = 0 on both sides <=> n•∂(inH2)/∂n = n•∂(exH2)/∂n + (miu_oi-1) * curvmn * exH2n &
!&                              + n•∂(exH1)/∂n + (miu_oi-1) * curvmn * exH1n - n•∂(inH1)/∂n
!
! Unknowns are ordered as follows:
![exE2t1, exE2t2, exE2nn, t1•∂(exE2)/∂n, t2•∂(exE2)/∂n, n•∂(exE2)/∂n]
![exH2t1, exH2t2, exH2nn, t1•∂(exH2)/∂n, t2•∂(exH2)/∂n, n•∂(exH2)/∂n]

        ithprtl = dmithprtl
        IF (BCType_EM(ithprtl) == '2SD') THEN
!$OMP PARALLEL PRIVATE (i,ii,jj,kk,NdA,id_tp,mdm_oi,mdm_io,ztp,ztp1,ztp2,ztp3)
!$OMP DO
            DO i = ndstaID(ithprtl), ndendID(ithprtl)
                id_tp = corelnkshell(ithprtl)
                IF (id_tp == 0) THEN
                    IF (dmcalfield == 'E') mdm_oi = exeps_EM/ineps_EM(ithprtl)
                    IF (dmcalfield == 'H') mdm_oi = exmiu_EM/inmiu_EM(ithprtl)
                ELSE
                    IF (dmcalfield == 'E') mdm_oi = ineps_EM(id_tp)/ineps_EM(ithprtl)
                    IF (dmcalfield == 'H') mdm_oi = inmiu_EM(id_tp)/inmiu_EM(ithprtl)
                END IF

                !(esp_oi-1) * curvmn * En_sct
                !(mu_oi -1) * curvmn * Hn_sct
                DivFrAA(i,3) = curvmn(i) * (mdm_oi - 1.0d0)
                DivFrAA(i,6) = 1.0d0

                !n•∂(exE1)/∂n + (esp_oi-1) * curvmn * exE1n - n•∂(inE1)/∂n
                IF (dmcalfield == 'E') THEN
                    ztp = nnx(i)*exE1x_EM(i) &
                    &    +nny(i)*exE1y_EM(i) &
                    &    +nnz(i)*exE1z_EM(i)
                    DivFrBB(i) = DivFrBB(i) + ztp * curvmn(i) * (mdm_oi - 1.0d0)

                    ztp = nnx(i)*exE1xdnn_EM(i) &
                    &    +nny(i)*exE1ydnn_EM(i) &
                    &    +nnz(i)*exE1zdnn_EM(i)
                    DivFrBB(i) = DivFrBB(i) + ztp

                    ztp = nnx(i)*inE1xdnn_EM(i) &
                    &    +nny(i)*inE1ydnn_EM(i) &
                    &    +nnz(i)*inE1zdnn_EM(i)
                    DivFrBB(i) = DivFrBB(i) - ztp
                END IF

                !n•∂(exH1)/∂n + (mu_oi-1)  * curvmn * exH1n - n•∂(inH1)/∂n
                IF (dmcalfield == 'H') THEN
                    ztp = nnx(i)*exH1x_EM(i) &
                    &    +nny(i)*exH1y_EM(i) &
                    &    +nnz(i)*exH1z_EM(i)
                    DivFrBB(i) = DivFrBB(i) + ztp * curvmn(i) * (mdm_oi - 1.0d0)

                    ztp = nnx(i)*exH1xdnn_EM(i) &
                    &    +nny(i)*exH1ydnn_EM(i) &
                    &    +nnz(i)*exH1zdnn_EM(i)
                    DivFrBB(i) = DivFrBB(i) + ztp

                    ztp = nnx(i)*inH1xdnn_EM(i) &
                    &    +nny(i)*inH1ydnn_EM(i) &
                    &    +nnz(i)*inH1zdnn_EM(i)
                    DivFrBB(i) = DivFrBB(i) - ztp
                END IF

            END DO
!$OMP END DO
!$OMP END PARALLEL
        END IF

    END SUBROUTINE

    ! Assemble tangential derivative relations implied by dielectric field
    ! continuity for the selected interface and field (E or H).

    SUBROUTINE Get_DEL_TanCnsc(dmithprtl, dmcalfield)

        INTEGER, INTENT(IN) :: dmithprtl
        CHARACTER(LEN=*), INTENT(IN) :: dmcalfield
        INTEGER :: i,j,ii,jj,kk,NdA,ithprtl,id_tp
        COMPLEX(KIND=KIND(1.0D0)) :: ztp,ztp1,ztp2,ztp3,mdm_io,mdm_oi

!*********
!d/dt:
!t1•∂(inE2)/∂n = miu_io*t1•∂(exE2)/∂n &
!&              +(eps_oi-miu_io)*∂(exE2n)/∂t1 + (1-miu_io)*n*∂(exE2t1*t1)/∂t1 + (1-miu_io)*n*∂(exE2t2*t2)/∂t1 &
!&              +(eps_oi-1)*∂(exE1n)/∂t1 + (1-miu_io)*n*∂(exE1)/∂t1 + miu_io*t1*∂(exE1)/∂n &
!&              -t1•∂(inE1)/∂n
!t2•∂(inE2)/∂n = miu_io*t2•∂(exE2)/∂n &
!&              +(eps_oi-miu_io)*∂(exE2n)/∂t2 + (1-miu_io)*n*∂(exE2t1*t1)/∂t2 + (1-miu_io)*n*∂(exE2t2*t2)/∂t2 &
!&              +(eps_oi-1)*∂(exE1n)/∂t2 + (1-miu_io)*n*∂(exE1)/∂t2 + miu_io*t2*∂(exE1)/∂n &
!&              -t2•∂(inE1)/∂n

!t1•∂(inH2)/∂n = eps_io*t1•∂(exH2)/∂n &
!&              +(miu_oi-eps_io)*∂(exH2n)/∂t1 + (1-eps_io)*n*∂(exH2t1*t1)/∂t1 + (1-eps_io)*n*∂(exH2t2*t2)/∂t1 &
!&              +(miu_oi-1)*∂(exH1n)/∂t1 + (1-eps_io)*n*∂(exH1)/∂t1 + eps_io*t1*∂(exH1)/∂n &
!&              -t1•∂(inH1)/∂n
!t2•∂(inH2)/∂n = eps_io*t2•∂(exH2)/∂n &
!&              +(miu_oi-eps_io)*∂(exH2n)/∂t2 + (1-eps_io)*n*∂(exH2t1*t1)/∂t2 + (1-eps_io)*n*∂(exH2t2*t2)/∂t2 &
!&              +(miu_oi-1)*∂(exH1n)/∂t2 + (1-eps_io)*n*∂(exH1)/∂t2 + eps_io*t2*∂(exH1)/∂n &
!&              -t2•∂(inH1)/∂n
!
! Unknowns are ordered as follows:
![exE2t1, exE2t2, exE2nn, t1•∂(exE2)/∂n, t2•∂(exE2)/∂n, n•∂(exE2)/∂n]
![exH2t1, exH2t2, exH2nn, t1•∂(exH2)/∂n, t2•∂(exH2)/∂n, n•∂(exH2)/∂n]

        ithprtl = dmithprtl
        IF (BCType_EM(ithprtl) == '2SD') THEN
!$OMP PARALLEL PRIVATE (i,ii,jj,kk,NdA,id_tp,mdm_oi,mdm_io,ztp,ztp1,ztp2,ztp3)
!$OMP DO
            DO i = ndstaID(ithprtl), ndendID(ithprtl)

                id_tp = corelnkshell(ithprtl)
                IF (id_tp == 0) THEN
                    IF (dmcalfield == 'E') THEN
                        mdm_io = inmiu_EM(ithprtl)/exmiu_EM
                        mdm_oi = exeps_EM/ineps_EM(ithprtl)
                    END IF
                    IF (dmcalfield == 'H') THEN
                        mdm_io = ineps_EM(ithprtl)/exeps_EM
                        mdm_oi = exmiu_EM/inmiu_EM(ithprtl)
                    END IF
                ELSE
                    IF (dmcalfield == 'E') THEN
                        mdm_io = inmiu_EM(ithprtl)/inmiu_EM(id_tp)
                        mdm_oi = ineps_EM(id_tp)/ineps_EM(ithprtl)
                    END IF
                    IF (dmcalfield == 'H') THEN
                        mdm_io = ineps_EM(ithprtl)/ineps_EM(id_tp)
                        mdm_oi = inmiu_EM(id_tp)/inmiu_EM(ithprtl)
                    END IF
                END IF

                !miu_io * t1•∂(exE2)/∂n
                !eps_io * t1•∂(exH2)/∂n
                TanAAt1(i,4) = mdm_io

                !miu_io * t2•∂(exE2)/∂n
                !eps_io * t2•∂(exH2)/∂n
                TanAAt2(i,5) = mdm_io

                !(eps_oi-miu_io)*∂(exE2n)/∂t1 + (1-miu_io)*n*∂(exE2t1*t1)/∂t1 + (1-miu_io)*n*∂(exE2t2*t2)/∂t1
                !(eps_oi-miu_io)*∂(exE2n)/∂t2 + (1-miu_io)*n*∂(exE2t1*t1)/∂t2 + (1-miu_io)*n*∂(exE2t2*t2)/∂t2
                !(miu_oi-eps_io)*∂(exE2n)/∂t1 + (1-eps_io)*n*∂(exE2t1*t1)/∂t1 + (1-eps_io)*n*∂(exE2t2*t2)/∂t1
                !(miu_oi-eps_io)*∂(exE2n)/∂t2 + (1-eps_io)*n*∂(exE2t1*t1)/∂t2 + (1-eps_io)*n*∂(exE2t2*t2)/∂t2

                ztp1 = ztpzero
                ztp2 = ztpzero
                DO kk = 1, mxnmbrndlnknd2ndslf
                    IF (ndlnknd2ndslf(i,kk) /= 0) THEN

                        NdA = ndlnknd2ndslf(i,kk)
                        ii = bcscl*(kk-1)

                        TanAAt1(i,ii+1) = TanAAt1(i,ii+1) &
                        &                +(1.0d0-mdm_io)* ( nnx(i)*t1x(NdA) &
                        &                                  +nny(i)*t1y(NdA) &
                        &                                  +nnz(i)*t1z(NdA) )*d_dt1(i,kk)
                        TanAAt1(i,ii+2) = TanAAt1(i,ii+2) &
                        &                +(1.0d0-mdm_io)* ( nnx(i)*t2x(NdA) &
                        &                                  +nny(i)*t2y(NdA) &
                        &                                  +nnz(i)*t2z(NdA) )*d_dt1(i,kk)
                        TanAAt1(i,ii+3) = TanAAt1(i,ii+3) &
                        &                +(mdm_oi-mdm_io)*d_dt1(i,kk)

                        TanAAt2(i,ii+1) = TanAAt2(i,ii+1) &
                        &                +(1.0d0-mdm_io)* ( nnx(i)*t1x(NdA) &
                        &                                  +nny(i)*t1y(NdA) &
                        &                                  +nnz(i)*t1z(NdA) )*d_dt2(i,kk)
                        TanAAt2(i,ii+2) = TanAAt2(i,ii+2) &
                        &                +(1.0d0-mdm_io)* ( nnx(i)*t2x(NdA) &
                        &                                  +nny(i)*t2y(NdA) &
                        &                                  +nnz(i)*t2z(NdA) )*d_dt2(i,kk)
                        TanAAt2(i,ii+3) = TanAAt2(i,ii+3) &
                        &                +(mdm_oi-mdm_io)*d_dt2(i,kk)

                        !(eps_oi-1)*∂(exE1n)/∂t1
                        !(eps_oi-1)*∂(exE1n)/∂t2
                        !(miu_oi-1)*∂(exH1n)/∂t1
                        !(miu_oi-1)*∂(exH1n)/∂t2
                        IF (dmcalfield == 'E') THEN
                            ztp = nnx(NdA)*exE1x_EM(NdA)+nny(NdA)*exE1y_EM(NdA)+nnz(NdA)*exE1z_EM(NdA)
                            ztp1 = ztp1 &
                            &     +(mdm_oi-1.0d0)*ztp*d_dt1(i,kk)
                            ztp2 = ztp2 &
                            &     +(mdm_oi-1.0d0)*ztp*d_dt2(i,kk)
                        END IF
                        IF (dmcalfield == 'H') THEN
                            ztp = nnx(NdA)*exH1x_EM(NdA)+nny(NdA)*exH1y_EM(NdA)+nnz(NdA)*exH1z_EM(NdA)
                            ztp1 = ztp1 &
                            &     +(mdm_oi-1.0d0)*ztp*d_dt1(i,kk)
                            ztp2 = ztp2 &
                            &     +(mdm_oi-1.0d0)*ztp*d_dt2(i,kk)
                        END IF



                    END IF
                END DO

                TanBBt1(i) = TanBBt1(i) + ztp1
                TanBBt2(i) = TanBBt2(i) + ztp2

!&              (1-miu_io)*n*∂(exE1)/∂t1 + miu_io*t1*∂(exE1)/∂n &
!&              -t1•∂(inE1)/∂n
!&              (1-eps_io)*n*∂(exH1)/∂t1 + eps_io*t1*∂(exH1)/∂n &
!&              -t1•∂(inH1)/∂n
                IF (dmcalfield == 'E') THEN
                    ztp1 = nnx(i)*exE1xdt1_EM(i) &
                    &     +nny(i)*exE1ydt1_EM(i) &
                    &     +nnz(i)*exE1zdt1_EM(i)
                    ztp2 = t1x(i)*exE1xdnn_EM(i) &
                    &     +t1y(i)*exE1ydnn_EM(i) &
                    &     +t1z(i)*exE1zdnn_EM(i)
                    ztp3 = t1x(i)*inE1xdnn_EM(i) &
                    &     +t1y(i)*inE1ydnn_EM(i) &
                    &     +t1z(i)*inE1zdnn_EM(i)
                END IF
                IF (dmcalfield == 'H') THEN
                    ztp1 = nnx(i)*exH1xdt1_EM(i) &
                    &     +nny(i)*exH1ydt1_EM(i) &
                    &     +nnz(i)*exH1zdt1_EM(i)
                    ztp2 = t1x(i)*exH1xdnn_EM(i) &
                    &     +t1y(i)*exH1ydnn_EM(i) &
                    &     +t1z(i)*exH1zdnn_EM(i)
                    ztp3 = t1x(i)*inH1xdnn_EM(i) &
                    &     +t1y(i)*inH1ydnn_EM(i) &
                    &     +t1z(i)*inH1zdnn_EM(i)
                END IF

                TanBBt1(i) = TanBBt1(i) + (1.0d0 - mdm_io) * ztp1 &
                &                       +          mdm_io  * ztp2 &
                &                       -                    ztp3

!&              (1-miu_io)*n*∂(exE1)/∂t2 + miu_io*t2*∂(exE1)/∂n &
!&              -t2•∂(inE2)/∂n
!&              (1-eps_io)*n*∂(exH1)/∂t2 + eps_io*t2*∂(exH1)/∂n &
!&              -t2•∂(inH1)/∂n
                IF (dmcalfield == 'E') THEN
                    ztp1 = nnx(i)*exE1xdt2_EM(i) &
                    &     +nny(i)*exE1ydt2_EM(i) &
                    &     +nnz(i)*exE1zdt2_EM(i)
                    ztp2 = t2x(i)*exE1xdnn_EM(i) &
                    &     +t2y(i)*exE1ydnn_EM(i) &
                    &     +t2z(i)*exE1zdnn_EM(i)
                    ztp3 = t2x(i)*inE1xdnn_EM(i) &
                    &     +t2y(i)*inE1ydnn_EM(i) &
                    &     +t2z(i)*inE1zdnn_EM(i)
                END IF
                IF (dmcalfield == 'H') THEN
                    ztp1 = nnx(i)*exH1xdt2_EM(i) &
                    &     +nny(i)*exH1ydt2_EM(i) &
                    &     +nnz(i)*exH1zdt2_EM(i)
                    ztp2 = t2x(i)*exH1xdnn_EM(i) &
                    &     +t2y(i)*exH1ydnn_EM(i) &
                    &     +t2z(i)*exH1zdnn_EM(i)
                    ztp3 = t2x(i)*inH1xdnn_EM(i) &
                    &     +t2y(i)*inH1ydnn_EM(i) &
                    &     +t2z(i)*inH1zdnn_EM(i)
                END IF

                TanBBt2(i) = TanBBt2(i) + (1.0d0 - mdm_io) * ztp1 &
                &                       +          mdm_io  * ztp2 &
                &                       -                    ztp3

            END DO
!$OMP END DO
!$OMP END PARALLEL
        END IF

    END SUBROUTINE




    ! Assemble the regularised Helmholtz G and H influence matrices for
    ! bounded material domains on the current linear or Q6 surface mesh.

    SUBROUTINE Get_BRIEFGGHHIN_EM

        INTEGER :: ithprtl, jthprtl, kthprtl, i, j, k
        INTEGER :: NdA, NdB, NdC, NdD, NdE, NdF, id_tp

        INTEGER :: Ttlbim

        DOUBLE PRECISION :: pr0x, pr0y, pr0z, pr0nnx, pr0nny, pr0nnz

        COMPLEX(KIND=KIND(1.0D0)) :: g1, g2, g3, g4, g5, g6
        COMPLEX(KIND=KIND(1.0D0)) :: h1, h2, h3, h4, h5, h6
        COMPLEX(KIND=KIND(1.0D0)) :: gnsbim, hnsbim

        COMPLEX(KIND=KIND(1.0D0)) :: bimGGINNS,bimHHINNS, ztp

!*********

        Ttlbim = ttlnmbrnd

!*********
!BRIEF - G & H

        ALLOCATE (bimGGIN(Ttlbim,Ttlbim))
        ALLOCATE (bimHHIN(Ttlbim,Ttlbim))

!$OMP PARALLEL PRIVATE (i,j)
!$OMP DO
        DO i = 1, Ttlbim
            DO j = 1, Ttlbim
                bimGGIN(i,j) = ztpzero
                bimHHIN(i,j) = ztpzero
            END DO
        END DO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE (i,k,ithprtl,jthprtl,id_tp) &
!$OMP & PRIVATE (NdA,NdB,NdC,NdD,NdE,NdF) &
!$OMP & PRIVATE (pr0x,pr0y,pr0z,pr0nnx,pr0nny,pr0nnz) &
!$OMP & PRIVATE (bimGGINNS,bimHHINNS,ztp) &
!$OMP & PRIVATE (g1,g2,g3,g4,g5,g6,h1,h2,h3,h4,h5,h6,gnsbim,hnsbim)
!$OMP DO
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

            pr0x = xnd(i)
            pr0y = ynd(i)
            pr0z = znd(i)
            pr0nnx = nnx(i)
            pr0nny = nny(i)
            pr0nnz = nnz(i)

!***
! Assemble the interior-domain G and H operators.
! Include the host surface and surfaces nested one level inside it.

            bimGGINNS = ztpzero
            bimHHINNS = ztpzero

            DO k = 1, ttlnmbrelmnt

                jthprtl = 1
                IF (nmbrprtl > 1 .AND. k > elendID(1)) THEN
                    DO id_tp = 2, nmbrprtl
                        IF (k >= elstaID(id_tp) .AND. k <= elendID(id_tp)) THEN
                            jthprtl = id_tp
                            EXIT
                        END IF
                    END DO
                END IF

                IF (corelnkshell(jthprtl) == ithprtl .OR. jthprtl == ithprtl) THEN

                    IF (MeshType == "L") THEN

                        NdA = elmntlnknd(k,1)
                        NdB = elmntlnknd(k,2)
                        NdC = elmntlnknd(k,3)

                        CALL CalGHLnrBRIEFLnrCOMP( ink_EM(ithprtl), k, &
                        &       pr0x, pr0y, pr0z, pr0nnx, pr0nny, pr0nnz, &
                        &       g1, g2, g3, h1, h2, h3, gnsbim, hnsbim)

                        bimGGIN(i,NdA) = bimGGIN(i,NdA) + g1
                        bimGGIN(i,NdB) = bimGGIN(i,NdB) + g2
                        bimGGIN(i,NdC) = bimGGIN(i,NdC) + g3

                        bimHHIN(i,NdA) = bimHHIN(i,NdA) + h1
                        bimHHIN(i,NdB) = bimHHIN(i,NdB) + h2
                        bimHHIN(i,NdC) = bimHHIN(i,NdC) + h3

                        bimGGINNS = bimGGINNS + gnsbim
                        bimHHINNS = bimHHINNS + hnsbim

                   END IF


                    IF (MeshType == "Q") THEN

                        NdA = elmntlnknd(k,1)
                        NdB = elmntlnknd(k,2)
                        NdC = elmntlnknd(k,3)
                        NdD = elmntlnknd(k,4)
                        NdE = elmntlnknd(k,5)
                        NdF = elmntlnknd(k,6)

                        CALL CalGHQdrBRIEFLnrCOMP( ink_EM(ithprtl), k, &
                        &       pr0x, pr0y, pr0z, pr0nnx, pr0nny, pr0nnz, &
                        &       g1, g2, g3, g4, g5, g6, &
                        &       h1, h2, h3, h4, h5, h6, gnsbim, hnsbim)

                        bimGGIN(i,NdA) = bimGGIN(i,NdA) + g1
                        bimGGIN(i,NdB) = bimGGIN(i,NdB) + g2
                        bimGGIN(i,NdC) = bimGGIN(i,NdC) + g3
                        bimGGIN(i,NdD) = bimGGIN(i,NdD) + g4
                        bimGGIN(i,NdE) = bimGGIN(i,NdE) + g5
                        bimGGIN(i,NdF) = bimGGIN(i,NdF) + g6

                        bimHHIN(i,NdA) = bimHHIN(i,NdA) + h1
                        bimHHIN(i,NdB) = bimHHIN(i,NdB) + h2
                        bimHHIN(i,NdC) = bimHHIN(i,NdC) + h3
                        bimHHIN(i,NdD) = bimHHIN(i,NdD) + h4
                        bimHHIN(i,NdE) = bimHHIN(i,NdE) + h5
                        bimHHIN(i,NdF) = bimHHIN(i,NdF) + h6

                        bimGGINNS = bimGGINNS + gnsbim
                        bimHHINNS = bimHHINNS + hnsbim

                    END IF

                END IF

            END DO

            bimGGIN(i,i) = bimGGIN(i,i) + bimGGINNS
            bimHHIN(i,i) = bimHHIN(i,i) + bimHHINNS

!End internal GG and HH
!***
        END DO
!$OMP END DO
!$OMP END PARALLEL

    END SUBROUTINE

    ! Assemble the regularised Helmholtz G and H influence matrices for the
    ! unbounded exterior domain, including all bounding interfaces.

    SUBROUTINE Get_BRIEFGGHHEX_EM

        INTEGER :: ithprtl, jthprtl, kthprtl, i, j, k
        INTEGER :: NdA, NdB, NdC, NdD, NdE, NdF, id_tp

        INTEGER :: Ttlbim

        DOUBLE PRECISION :: pr0x, pr0y, pr0z, pr0nnx, pr0nny, pr0nnz

        COMPLEX(KIND=KIND(1.0D0)) :: g1, g2, g3, g4, g5, g6
        COMPLEX(KIND=KIND(1.0D0)) :: h1, h2, h3, h4, h5, h6
        COMPLEX(KIND=KIND(1.0D0)) :: gnsbim, hnsbim

        COMPLEX(KIND=KIND(1.0D0)) :: bimGGEXNS,bimHHEXNS, ztp

!*********

        Ttlbim = ttlnmbrnd

!*********
!BRIEF - G & H

        ALLOCATE (bimGGEX(Ttlbim,Ttlbim))
        ALLOCATE (bimHHEX(Ttlbim,Ttlbim))

!$OMP PARALLEL PRIVATE (i,j)
!$OMP DO
        DO j = 1, Ttlbim
            DO i = 1, Ttlbim
                bimGGEX(i,j) = ztpzero
                bimHHEX(i,j) = ztpzero
            END DO
        END DO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE (i,k,ithprtl,jthprtl,id_tp) &
!$OMP & PRIVATE (NdA,NdB,NdC,NdD,NdE,NdF) &
!$OMP & PRIVATE (pr0x,pr0y,pr0z,pr0nnx,pr0nny,pr0nnz) &
!$OMP & PRIVATE (bimGGEXNS,bimHHEXNS,ztp) &
!$OMP & PRIVATE (g1,g2,g3,g4,g5,g6,h1,h2,h3,h4,h5,h6,gnsbim,hnsbim)
!$OMP DO
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

            pr0x = xnd(i)
            pr0y = ynd(i)
            pr0z = znd(i)
            pr0nnx = nnx(i)
            pr0nny = nny(i)
            pr0nnz = nnz(i)

!***
! Assemble the exterior-domain G and H operators.
!---
!There are two cases:
! Case 1: the unbounded exterior domain with surface at infinity
! Case 2: a bounded material domain: corelnkshell(ithprtl)
! For both cases, the surfaces with the same 'corelnkshell(ithprtl)' should be involved

            IF (corelnkshell(ithprtl) == 0) THEN         !Case 1

                IF (NrmlInOut(ithprtl) == 1) THEN
                    bimGGEXNS = ztpzero
                    bimHHEXNS = DCMPLX(4.0d0*pai,0.0d0)     !contribution from infinity
                END IF

                IF (NrmlInOut(ithprtl) ==-1) THEN
                    bimGGEXNS = ztpzero
                    bimHHEXNS =-DCMPLX(4.0d0*pai,0.0d0)     !contribution from infinity
                END IF

                ztp = exk_EM

            END IF

            IF (corelnkshell(ithprtl) > 0) THEN          !Case 2

                bimGGEXNS = ztpzero
                bimHHEXNS = ztpzero

                ztp = ink_EM(corelnkshell(ithprtl))

            END IF

            DO k = 1, ttlnmbrelmnt

                jthprtl = 1
                IF (nmbrprtl > 1 .AND. k > elendID(1)) THEN
                    DO id_tp = 2, nmbrprtl
                        IF (k >= elstaID(id_tp) .AND. k <= elendID(id_tp)) THEN
                            jthprtl = id_tp
                            EXIT
                        END IF
                    END DO
                END IF

                IF (  corelnkshell(jthprtl) == corelnkshell(ithprtl) &
                &.OR. jthprtl == corelnkshell(ithprtl)) THEN

                    IF (MeshType == "L") THEN

                        NdA = elmntlnknd(k,1)
                        NdB = elmntlnknd(k,2)
                        NdC = elmntlnknd(k,3)

                        CALL CalGHLnrBRIEFLnrCOMP( ztp, k, &
                        &       pr0x, pr0y, pr0z, pr0nnx, pr0nny, pr0nnz, &
                        &       g1, g2, g3, h1, h2, h3, gnsbim, hnsbim)

                        bimGGEX(i,NdA)=bimGGEX(i,NdA)+g1
                        bimGGEX(i,NdB)=bimGGEX(i,NdB)+g2
                        bimGGEX(i,NdC)=bimGGEX(i,NdC)+g3

                        bimHHEX(i,NdA)=bimHHEX(i,NdA)+h1
                        bimHHEX(i,NdB)=bimHHEX(i,NdB)+h2
                        bimHHEX(i,NdC)=bimHHEX(i,NdC)+h3

                        bimGGEXNS = bimGGEXNS + gnsbim
                        bimHHEXNS = bimHHEXNS + hnsbim

                    END IF


                    IF (MeshType == "Q") THEN

                        NdA = elmntlnknd(k,1)
                        NdB = elmntlnknd(k,2)
                        NdC = elmntlnknd(k,3)
                        NdD = elmntlnknd(k,4)
                        NdE = elmntlnknd(k,5)
                        NdF = elmntlnknd(k,6)

                        CALL CalGHQdrBRIEFLnrCOMP( ztp, k, &
                        &       pr0x, pr0y, pr0z, pr0nnx, pr0nny, pr0nnz, &
                        &       g1, g2, g3, g4, g5, g6, &
                        &       h1, h2, h3, h4, h5, h6, gnsbim, hnsbim)

                        bimGGEX(i,NdA)=bimGGEX(i,NdA)+g1
                        bimGGEX(i,NdB)=bimGGEX(i,NdB)+g2
                        bimGGEX(i,NdC)=bimGGEX(i,NdC)+g3
                        bimGGEX(i,NdD)=bimGGEX(i,NdD)+g4
                        bimGGEX(i,NdE)=bimGGEX(i,NdE)+g5
                        bimGGEX(i,NdF)=bimGGEX(i,NdF)+g6

                        bimHHEX(i,NdA)=bimHHEX(i,NdA)+h1
                        bimHHEX(i,NdB)=bimHHEX(i,NdB)+h2
                        bimHHEX(i,NdC)=bimHHEX(i,NdC)+h3
                        bimHHEX(i,NdD)=bimHHEX(i,NdD)+h4
                        bimHHEX(i,NdE)=bimHHEX(i,NdE)+h5
                        bimHHEX(i,NdF)=bimHHEX(i,NdF)+h6

                        bimGGEXNS = bimGGEXNS + gnsbim
                        bimHHEXNS = bimHHEXNS + hnsbim


                    END IF

                END IF

            END DO

            bimGGEX(i,i) = bimGGEX(i,i) + bimGGEXNS
            bimHHEX(i,i) = bimHHEX(i,i) + bimHHEXNS

!---
!End calculating external GG and HH
!***

        END DO
!$OMP END DO
!$OMP END PARALLEL

    END SUBROUTINE

    ! Assemble and solve the coupled boundary system using electric-field
    ! scattered traces as the primary representation. Results are written to
    ! the exterior/interior E arrays and their normal-derivative arrays.

    SUBROUTINE SlvPrblm_EM_Esc

        INTEGER :: ithprtl, jthprtl, kthprtl, i, j, k, ii, jj, kk, &
        &          icmpt, id_tp, slf_i, slf_j, slf_k, GLQi, icnt

        INTEGER :: NdA, NdB, NdC, NdD, NdE, NdF, ndst, nded, elst, eled

        INTEGER :: Ttlbim, MatSize
        INTEGER, ALLOCATABLE, DIMENSION (:) :: mtidst

        COMPLEX(KIND=KIND(1.0D0)) :: mdm_oi, mdm_io

        DOUBLE PRECISION :: tp, tp1, tp2, tp3, tp4, tp5, tp6, tp7, tp8, tp9
        COMPLEX(KIND=KIND(1.0D0)) :: ztp,ztp1,ztp2,ztp3,ztp4,ztp5,ztp6,ztp7,ztp8,ztp9

        ALLOCATE (DivFrAA(ttlnmbrnd,6*ttlddtnd))
        ALLOCATE (DivFrBB(ttlnmbrnd))

!$OMP PARALLEL PRIVATE (i,j)
!$OMP DO
        DO i = 1, ttlnmbrnd
            DO j = 1, 6*ttlddtnd
                DivFrAA(i,j) = ztpzero
            END DO
            DivFrBB(i) = ztpzero
        END DO
!$OMP END DO
!$OMP END PARALLEL

        ALLOCATE (TanAAt1(ttlnmbrnd,6*ttlddtnd))
        ALLOCATE (TanAAt2(ttlnmbrnd,6*ttlddtnd))
        ALLOCATE (TanBBt1(ttlnmbrnd))
        ALLOCATE (TanBBt2(ttlnmbrnd))

!$OMP PARALLEL PRIVATE (i,j)
!$OMP DO
        DO i = 1, ttlnmbrnd
            DO j = 1, 6*ttlddtnd
                TanAAt1(i,j) = ztpzero
                TanAAt2(i,j) = ztpzero
            END DO
            TanBBt1(i) = ztpzero
            TanBBt2(i) = ztpzero
        END DO
!$OMP END DO
!$OMP END PARALLEL

        DO ithprtl = 1, nmbrprtl
            IF (BCType_EM(ithprtl) == 'PEC') THEN
                CALL Get_PEC_DivFr(ithprtl, 'E')
                CALL Get_PEC_TanCn(ithprtl, 'E')
            END IF
            IF (BCType_EM(ithprtl) == '2SD') THEN
                CALL Get_DEL_DivFrsc(ithprtl, 'E')
                CALL Get_DEL_TanCnsc(ithprtl, 'E')
            END IF
        END DO


        Ttlbim = ttlnmbrnd
        ALLOCATE (mtidst(ttlnmbrnd))

        MatSize = 0
        DO ithprtl = 1, nmbrprtl
            IF (BCType_EM(ithprtl) == 'PEC') THEN
                id_tp = 1
                DO i = ndstaID(ithprtl), ndendID(ithprtl)
                    mtidst(i) = MatSize + id_tp
                    id_tp = id_tp + 3
                END DO
                MatSize = MatSize + 3*nmbrnd(ithprtl)
            END IF
            IF (BCType_EM(ithprtl) == '2SD') THEN
                id_tp = 1
                DO i = ndstaID(ithprtl), ndendID(ithprtl)
                    mtidst(i) = MatSize + id_tp
                    id_tp = id_tp + 6
                END DO
                MatSize = MatSize + 6*nmbrnd(ithprtl)
            END IF
        END DO

        ALLOCATE(bimAA(MatSize,MatSize))
        ALLOCATE(bimBB(MatSize))
        ALLOCATE(bimX(MatSize))

!$OMP PARALLEL PRIVATE (i,j)
!$OMP DO
        DO i = 1, MatSize
            DO j = 1, MatSize
                bimAA(i,j) = ztpzero
            END DO
            bimBB(i) = ztpzero
            bimX(i) = ztpzero
        END DO
!$OMP END DO
!$OMP END PARALLEL

        CALL Get_BRIEFGGHHEX_EM
        CALL Get_BRIEFGGHHIN_EM

!$OMP PARALLEL PRIVATE (i,j,ii,jj,kk,NdA) &
!$OMP & PRIVATE (ithprtl,jthprtl,id_tp) &
!$OMP & PRIVATE (tp,tp1,tp2,tp3,ztp,ztp1,ztp2,ztp3,ztp4,ztp5,mdm_oi,mdm_io)
!$OMP DO
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

!***
! Exterior-field equations.
            IF (     BCType_EM(ithprtl) == 'PEC' &
            &   .OR. BCType_EM(ithprtl) == '2SD') THEN

!---
! Helmholtz equations for E, rows 1--3 of the assembled system.
!Part 1: surfaces share the same 'corelnkshell'
!!!NOTE: 'E' is 'exE2'

                DO jthprtl = 1, nmbrprtl

                    id_tp = corelnkshell(jthprtl)

                    IF (corelnkshell(jthprtl) == corelnkshell(ithprtl)) THEN

                        DO j = ndstaID(jthprtl), ndendID(jthprtl)

                            IF (BCType_EM(jthprtl) == 'PEC') THEN

                                !exE2nn
                                bimAA(mtidst(i)  ,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)  ,mtidst(j)  ) + bimHHEX(i,j)*nnx(j)
                                bimAA(mtidst(i)+1,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)+1,mtidst(j)  ) + bimHHEX(i,j)*nny(j)
                                bimAA(mtidst(i)+2,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)+2,mtidst(j)  ) + bimHHEX(i,j)*nnz(j)

                                !t1•∂(exE2)/∂n
                                bimAA(mtidst(i)  ,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)  ,mtidst(j)+1) - bimGGEX(i,j)*t1x(j)
                                bimAA(mtidst(i)+1,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)+1,mtidst(j)+1) - bimGGEX(i,j)*t1y(j)
                                bimAA(mtidst(i)+2,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)+2,mtidst(j)+1) - bimGGEX(i,j)*t1z(j)

                                !t2•∂(exE2)/∂n
                                bimAA(mtidst(i)  ,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)  ,mtidst(j)+2) - bimGGEX(i,j)*t2x(j)
                                bimAA(mtidst(i)+1,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)+1,mtidst(j)+2) - bimGGEX(i,j)*t2y(j)
                                bimAA(mtidst(i)+2,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)+2,mtidst(j)+2) - bimGGEX(i,j)*t2z(j)

                                !n•∂(exE2)/∂n by BC: [exE2nn]
                                DO kk = 1, mxnmbrndlnknd2ndslf
                                    IF (ndlnknd2ndslf(j,kk) /= 0) THEN

                                        NdA = ndlnknd2ndslf(j,kk)
                                        ii = 6*(kk-1)

                                        bimAA(mtidst(i)  ,mtidst(NdA)  ) &
                                    &=  bimAA(mtidst(i)  ,mtidst(NdA)  ) &
                                    &       -bimGGEX(i,j)*nnx(j)*DivFrAA(j,ii+3)
                                        bimAA(mtidst(i)+1,mtidst(NdA)  ) &
                                    &=  bimAA(mtidst(i)+1,mtidst(NdA)  ) &
                                    &       -bimGGEX(i,j)*nny(j)*DivFrAA(j,ii+3)
                                        bimAA(mtidst(i)+2,mtidst(NdA)  ) &
                                    &=  bimAA(mtidst(i)+2,mtidst(NdA)  ) &
                                    &       -bimGGEX(i,j)*nnz(j)*DivFrAA(j,ii+3)

                                    END IF
                                END DO

                                !n•∂(exE2)/∂n by BC: [exE1t1, exE1t2]
                                bimBB(mtidst(i)  ) = bimBB(mtidst(i)  ) &
                                &                   +bimGGEX(i,j)*nnx(j)*DivFrBB(j)
                                bimBB(mtidst(i)+1) = bimBB(mtidst(i)+1) &
                                &                   +bimGGEX(i,j)*nny(j)*DivFrBB(j)
                                bimBB(mtidst(i)+2) = bimBB(mtidst(i)+2) &
                                &                   +bimGGEX(i,j)*nnz(j)*DivFrBB(j)

                                ztp1 = exE1x_EM(j)*t1x(j)&
                                &     +exE1y_EM(j)*t1y(j)&
                                &     +exE1z_EM(j)*t1z(j)
                                ztp2 = exE1x_EM(j)*t2x(j)&
                                &     +exE1y_EM(j)*t2y(j)&
                                &     +exE1z_EM(j)*t2z(j)

                                bimBB(mtidst(i)  ) = bimBB(mtidst(i)  ) &
                                &       +bimHHEX(i,j) * (ztp1*t1x(j) + ztp2*t2x(j))

                                bimBB(mtidst(i)+1) = bimBB(mtidst(i)+1) &
                                &       +bimHHEX(i,j) * (ztp1*t1y(j) + ztp2*t2y(j))

                                bimBB(mtidst(i)+2) = bimBB(mtidst(i)+2) &
                                &       +bimHHEX(i,j) * (ztp1*t1z(j) + ztp2*t2z(j))

                            END IF

                            IF (BCType_EM(jthprtl) == '2SD') THEN

                                !exE2t1
                                bimAA(mtidst(i)  ,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)  ,mtidst(j)  ) + bimHHEX(i,j)*t1x(j)
                                bimAA(mtidst(i)+1,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)+1,mtidst(j)  ) + bimHHEX(i,j)*t1y(j)
                                bimAA(mtidst(i)+2,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)+2,mtidst(j)  ) + bimHHEX(i,j)*t1z(j)

                                !exE2t2
                                bimAA(mtidst(i)  ,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)  ,mtidst(j)+1) + bimHHEX(i,j)*t2x(j)
                                bimAA(mtidst(i)+1,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)+1,mtidst(j)+1) + bimHHEX(i,j)*t2y(j)
                                bimAA(mtidst(i)+2,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)+2,mtidst(j)+1) + bimHHEX(i,j)*t2z(j)

                                !exE2nn
                                bimAA(mtidst(i)  ,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)  ,mtidst(j)+2) + bimHHEX(i,j)*nnx(j)
                                bimAA(mtidst(i)+1,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)+1,mtidst(j)+2) + bimHHEX(i,j)*nny(j)
                                bimAA(mtidst(i)+2,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)+2,mtidst(j)+2) + bimHHEX(i,j)*nnz(j)

                                !t1•∂(exE2)/∂n
                                bimAA(mtidst(i)  ,mtidst(j)+3) &
                            &=  bimAA(mtidst(i)  ,mtidst(j)+3) - bimGGEX(i,j)*t1x(j)
                                bimAA(mtidst(i)+1,mtidst(j)+3) &
                            &=  bimAA(mtidst(i)+1,mtidst(j)+3) - bimGGEX(i,j)*t1y(j)
                                bimAA(mtidst(i)+2,mtidst(j)+3) &
                            &=  bimAA(mtidst(i)+2,mtidst(j)+3) - bimGGEX(i,j)*t1z(j)

                                !t2•∂(exE2)/∂n
                                bimAA(mtidst(i)  ,mtidst(j)+4) &
                            &=  bimAA(mtidst(i)  ,mtidst(j)+4) - bimGGEX(i,j)*t2x(j)
                                bimAA(mtidst(i)+1,mtidst(j)+4) &
                            &=  bimAA(mtidst(i)+1,mtidst(j)+4) - bimGGEX(i,j)*t2y(j)
                                bimAA(mtidst(i)+2,mtidst(j)+4) &
                            &=  bimAA(mtidst(i)+2,mtidst(j)+4) - bimGGEX(i,j)*t2z(j)

                                !n•∂(exE2)/∂n
                                bimAA(mtidst(i)  ,mtidst(j)+5) &
                            &=  bimAA(mtidst(i)  ,mtidst(j)+5) - bimGGEX(i,j)*nnx(j)
                                bimAA(mtidst(i)+1,mtidst(j)+5) &
                            &=  bimAA(mtidst(i)+1,mtidst(j)+5) - bimGGEX(i,j)*nny(j)
                                bimAA(mtidst(i)+2,mtidst(j)+5) &
                            &=  bimAA(mtidst(i)+2,mtidst(j)+5) - bimGGEX(i,j)*nnz(j)

                            END IF

                        END DO

                    END IF

                END DO

!---

!---
! Continue the E rows 1--3 for bounded exterior-domain surfaces.
!!!NOTE: on this bounded surface, 'E' is 'inE2'
!!!NOTE: ROW 4 to 6 with GGEX and HHEX in the big Matrix (Table 1) in Die_MS II
!!!NOTE: Eq IDs are still 1 to 3, such as mtidst(i), mtidst(i)+1, mtidst(i)+2

                IF (corelnkshell(ithprtl) > 0) THEN

                    jthprtl = corelnkshell(ithprtl)

                    id_tp = corelnkshell(jthprtl)
                    IF (id_tp == 0) THEN
                        mdm_oi = exeps_EM/ineps_EM(jthprtl)
                    ELSE
                        mdm_oi = ineps_EM(id_tp)/ineps_EM(jthprtl)
                    END IF

                    DO j = ndstaID(jthprtl), ndendID(jthprtl)

                        IF (BCType_EM(jthprtl) == '2SD') THEN

                            bimAA(mtidst(i)  ,mtidst(j)  ) &
                        &=  bimAA(mtidst(i)  ,mtidst(j)  ) + bimHHEX(i,j)*t1x(j)
                            bimAA(mtidst(i)+1,mtidst(j)  ) &
                        &=  bimAA(mtidst(i)+1,mtidst(j)  ) + bimHHEX(i,j)*t1y(j)
                            bimAA(mtidst(i)+2,mtidst(j)  ) &
                        &=  bimAA(mtidst(i)+2,mtidst(j)  ) + bimHHEX(i,j)*t1z(j)

                            bimAA(mtidst(i)  ,mtidst(j)+1) &
                        &=  bimAA(mtidst(i)  ,mtidst(j)+1) + bimHHEX(i,j)*t2x(j)
                            bimAA(mtidst(i)+1,mtidst(j)+1) &
                        &=  bimAA(mtidst(i)+1,mtidst(j)+1) + bimHHEX(i,j)*t2y(j)
                            bimAA(mtidst(i)+2,mtidst(j)+1) &
                        &=  bimAA(mtidst(i)+2,mtidst(j)+1) + bimHHEX(i,j)*t2z(j)

                            bimAA(mtidst(i)  ,mtidst(j)+2) &
                        &=  bimAA(mtidst(i)  ,mtidst(j)+2) + bimHHEX(i,j)*nnx(j)*mdm_oi
                            bimAA(mtidst(i)+1,mtidst(j)+2) &
                        &=  bimAA(mtidst(i)+1,mtidst(j)+2) + bimHHEX(i,j)*nny(j)*mdm_oi
                            bimAA(mtidst(i)+2,mtidst(j)+2) &
                        &=  bimAA(mtidst(i)+2,mtidst(j)+2) + bimHHEX(i,j)*nnz(j)*mdm_oi

![t1•∂(inE2)/∂n,t2•∂(inE2)/∂n,n•∂(inE2)/∂n] by BC:
![exE2t1,exE1t2,exE2nn,t1•∂(exE2)/∂n,t2•∂(exE2)/∂n]
                            DO kk = 1, mxnmbrndlnknd2ndslf
                                IF (ndlnknd2ndslf(j,kk) /= 0) THEN

                                    NdA = ndlnknd2ndslf(j,kk)
                                    ii = 6*(kk-1)

                                    DO id_tp = 1, 6

                                        bimAA(mtidst(i)  ,mtidst(NdA)+id_tp-1) &
                                    &=  bimAA(mtidst(i)  ,mtidst(NdA)+id_tp-1) &
                                    &       -bimGGEX(i,j)*( t1x(j)*TanAAt1(j,ii+id_tp) &
                                    &                      +t2x(j)*TanAAt2(j,ii+id_tp) &
                                    &                      +nnx(j)*DivFrAA(j,ii+id_tp) )

                                        bimAA(mtidst(i)+1,mtidst(NdA)+id_tp-1) &
                                    &=  bimAA(mtidst(i)+1,mtidst(NdA)+id_tp-1) &
                                    &       -bimGGEX(i,j)*( t1y(j)*TanAAt1(j,ii+id_tp) &
                                    &                      +t2y(j)*TanAAt2(j,ii+id_tp) &
                                    &                      +nny(j)*DivFrAA(j,ii+id_tp) )

                                        bimAA(mtidst(i)+2,mtidst(NdA)+id_tp-1) &
                                    &=  bimAA(mtidst(i)+2,mtidst(NdA)+id_tp-1) &
                                    &       -bimGGEX(i,j)*( t1z(j)*TanAAt1(j,ii+id_tp) &
                                    &                      +t2z(j)*TanAAt2(j,ii+id_tp) &
                                    &                      +nnz(j)*DivFrAA(j,ii+id_tp) )

                                    END DO

                                END IF
                            END DO

                            ztp1 = t1x(j)*(exE1x_EM(j)-inE1x_EM(j)) &
                            &     +t1y(j)*(exE1y_EM(j)-inE1y_EM(j)) &
                            &     +t1z(j)*(exE1z_EM(j)-inE1z_EM(j))
                            ztp2 = t2x(j)*(exE1x_EM(j)-inE1x_EM(j)) &
                            &     +t2y(j)*(exE1y_EM(j)-inE1y_EM(j)) &
                            &     +t2z(j)*(exE1z_EM(j)-inE1z_EM(j))
                            ztp3 = nnx(j)* exE1x_EM(j) &
                            &     +nny(j)* exE1y_EM(j) &
                            &     +nnz(j)* exE1z_EM(j)
                            ztp3 = ztp3 * mdm_oi
                            ztp  = nnx(j)* inE1x_EM(j) &
                            &     +nny(j)* inE1y_EM(j) &
                            &     +nnz(j)* inE1z_EM(j)

                            bimBB(mtidst(i)  ) = bimBB(mtidst(i)  ) &
                            &       -bimHHEX(i,j) * ( t1x(j)*ztp1+t2x(j)*ztp2 &
                            &                        +nnx(j)*ztp3-nnx(j)*ztp  )
                            bimBB(mtidst(i)+1) = bimBB(mtidst(i)+1) &
                            &       -bimHHEX(i,j) * ( t1y(j)*ztp1+t2y(j)*ztp2 &
                            &                        +nny(j)*ztp3-nny(j)*ztp  )
                            bimBB(mtidst(i)+2) = bimBB(mtidst(i)+2) &
                            &       -bimHHEX(i,j) * ( t1z(j)*ztp1+t2z(j)*ztp2 &
                            &                        +nnz(j)*ztp3-nnz(j)*ztp  )

![t1•∂(inE2)/∂n,t2•∂(inE2)/∂n,n•∂(inE2)/∂n] by BC:
![exE2t1,exE1t2,exE2nn,t1•∂(exE2)/∂n,t2•∂(exE2)/∂n]
                            ztp1 = t1x(j)*TanBBt1(j) &
                            &     +t2x(j)*TanBBt2(j) &
                            &     +nnx(j)*DivFrBB(j)
                            ztp2 = t1y(j)*TanBBt1(j) &
                            &     +t2y(j)*TanBBt2(j) &
                            &     +nny(j)*DivFrBB(j)
                            ztp3 = t1z(j)*TanBBt1(j) &
                            &     +t2z(j)*TanBBt2(j) &
                            &     +nnz(j)*DivFrBB(j)
                            bimBB(mtidst(i)  ) = bimBB(mtidst(i)  ) + bimGGEX(i,j) * ztp1
                            bimBB(mtidst(i)+1) = bimBB(mtidst(i)+1) + bimGGEX(i,j) * ztp2
                            bimBB(mtidst(i)+2) = bimBB(mtidst(i)+2) + bimGGEX(i,j) * ztp3

                        END IF

                    END DO

                END IF

!---

            END IF

!End External field
!***

!***
! Interior-field equations.
            IF (BCType_EM(ithprtl) == '2SD') THEN

!---
! Helmholtz equations for E, rows 4--6 of the assembled system.

                DO jthprtl = 1, nmbrprtl

                    IF (jthprtl == ithprtl) THEN

                        id_tp = corelnkshell(jthprtl)
                        IF (id_tp == 0) THEN
                            mdm_oi = exeps_EM/ineps_EM(jthprtl)
                        ELSE
                            mdm_oi = ineps_EM(id_tp)/ineps_EM(jthprtl)
                        END IF

                        DO j = ndstaID(jthprtl), ndendID(jthprtl)

                            bimAA(mtidst(i)+3,mtidst(j)  ) &
                        &=  bimAA(mtidst(i)+3,mtidst(j)  ) + bimHHIN(i,j)*t1x(j)
                            bimAA(mtidst(i)+4,mtidst(j)  ) &
                        &=  bimAA(mtidst(i)+4,mtidst(j)  ) + bimHHIN(i,j)*t1y(j)
                            bimAA(mtidst(i)+5,mtidst(j)  ) &
                        &=  bimAA(mtidst(i)+5,mtidst(j)  ) + bimHHIN(i,j)*t1z(j)

                            bimAA(mtidst(i)+3,mtidst(j)+1) &
                        &=  bimAA(mtidst(i)+3,mtidst(j)+1) + bimHHIN(i,j)*t2x(j)
                            bimAA(mtidst(i)+4,mtidst(j)+1) &
                        &=  bimAA(mtidst(i)+4,mtidst(j)+1) + bimHHIN(i,j)*t2y(j)
                            bimAA(mtidst(i)+5,mtidst(j)+1) &
                        &=  bimAA(mtidst(i)+5,mtidst(j)+1) + bimHHIN(i,j)*t2z(j)

                            bimAA(mtidst(i)+3,mtidst(j)+2) &
                        &=  bimAA(mtidst(i)+3,mtidst(j)+2) + bimHHIN(i,j)*nnx(j)*mdm_oi
                            bimAA(mtidst(i)+4,mtidst(j)+2) &
                        &=  bimAA(mtidst(i)+4,mtidst(j)+2) + bimHHIN(i,j)*nny(j)*mdm_oi
                            bimAA(mtidst(i)+5,mtidst(j)+2) &
                        &=  bimAA(mtidst(i)+5,mtidst(j)+2) + bimHHIN(i,j)*nnz(j)*mdm_oi

![t1•∂(inE2)/∂n,t2•∂(inE2)/∂n,n•∂(inE2)/∂n] by BC:
![exE2t1,exE1t2,exE2nn,t1•∂(exE2)/∂n,t2•∂(exE2)/∂n]
                            DO kk = 1, mxnmbrndlnknd2ndslf
                                IF (ndlnknd2ndslf(j,kk) /= 0) THEN

                                    NdA = ndlnknd2ndslf(j,kk)
                                    ii = 6*(kk-1)

                                    DO id_tp = 1, 6

                                        bimAA(mtidst(i)+3,mtidst(NdA)+id_tp-1) &
                                    &=  bimAA(mtidst(i)+3,mtidst(NdA)+id_tp-1) &
                                    &       -bimGGIN(i,j)*( t1x(j)*TanAAt1(j,ii+id_tp) &
                                    &                      +t2x(j)*TanAAt2(j,ii+id_tp) &
                                    &                      +nnx(j)*DivFrAA(j,ii+id_tp) )

                                        bimAA(mtidst(i)+4,mtidst(NdA)+id_tp-1) &
                                    &=  bimAA(mtidst(i)+4,mtidst(NdA)+id_tp-1) &
                                    &       -bimGGIN(i,j)*( t1y(j)*TanAAt1(j,ii+id_tp) &
                                    &                      +t2y(j)*TanAAt2(j,ii+id_tp) &
                                    &                      +nny(j)*DivFrAA(j,ii+id_tp) )

                                        bimAA(mtidst(i)+5,mtidst(NdA)+id_tp-1) &
                                    &=  bimAA(mtidst(i)+5,mtidst(NdA)+id_tp-1) &
                                    &       -bimGGIN(i,j)*( t1z(j)*TanAAt1(j,ii+id_tp) &
                                    &                      +t2z(j)*TanAAt2(j,ii+id_tp) &
                                    &                      +nnz(j)*DivFrAA(j,ii+id_tp) )

                                    END DO

                                END IF
                            END DO

                            ztp1 = t1x(j)*(exE1x_EM(j)-inE1x_EM(j)) &
                            &     +t1y(j)*(exE1y_EM(j)-inE1y_EM(j)) &
                            &     +t1z(j)*(exE1z_EM(j)-inE1z_EM(j))
                            ztp2 = t2x(j)*(exE1x_EM(j)-inE1x_EM(j)) &
                            &     +t2y(j)*(exE1y_EM(j)-inE1y_EM(j)) &
                            &     +t2z(j)*(exE1z_EM(j)-inE1z_EM(j))
                            ztp3 = nnx(j)* exE1x_EM(j) &
                            &     +nny(j)* exE1y_EM(j) &
                            &     +nnz(j)* exE1z_EM(j)
                            ztp3 = ztp3 * mdm_oi
                            ztp  = nnx(j)* inE1x_EM(j) &
                            &     +nny(j)* inE1y_EM(j) &
                            &     +nnz(j)* inE1z_EM(j)

                            bimBB(mtidst(i)+3) = bimBB(mtidst(i)+3) &
                            &       -bimHHIN(i,j) * ( t1x(j)*ztp1+t2x(j)*ztp2 &
                            &                        +nnx(j)*ztp3-nnx(j)*ztp  )
                            bimBB(mtidst(i)+4) = bimBB(mtidst(i)+4) &
                            &       -bimHHIN(i,j) * ( t1y(j)*ztp1+t2y(j)*ztp2 &
                            &                        +nny(j)*ztp3-nny(j)*ztp  )
                            bimBB(mtidst(i)+5) = bimBB(mtidst(i)+5) &
                            &       -bimHHIN(i,j) * ( t1z(j)*ztp1+t2z(j)*ztp2 &
                            &                        +nnz(j)*ztp3-nnz(j)*ztp  )

![t1•∂(inE2)/∂n,t2•∂(inE2)/∂n,n•∂(inE2)/∂n] by BC:
![exE2t1,exE1t2,exE2nn,t1•∂(exE2)/∂n,t2•∂(exE2)/∂n]
                            ztp1 = t1x(j)*TanBBt1(j) &
                            &     +t2x(j)*TanBBt2(j) &
                            &     +nnx(j)*DivFrBB(j)
                            ztp2 = t1y(j)*TanBBt1(j) &
                            &     +t2y(j)*TanBBt2(j) &
                            &     +nny(j)*DivFrBB(j)
                            ztp3 = t1z(j)*TanBBt1(j) &
                            &     +t2z(j)*TanBBt2(j) &
                            &     +nnz(j)*DivFrBB(j)
                            bimBB(mtidst(i)+3) = bimBB(mtidst(i)+3) + bimGGIN(i,j) * ztp1
                            bimBB(mtidst(i)+4) = bimBB(mtidst(i)+4) + bimGGIN(i,j) * ztp2
                            bimBB(mtidst(i)+5) = bimBB(mtidst(i)+5) + bimGGIN(i,j) * ztp3
                        END DO

                    END IF

!---

!---
!Part2: the surfaces that enclosed by the surface where x0 locates
!(just one level down needed)
!!NOTE: 'E' is 'exE2'
!!NOTE: ROW 1 to 3 with GGIN and HHIN in the big Matrix (Table 1) in Die_MS II
!!NOTE: Eq IDs are still 4 to 6, such as mtidst(i)+3, mtidst(i)+4, mtidst(i)+5

                    IF (corelnkshell(jthprtl) == ithprtl) THEN

                        DO j = ndstaID(jthprtl), ndendID(jthprtl)

                            IF (BCType_EM(jthprtl) == 'PEC') THEN

                                !exE2nn
                                bimAA(mtidst(i)+3,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)+3,mtidst(j)  ) + bimHHIN(i,j)*nnx(j)
                                bimAA(mtidst(i)+4,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)+4,mtidst(j)  ) + bimHHIN(i,j)*nny(j)
                                bimAA(mtidst(i)+5,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)+5,mtidst(j)  ) + bimHHIN(i,j)*nnz(j)

                                !t1•∂(exE2)/∂n
                                bimAA(mtidst(i)+3,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)+3,mtidst(j)+1) - bimGGIN(i,j)*t1x(j)
                                bimAA(mtidst(i)+4,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)+4,mtidst(j)+1) - bimGGIN(i,j)*t1y(j)
                                bimAA(mtidst(i)+5,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)+5,mtidst(j)+1) - bimGGIN(i,j)*t1z(j)

                                !t2•∂(exE2)/∂n
                                bimAA(mtidst(i)+3,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)+3,mtidst(j)+2) - bimGGIN(i,j)*t2x(j)
                                bimAA(mtidst(i)+4,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)+4,mtidst(j)+2) - bimGGIN(i,j)*t2y(j)
                                bimAA(mtidst(i)+5,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)+5,mtidst(j)+2) - bimGGIN(i,j)*t2z(j)

                                !n•∂(exE2)/∂n by BC: [exE2nn]
                                DO kk = 1, mxnmbrndlnknd2ndslf
                                    IF (ndlnknd2ndslf(j,kk) /= 0) THEN

                                        NdA = ndlnknd2ndslf(j,kk)
                                        ii = 6*(kk-1)

                                        bimAA(mtidst(i)+3,mtidst(NdA)  ) &
                                    &=  bimAA(mtidst(i)+3,mtidst(NdA)  ) &
                                    &       -bimGGIN(i,j)*nnx(j)*DivFrAA(j,ii+3)
                                        bimAA(mtidst(i)+4,mtidst(NdA)  ) &
                                    &=  bimAA(mtidst(i)+4,mtidst(NdA)  ) &
                                    &       -bimGGIN(i,j)*nny(j)*DivFrAA(j,ii+3)
                                        bimAA(mtidst(i)+5,mtidst(NdA)  ) &
                                    &=  bimAA(mtidst(i)+5,mtidst(NdA)  ) &
                                    &       -bimGGIN(i,j)*nnz(j)*DivFrAA(j,ii+3)

                                    END IF
                                END DO

                                !n•∂(exE2)/∂n by BC: [exE1t1, exE1t2]
                                bimBB(mtidst(i)+3) = bimBB(mtidst(i)+3) &
                                &                   +bimGGIN(i,j)*nnx(j)*DivFrBB(j)
                                bimBB(mtidst(i)+4) = bimBB(mtidst(i)+4) &
                                &                   +bimGGIN(i,j)*nny(j)*DivFrBB(j)
                                bimBB(mtidst(i)+5) = bimBB(mtidst(i)+5) &
                                &                   +bimGGIN(i,j)*nnz(j)*DivFrBB(j)

                                ztp1 = exE1x_EM(j)*t1x(j)&
                                &     +exE1y_EM(j)*t1y(j)&
                                &     +exE1z_EM(j)*t1z(j)
                                ztp2 = exE1x_EM(j)*t2x(j)&
                                &     +exE1y_EM(j)*t2y(j)&
                                &     +exE1z_EM(j)*t2z(j)

                                bimBB(mtidst(i)+3) = bimBB(mtidst(i)+3) &
                                &       +bimHHIN(i,j) * (ztp1*t1x(j) + ztp2*t2x(j))

                                bimBB(mtidst(i)+4) = bimBB(mtidst(i)+4) &
                                &       +bimHHIN(i,j) * (ztp1*t1y(j) + ztp2*t2y(j))

                                bimBB(mtidst(i)+5) = bimBB(mtidst(i)+5) &
                                &       +bimHHIN(i,j) * (ztp1*t1z(j) + ztp2*t2z(j))

                            END IF

                            IF (BCType_EM(jthprtl) == '2SD') THEN

                                !exE2t1
                                bimAA(mtidst(i)+3,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)+3,mtidst(j)  ) + bimHHIN(i,j)*t1x(j)
                                bimAA(mtidst(i)+4,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)+4,mtidst(j)  ) + bimHHIN(i,j)*t1y(j)
                                bimAA(mtidst(i)+5,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)+5,mtidst(j)  ) + bimHHIN(i,j)*t1z(j)

                                !exE2t2
                                bimAA(mtidst(i)+3,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)+3,mtidst(j)+1) + bimHHIN(i,j)*t2x(j)
                                bimAA(mtidst(i)+4,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)+4,mtidst(j)+1) + bimHHIN(i,j)*t2y(j)
                                bimAA(mtidst(i)+5,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)+5,mtidst(j)+1) + bimHHIN(i,j)*t2z(j)

                                !exE2nn
                                bimAA(mtidst(i)+3,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)+3,mtidst(j)+2) + bimHHIN(i,j)*nnx(j)
                                bimAA(mtidst(i)+4,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)+4,mtidst(j)+2) + bimHHIN(i,j)*nny(j)
                                bimAA(mtidst(i)+5,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)+5,mtidst(j)+2) + bimHHIN(i,j)*nnz(j)

                                !t1•∂(exE2)/∂n
                                bimAA(mtidst(i)+3,mtidst(j)+3) &
                            &=  bimAA(mtidst(i)+3,mtidst(j)+3) - bimGGIN(i,j)*t1x(j)
                                bimAA(mtidst(i)+4,mtidst(j)+3) &
                            &=  bimAA(mtidst(i)+4,mtidst(j)+3) - bimGGIN(i,j)*t1y(j)
                                bimAA(mtidst(i)+5,mtidst(j)+3) &
                            &=  bimAA(mtidst(i)+5,mtidst(j)+3) - bimGGIN(i,j)*t1z(j)

                                !t2•∂(exE2)/∂n
                                bimAA(mtidst(i)+3,mtidst(j)+4) &
                            &=  bimAA(mtidst(i)+3,mtidst(j)+4) - bimGGIN(i,j)*t2x(j)
                                bimAA(mtidst(i)+4,mtidst(j)+4) &
                            &=  bimAA(mtidst(i)+4,mtidst(j)+4) - bimGGIN(i,j)*t2y(j)
                                bimAA(mtidst(i)+5,mtidst(j)+4) &
                            &=  bimAA(mtidst(i)+5,mtidst(j)+4) - bimGGIN(i,j)*t2z(j)

                                !n•∂(exE2)/∂n
                                bimAA(mtidst(i)+3,mtidst(j)+5) &
                            &=  bimAA(mtidst(i)+3,mtidst(j)+5) - bimGGIN(i,j)*nnx(j)
                                bimAA(mtidst(i)+4,mtidst(j)+5) &
                            &=  bimAA(mtidst(i)+4,mtidst(j)+5) - bimGGIN(i,j)*nny(j)
                                bimAA(mtidst(i)+5,mtidst(j)+5) &
                            &=  bimAA(mtidst(i)+5,mtidst(j)+5) - bimGGIN(i,j)*nnz(j)

                            END IF

                        END DO

                    END IF

                END DO

!---

            END IF

!End Internal field
!***

        END DO
!$OMP END DO
!$OMP END PARALLEL

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!second argument = 0 -> LU decomposition
        CALL Slv_LinearMatrixCOMP_LU(MatSize, 0)

!$OMP PARALLEL PRIVATE (i,j,ii,jj,kk,NdA) &
!$OMP & PRIVATE (ithprtl,jthprtl,id_tp) &
!$OMP & PRIVATE (tp,tp1,tp2,tp3,ztp,ztp1,ztp2,ztp3,ztp4,ztp5,mdm_oi,mdm_io)
!$OMP DO

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

            id_tp = corelnkshell(ithprtl)

            IF (BCType_EM(ithprtl) == 'PEC') THEN

                ztp1 = exE1x_EM(i)*t1x(i)&
                &     +exE1y_EM(i)*t1y(i)&
                &     +exE1z_EM(i)*t1z(i)
                ztp1 =-ztp1
                ztp2 = exE1x_EM(i)*t2x(i)&
                &     +exE1y_EM(i)*t2y(i)&
                &     +exE1z_EM(i)*t2z(i)
                ztp2 =-ztp2

                exE2x_EM(i) = bimX(mtidst(i)  )*nnx(i) &
                &            +             ztp1*t1x(i) &
                &            +             ztp2*t2x(i)
                exE2y_EM(i) = bimX(mtidst(i)  )*nny(i) &
                &            +             ztp1*t1y(i) &
                &            +             ztp2*t2y(i)
                exE2z_EM(i) = bimX(mtidst(i)  )*nnz(i) &
                &            +             ztp1*t1z(i) &
                &            +             ztp2*t2z(i)

                exE3x_EM(i) = exE1x_EM(i) + exE2x_EM(i)
                exE3y_EM(i) = exE1y_EM(i) + exE2y_EM(i)
                exE3z_EM(i) = exE1z_EM(i) + exE2z_EM(i)

                ztp3 = ztpzero
                DO kk = 1, mxnmbrndlnknd2ndslf
                    IF (ndlnknd2ndslf(i,kk) /= 0) THEN

                        NdA = ndlnknd2ndslf(i,kk)
                        ii = bcscl*(kk-1)

                        ztp3 = ztp3 + bimX(mtidst(NdA)  ) * DivFrAA(i,ii+3)

                    END IF
                END DO

                ztp3 = ztp3 + DivFrBB(i)

                exE2xdnn_EM(i) = ztp3              *nnx(i) &
                &                +bimX(mtidst(i)+1)*t1x(i) &
                &                +bimX(mtidst(i)+2)*t2x(i)
                exE2ydnn_EM(i) = ztp3              *nny(i) &
                &                +bimX(mtidst(i)+1)*t1y(i) &
                &                +bimX(mtidst(i)+2)*t2y(i)
                exE2zdnn_EM(i) = ztp3              *nnz(i) &
                &                +bimX(mtidst(i)+1)*t1z(i) &
                &                +bimX(mtidst(i)+2)*t2z(i)

                exE3xdnn_EM(i) = exE1xdnn_EM(i) + exE2xdnn_EM(i)
                exE3ydnn_EM(i) = exE1ydnn_EM(i) + exE2ydnn_EM(i)
                exE3zdnn_EM(i) = exE1zdnn_EM(i) + exE2zdnn_EM(i)

            END IF

            IF (BCType_EM(ithprtl) == '2SD') THEN

                id_tp = corelnkshell(ithprtl)
                IF (id_tp == 0) THEN
                    mdm_oi = exeps_EM/ineps_EM(ithprtl)
                ELSE
                    mdm_oi = ineps_EM(id_tp)/ineps_EM(ithprtl)
                END IF

                exE2x_EM(i) = bimX(mtidst(i)  )*t1x(i) &
                &            +bimX(mtidst(i)+1)*t2x(i) &
                &            +bimX(mtidst(i)+2)*nnx(i)
                exE2y_EM(i) = bimX(mtidst(i)  )*t1y(i) &
                &            +bimX(mtidst(i)+1)*t2y(i) &
                &            +bimX(mtidst(i)+2)*nny(i)
                exE2z_EM(i) = bimX(mtidst(i)  )*t1z(i) &
                &            +bimX(mtidst(i)+1)*t2z(i) &
                &            +bimX(mtidst(i)+2)*nnz(i)

                exE3x_EM(i) = exE1x_EM(i) + exE2x_EM(i)
                exE3y_EM(i) = exE1y_EM(i) + exE2y_EM(i)
                exE3z_EM(i) = exE1z_EM(i) + exE2z_EM(i)

                exE2xdnn_EM(i) = bimX(mtidst(i)+3)*t1x(i) &
                &               +bimX(mtidst(i)+4)*t2x(i) &
                &               +bimX(mtidst(i)+5)*nnx(i)
                exE2ydnn_EM(i) = bimX(mtidst(i)+3)*t1y(i) &
                &               +bimX(mtidst(i)+4)*t2y(i) &
                &               +bimX(mtidst(i)+5)*nny(i)
                exE2zdnn_EM(i) = bimX(mtidst(i)+3)*t1z(i) &
                &               +bimX(mtidst(i)+4)*t2z(i) &
                &               +bimX(mtidst(i)+5)*nnz(i)

                exE3xdnn_EM(i) = exE1xdnn_EM(i) + exE2xdnn_EM(i)
                exE3ydnn_EM(i) = exE1ydnn_EM(i) + exE2ydnn_EM(i)
                exE3zdnn_EM(i) = exE1zdnn_EM(i) + exE2zdnn_EM(i)

                ztp1 = (exE3x_EM(i)       -inE1x_EM(i))*t1x(i)&
                &     +(exE3y_EM(i)       -inE1y_EM(i))*t1y(i)&
                &     +(exE3z_EM(i)       -inE1z_EM(i))*t1z(i)
                ztp2 = (exE3x_EM(i)       -inE1x_EM(i))*t2x(i)&
                &     +(exE3y_EM(i)       -inE1y_EM(i))*t2y(i)&
                &     +(exE3z_EM(i)       -inE1z_EM(i))*t2z(i)
                ztp3 = (exE3x_EM(i)*mdm_oi-inE1x_EM(i))*nnx(i)&
                &     +(exE3y_EM(i)*mdm_oi-inE1y_EM(i))*nny(i)&
                &     +(exE3z_EM(i)*mdm_oi-inE1z_EM(i))*nnz(i)

                inE2x_EM(i) = ztp1*t1x(i) &
                &            +ztp2*t2x(i) &
                &            +ztp3*nnx(i)
                inE2y_EM(i) = ztp1*t1y(i) &
                &            +ztp2*t2y(i) &
                &            +ztp3*nny(i)
                inE2z_EM(i) = ztp1*t1z(i) &
                &            +ztp2*t2z(i) &
                &            +ztp3*nnz(i)

                inE3x_EM(i) = inE1x_EM(i) + inE2x_EM(i)
                inE3y_EM(i) = inE1y_EM(i) + inE2y_EM(i)
                inE3z_EM(i) = inE1z_EM(i) + inE2z_EM(i)

                ztp1 = ztpzero
                ztp2 = ztpzero
                ztp3 = ztpzero
                DO kk = 1, mxnmbrndlnknd2ndslf
                    IF (ndlnknd2ndslf(i,kk) /= 0) THEN

                        NdA = ndlnknd2ndslf(i,kk)
                        ii = 6*(kk-1)

                        DO id_tp = 1,  6

                            ztp1 = ztp1 + bimX(mtidst(NdA)+id_tp-1) * TanAAt1(i,ii+id_tp)
                            ztp2 = ztp2 + bimX(mtidst(NdA)+id_tp-1) * TanAAt2(i,ii+id_tp)
                            ztp3 = ztp3 + bimX(mtidst(NdA)+id_tp-1) * DivFrAA(i,ii+id_tp)

                        END DO

                    END IF
                END DO

                ztp1 = ztp1 + TanBBt1(i)
                ztp2 = ztp2 + TanBBt2(i)
                ztp3 = ztp3 + DivFrBB(i)

                inE2xdnn_EM(i) = ztp1*t1x(i) &
                &               +ztp2*t2x(i) &
                &               +ztp3*nnx(i)
                inE2ydnn_EM(i) = ztp1*t1y(i) &
                &               +ztp2*t2y(i) &
                &               +ztp3*nny(i)
                inE2zdnn_EM(i) = ztp1*t1z(i) &
                &               +ztp2*t2z(i) &
                &               +ztp3*nnz(i)

                inE3xdnn_EM(i) = inE1xdnn_EM(i) + inE2xdnn_EM(i)
                inE3ydnn_EM(i) = inE1ydnn_EM(i) + inE2ydnn_EM(i)
                inE3zdnn_EM(i) = inE1zdnn_EM(i) + inE2zdnn_EM(i)

            END IF

        END DO

!$OMP END DO
!$OMP END PARALLEL

        DEALLOCATE (bimX,TanAAt1,TanAAt2,TanBBt1,TanBBt2,DivFrAA,DivFrBB,mtidst)
        PRINT *, 'E OK by Esc'

        IF (PostProcdn_on == 1) THEN

            DO icmpt = 1, 3

                Ttlbim = ttlnmbrnd
                ALLOCATE (mtidst(ttlnmbrnd))

                MatSize = 0
                DO ithprtl = 1, nmbrprtl
                    IF (BCType_EM(ithprtl) == 'PEC') THEN
                        id_tp = 1
                        DO i = ndstaID(ithprtl), ndendID(ithprtl)
                            mtidst(i) = MatSize + id_tp
                            id_tp = id_tp + 1
                        END DO
                        MatSize = MatSize + 1*nmbrnd(ithprtl)
                    END IF
                    IF (BCType_EM(ithprtl) == '2SD') THEN
                        id_tp = 1
                        DO i = ndstaID(ithprtl), ndendID(ithprtl)
                            mtidst(i) = MatSize + id_tp
                            id_tp = id_tp + 2
                        END DO
                        MatSize = MatSize + 2*nmbrnd(ithprtl)
                    END IF
                END DO

                ALLOCATE(bimAA(MatSize,MatSize))
                ALLOCATE(bimBB(MatSize))
                ALLOCATE(bimX(MatSize))

!$OMP PARALLEL PRIVATE (i,j)
!$OMP DO
                DO i = 1, MatSize
                    DO j = 1, MatSize
                        bimAA(i,j) = ztpzero
                    END DO
                    bimBB(i) = ztpzero
                    bimX(i) = ztpzero
                END DO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE (i,j,ii,jj,kk,NdA) &
!$OMP & PRIVATE (ithprtl,jthprtl,id_tp) &
!$OMP & PRIVATE (tp,tp1,tp2,tp3,ztp,ztp1,ztp2,ztp3,mdm_oi,mdm_io)
!$OMP DO
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

!***
! Exterior-field equations.
                    IF (     BCType_EM(ithprtl) == 'PEC' &
                    &   .OR. BCType_EM(ithprtl) == '2SD') THEN

!---
! Helmholtz equations for E, rows 1--3 of the assembled system.
!Part 1: surfaces share the same 'corelnkshell'
!!!NOTE: 'E' is 'exE2'

                        DO jthprtl = 1, nmbrprtl

                            id_tp = corelnkshell(jthprtl)

                            IF (corelnkshell(jthprtl) == corelnkshell(ithprtl)) THEN

                                DO j = ndstaID(jthprtl), ndendID(jthprtl)

                                    IF (     BCType_EM(jthprtl) == 'PEC' &
                                    &   .OR. BCType_EM(jthprtl) == '2SD') THEN

                                        !∂(exH2x)/∂n,∂(exH2y)/∂n,∂(exH2z)/∂n
                                        bimAA(mtidst(i)  ,mtidst(j)  ) &
                                    &=  bimAA(mtidst(i)  ,mtidst(j)  ) + bimGGEX(i,j)

                                        IF (icmpt == 1) &
                                        &   bimBB(mtidst(i)  ) = bimBB(mtidst(i)  ) + bimHHEX(i,j)*exE2x_EM(j)
                                        IF (icmpt == 2) &
                                        &   bimBB(mtidst(i)  ) = bimBB(mtidst(i)  ) + bimHHEX(i,j)*exE2y_EM(j)
                                        IF (icmpt == 3) &
                                        &   bimBB(mtidst(i)  ) = bimBB(mtidst(i)  ) + bimHHEX(i,j)*exE2z_EM(j)

                                    END IF



                                END DO

                            END IF

                        END DO

!---

!---
! Continue the E rows 1--3 for bounded exterior-domain surfaces.
!!!NOTE: on this bounded surface, 'E' is 'inE2'
!!!NOTE: ROW 4 to 6 with GGEX and HHEX in the big Matrix (Table 1) in Die_MS II
!!!NOTE: Eq IDs are still 1 to 3, such as mtidst(i), mtidst(i)+1, mtidst(i)+2

                        IF (corelnkshell(ithprtl) > 0) THEN

                            jthprtl = corelnkshell(ithprtl)

                            id_tp = corelnkshell(jthprtl)
                            IF (id_tp == 0) THEN
                                mdm_oi = exeps_EM/ineps_EM(jthprtl)
                            ELSE
                                mdm_oi = ineps_EM(id_tp)/ineps_EM(jthprtl)
                            END IF

                            DO j = ndstaID(jthprtl), ndendID(jthprtl)

                                IF (BCType_EM(jthprtl) == '2SD') THEN

                                    !∂(exH2x)/∂n,∂(exH2y)/∂n,∂(exH2z)/∂n
                                    bimAA(mtidst(i)  ,mtidst(j)  ) &
                                &=  bimAA(mtidst(i)  ,mtidst(j)  ) + bimGGEX(i,j)

                                    IF (icmpt == 1) &
                                    &   bimBB(mtidst(i)  ) = bimBB(mtidst(i)  ) + bimHHEX(i,j)*inE2x_EM(j)
                                    IF (icmpt == 2) &
                                    &   bimBB(mtidst(i)  ) = bimBB(mtidst(i)  ) + bimHHEX(i,j)*inE2y_EM(j)
                                    IF (icmpt == 3) &
                                    &   bimBB(mtidst(i)  ) = bimBB(mtidst(i)  ) + bimHHEX(i,j)*inE2z_EM(j)

                                END IF

                            END DO

                        END IF

!---

                    END IF

!End External field
!***

!***
! Interior-field equations.
                    IF (BCType_EM(ithprtl) == '2SD') THEN

!---
! Helmholtz equations for E, rows 4--6 of the assembled system.

                        DO jthprtl = 1, nmbrprtl

                            IF (jthprtl == ithprtl) THEN

                                id_tp = corelnkshell(jthprtl)
                                IF (id_tp == 0) THEN
                                    mdm_oi = exeps_EM/ineps_EM(jthprtl)
                                ELSE
                                    mdm_oi = ineps_EM(id_tp)/ineps_EM(jthprtl)
                                END IF

                                DO j = ndstaID(jthprtl), ndendID(jthprtl)

                                    !∂(inH2x)/∂n,∂(inH2y)/∂n,∂(inH2z)/∂n
                                    bimAA(mtidst(i)+1,mtidst(j)+1) &
                                &=  bimAA(mtidst(i)+1,mtidst(j)+1) + bimGGIN(i,j)

                                    IF (icmpt == 1) &
                                    &   bimBB(mtidst(i)+1) = bimBB(mtidst(i)+1) + bimHHIN(i,j)*inE2x_EM(j)
                                    IF (icmpt == 2) &
                                    &   bimBB(mtidst(i)+1) = bimBB(mtidst(i)+1) + bimHHIN(i,j)*inE2y_EM(j)
                                    IF (icmpt == 3) &
                                    &   bimBB(mtidst(i)+1) = bimBB(mtidst(i)+1) + bimHHIN(i,j)*inE2z_EM(j)

                                END DO

                            END IF

!---

!---
!Part2: the surfaces that enclosed by the surface where x0 locates
!(just one level down needed)
!!NOTE: 'E' is 'exE2'
!!NOTE: ROW 1 to 3 with GGIN and HHIN in the big Matrix (Table 1) in Die_MS II
!!NOTE: Eq IDs are still 4 to 6, such as mtidst(i)+3, mtidst(i)+4, mtidst(i)+5

                            IF (corelnkshell(jthprtl) == ithprtl) THEN

                                DO j = ndstaID(jthprtl), ndendID(jthprtl)

                                    IF (     BCType_EM(jthprtl) == 'PEC' &
                                    &   .OR. BCType_EM(jthprtl) == '2SD') THEN

                                        !∂(inH2x)/∂n,∂(inH2y)/∂n,∂(inH2z)/∂n
                                        bimAA(mtidst(i)+1,mtidst(j)+1) &
                                    &=  bimAA(mtidst(i)+1,mtidst(j)+1) + bimGGIN(i,j)

                                        IF (icmpt == 1) &
                                        &   bimBB(mtidst(i)+1) = bimBB(mtidst(i)+1) + bimHHIN(i,j)*exE2x_EM(j)
                                        IF (icmpt == 2) &
                                        &   bimBB(mtidst(i)+1) = bimBB(mtidst(i)+1) + bimHHIN(i,j)*exE2y_EM(j)
                                        IF (icmpt == 3) &
                                        &   bimBB(mtidst(i)+1) = bimBB(mtidst(i)+1) + bimHHIN(i,j)*exE2z_EM(j)

                                    END IF

                                END DO

                            END IF

                        END DO

!---

                    END IF

!End Internal field
!***

                END DO
!$OMP END DO
!$OMP END PARALLEL

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!second argument = 0 -> LU decomposition
                CALL Slv_LinearMatrixCOMP_LU(MatSize, 0)

!$OMP PARALLEL PRIVATE (i,j,ii,jj,kk,NdA) &
!$OMP & PRIVATE (ithprtl,jthprtl,id_tp) &
!$OMP & PRIVATE (tp,tp1,tp2,tp3,ztp,ztp1,ztp2,ztp3,mdm_oi,mdm_io)
!$OMP DO

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

                    IF (     BCType_EM(ithprtl) == 'PEC' &
                    &   .OR. BCType_EM(ithprtl) == '2SD') THEN

                        IF (icmpt == 1) THEN
                            exE2xdnn_EM(i) = bimX(mtidst(i)  )
                            exE3xdnn_EM(i) = exE1xdnn_EM(i) + exE2xdnn_EM(i)
                        END IF
                        IF (icmpt == 2) THEN
                            exE2ydnn_EM(i) = bimX(mtidst(i)  )
                            exE3ydnn_EM(i) = exE1ydnn_EM(i) + exE2ydnn_EM(i)
                        END IF
                        IF (icmpt == 3) THEN
                            exE2zdnn_EM(i) = bimX(mtidst(i)  )
                            exE3zdnn_EM(i) = exE1zdnn_EM(i) + exE2zdnn_EM(i)
                        END IF

                    END IF

                    IF (BCType_EM(ithprtl) == '2SD') THEN

                        IF (icmpt == 1) THEN
                            inE2xdnn_EM(i) = bimX(mtidst(i)+1)
                            inE3xdnn_EM(i) = inE1xdnn_EM(i) + inE2xdnn_EM(i)
                        END IF
                        IF (icmpt == 2) THEN
                            inE2ydnn_EM(i) = bimX(mtidst(i)+1)
                            inE3ydnn_EM(i) = inE1ydnn_EM(i) + inE2ydnn_EM(i)
                        END IF
                        IF (icmpt == 3) THEN
                            inE2zdnn_EM(i) = bimX(mtidst(i)+1)
                            inE3zdnn_EM(i) = inE1zdnn_EM(i) + inE2zdnn_EM(i)
                        END IF

                    END IF

                END DO

!$OMP END DO
!$OMP END PARALLEL

                DEALLOCATE (bimX,mtidst)
                PRINT *, 'Edn OK by post-processing'

            END DO

        END IF

        DEALLOCATE (bimGGEX, bimHHEX, bimGGIN, bimHHIN)

    END SUBROUTINE



    ! Assemble and solve the complementary magnetic-field boundary system.
    ! Results are written to the exterior/interior H arrays and their
    ! normal-derivative arrays.

    SUBROUTINE SlvPrblm_EM_Hsc

        INTEGER :: ithprtl, jthprtl, kthprtl, i, j, k, ii, jj, kk, &
        &          icmpt, id_tp, slf_i, slf_j, slf_k, GLQi, icnt

        INTEGER :: NdA, NdB, NdC, NdD, NdE, NdF, ndst, nded, elst, eled

        INTEGER :: Ttlbim, MatSize
        INTEGER, ALLOCATABLE, DIMENSION (:) :: mtidst

        COMPLEX(KIND=KIND(1.0D0)) :: mdm_oi, mdm_io

        DOUBLE PRECISION :: tp, tp1, tp2, tp3, tp4, tp5, tp6, tp7, tp8, tp9
        COMPLEX(KIND=KIND(1.0D0)) :: ztp,ztp1,ztp2,ztp3,ztp4,ztp5,ztp6,ztp7,ztp8,ztp9

        ALLOCATE (DivFrAA(ttlnmbrnd,6*ttlddtnd))
        ALLOCATE (DivFrBB(ttlnmbrnd))

!$OMP PARALLEL PRIVATE (i,j)
!$OMP DO
        DO i = 1, ttlnmbrnd
            DO j = 1, 6*ttlddtnd
                DivFrAA(i,j) = ztpzero
            END DO
            DivFrBB(i) = ztpzero
        END DO
!$OMP END DO
!$OMP END PARALLEL

        ALLOCATE (TanAAt1(ttlnmbrnd,6*ttlddtnd))
        ALLOCATE (TanAAt2(ttlnmbrnd,6*ttlddtnd))
        ALLOCATE (TanBBt1(ttlnmbrnd))
        ALLOCATE (TanBBt2(ttlnmbrnd))

!$OMP PARALLEL PRIVATE (i,j)
!$OMP DO
        DO i = 1, ttlnmbrnd
            DO j = 1, 6*ttlddtnd
                TanAAt1(i,j) = ztpzero
                TanAAt2(i,j) = ztpzero
            END DO
            TanBBt1(i) = ztpzero
            TanBBt2(i) = ztpzero
        END DO
!$OMP END DO
!$OMP END PARALLEL

        DO ithprtl = 1, nmbrprtl
            IF (BCType_EM(ithprtl) == 'PEC') THEN
                CALL Get_PEC_DivFr(ithprtl, 'H')
                CALL Get_PEC_TanCn(ithprtl, 'H')
            END IF
            IF (BCType_EM(ithprtl) == '2SD') THEN
                CALL Get_DEL_DivFrsc(ithprtl, 'H')
                CALL Get_DEL_TanCnsc(ithprtl, 'H')
            END IF
        END DO


        Ttlbim = ttlnmbrnd
        ALLOCATE (mtidst(ttlnmbrnd))

        MatSize = 0
        DO ithprtl = 1, nmbrprtl
            IF (BCType_EM(ithprtl) == 'PEC') THEN
                id_tp = 1
                DO i = ndstaID(ithprtl), ndendID(ithprtl)
                    mtidst(i) = MatSize + id_tp
                    id_tp = id_tp + 3
                END DO
                MatSize = MatSize + 3*nmbrnd(ithprtl)
            END IF
            IF (BCType_EM(ithprtl) == '2SD') THEN
                id_tp = 1
                DO i = ndstaID(ithprtl), ndendID(ithprtl)
                    mtidst(i) = MatSize + id_tp
                    id_tp = id_tp + 6
                END DO
                MatSize = MatSize + 6*nmbrnd(ithprtl)
            END IF
        END DO

        ALLOCATE(bimAA(MatSize,MatSize))
        ALLOCATE(bimBB(MatSize))
        ALLOCATE(bimX(MatSize))

!$OMP PARALLEL PRIVATE (i,j)
!$OMP DO
        DO i = 1, MatSize
            DO j = 1, MatSize
                bimAA(i,j) = ztpzero
            END DO
            bimBB(i) = ztpzero
            bimX(i) = ztpzero
        END DO
!$OMP END DO
!$OMP END PARALLEL

        CALL Get_BRIEFGGHHEX_EM
        CALL Get_BRIEFGGHHIN_EM

!$OMP PARALLEL PRIVATE (i,j,ii,jj,kk,NdA) &
!$OMP & PRIVATE (ithprtl,jthprtl,id_tp) &
!$OMP & PRIVATE (tp,tp1,tp2,tp3,ztp,ztp1,ztp2,ztp3,mdm_oi,mdm_io)
!$OMP DO
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

!***
! Exterior-field equations.
            IF (     BCType_EM(ithprtl) == 'PEC' &
            &   .OR. BCType_EM(ithprtl) == '2SD') THEN

!---
! Helmholtz equations for H, rows 1--3 of the assembled system.
!Part 1: surfaces share the same 'corelnkshell'
!!!NOTE: 'H' is 'exH2'

                DO jthprtl = 1, nmbrprtl

                    id_tp = corelnkshell(jthprtl)

                    IF (corelnkshell(jthprtl) == corelnkshell(ithprtl)) THEN

                        DO j = ndstaID(jthprtl), ndendID(jthprtl)

                            IF (BCType_EM(jthprtl) == 'PEC') THEN

                                !exH2t1
                                bimAA(mtidst(i)  ,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)  ,mtidst(j)  ) + bimHHEX(i,j)*t1x(j)
                                bimAA(mtidst(i)+1,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)+1,mtidst(j)  ) + bimHHEX(i,j)*t1y(j)
                                bimAA(mtidst(i)+2,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)+2,mtidst(j)  ) + bimHHEX(i,j)*t1z(j)

                                !exH2t2
                                bimAA(mtidst(i)  ,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)  ,mtidst(j)+1) + bimHHEX(i,j)*t2x(j)
                                bimAA(mtidst(i)+1,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)+1,mtidst(j)+1) + bimHHEX(i,j)*t2y(j)
                                bimAA(mtidst(i)+2,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)+2,mtidst(j)+1) + bimHHEX(i,j)*t2z(j)

                                !n •∂(exH2)/∂n
                                bimAA(mtidst(i)  ,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)  ,mtidst(j)+2) - bimGGEX(i,j)*nnx(j)
                                bimAA(mtidst(i)+1,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)+1,mtidst(j)+2) - bimGGEX(i,j)*nny(j)
                                bimAA(mtidst(i)+2,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)+2,mtidst(j)+2) - bimGGEX(i,j)*nnz(j)

                                ![t1•∂(exH2)/∂n,t2•∂(exH2)/∂n] by BC: [exH2t1,exH2t2]
                                DO kk = 1, mxnmbrndlnknd2ndslf
                                    IF (ndlnknd2ndslf(j,kk) /= 0) THEN

                                        NdA = ndlnknd2ndslf(j,kk)
                                        ii = 6*(kk-1)

                                        bimAA(mtidst(i)  ,mtidst(NdA)  ) &
                                    &=  bimAA(mtidst(i)  ,mtidst(NdA)  ) &
                                    &       -bimGGEX(i,j)*t1x(j)*TanAAt1(j,ii+1)
                                        bimAA(mtidst(i)+1,mtidst(NdA)  ) &
                                    &=  bimAA(mtidst(i)+1,mtidst(NdA)  ) &
                                    &       -bimGGEX(i,j)*t1y(j)*TanAAt1(j,ii+1)
                                        bimAA(mtidst(i)+2,mtidst(NdA)  ) &
                                    &=  bimAA(mtidst(i)+2,mtidst(NdA)  ) &
                                    &       -bimGGEX(i,j)*t1z(j)*TanAAt1(j,ii+1)
                                        bimAA(mtidst(i)  ,mtidst(NdA)+1) &
                                    &=  bimAA(mtidst(i)  ,mtidst(NdA)+1) &
                                    &       -bimGGEX(i,j)*t1x(j)*TanAAt1(j,ii+2)
                                        bimAA(mtidst(i)+1,mtidst(NdA)+1) &
                                    &=  bimAA(mtidst(i)+1,mtidst(NdA)+1) &
                                    &       -bimGGEX(i,j)*t1y(j)*TanAAt1(j,ii+2)
                                        bimAA(mtidst(i)+2,mtidst(NdA)+1) &
                                    &=  bimAA(mtidst(i)+2,mtidst(NdA)+1) &
                                    &       -bimGGEX(i,j)*t1z(j)*TanAAt1(j,ii+2)

                                        bimAA(mtidst(i)  ,mtidst(NdA)  ) &
                                    &=  bimAA(mtidst(i)  ,mtidst(NdA)  ) &
                                    &       -bimGGEX(i,j)*t2x(j)*TanAAt2(j,ii+1)
                                        bimAA(mtidst(i)+1,mtidst(NdA)  ) &
                                    &=  bimAA(mtidst(i)+1,mtidst(NdA)  ) &
                                    &       -bimGGEX(i,j)*t2y(j)*TanAAt2(j,ii+1)
                                        bimAA(mtidst(i)+2,mtidst(NdA)  ) &
                                    &=  bimAA(mtidst(i)+2,mtidst(NdA)  ) &
                                    &       -bimGGEX(i,j)*t2z(j)*TanAAt2(j,ii+1)
                                        bimAA(mtidst(i)  ,mtidst(NdA)+1) &
                                    &=  bimAA(mtidst(i)  ,mtidst(NdA)+1) &
                                    &       -bimGGEX(i,j)*t2x(j)*TanAAt2(j,ii+2)
                                        bimAA(mtidst(i)+1,mtidst(NdA)+1) &
                                    &=  bimAA(mtidst(i)+1,mtidst(NdA)+1) &
                                    &       -bimGGEX(i,j)*t2y(j)*TanAAt2(j,ii+2)
                                        bimAA(mtidst(i)+2,mtidst(NdA)+1) &
                                    &=  bimAA(mtidst(i)+2,mtidst(NdA)+1) &
                                    &       -bimGGEX(i,j)*t2z(j)*TanAAt2(j,ii+2)

                                    END IF
                                END DO

                                ![t1•∂(exH2)/∂n,t2•∂(exH2)/∂n] by BC: [exH1t1,exH1t2,exH1nn]
                                bimBB(mtidst(i)  ) = bimBB(mtidst(i)  ) &
                                &                   +bimGGEX(i,j)*( t1x(j)*TanBBt1(j) &
                                &                                  +t2x(j)*TanBBt2(j) )
                                bimBB(mtidst(i)+1) = bimBB(mtidst(i)+1) &
                                &                   +bimGGEX(i,j)*( t1y(j)*TanBBt1(j) &
                                &                                  +t2y(j)*TanBBt2(j) )
                                bimBB(mtidst(i)+2) = bimBB(mtidst(i)+2) &
                                &                   +bimGGEX(i,j)*( t1z(j)*TanBBt1(j) &
                                &                                  +t2z(j)*TanBBt2(j) )

                                ztp1 = exH1x_EM(j)*nnx(j)&
                                &     +exH1y_EM(j)*nny(j)&
                                &     +exH1z_EM(j)*nnz(j)

                                bimBB(mtidst(i)  ) = bimBB(mtidst(i)  ) &
                                &       +bimHHEX(i,j) * ztp1*nnx(j)

                                bimBB(mtidst(i)+1) = bimBB(mtidst(i)+1) &
                                &       +bimHHEX(i,j) * ztp1*nny(j)

                                bimBB(mtidst(i)+2) = bimBB(mtidst(i)+2) &
                                &       +bimHHEX(i,j) * ztp1*nnz(j)

                            END IF

                            IF (BCType_EM(jthprtl) == '2SD') THEN

                                !exH2t1
                                bimAA(mtidst(i)  ,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)  ,mtidst(j)  ) + bimHHEX(i,j)*t1x(j)
                                bimAA(mtidst(i)+1,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)+1,mtidst(j)  ) + bimHHEX(i,j)*t1y(j)
                                bimAA(mtidst(i)+2,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)+2,mtidst(j)  ) + bimHHEX(i,j)*t1z(j)

                                !exH2t2
                                bimAA(mtidst(i)  ,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)  ,mtidst(j)+1) + bimHHEX(i,j)*t2x(j)
                                bimAA(mtidst(i)+1,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)+1,mtidst(j)+1) + bimHHEX(i,j)*t2y(j)
                                bimAA(mtidst(i)+2,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)+2,mtidst(j)+1) + bimHHEX(i,j)*t2z(j)

                                !exH2nn
                                bimAA(mtidst(i)  ,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)  ,mtidst(j)+2) + bimHHEX(i,j)*nnx(j)
                                bimAA(mtidst(i)+1,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)+1,mtidst(j)+2) + bimHHEX(i,j)*nny(j)
                                bimAA(mtidst(i)+2,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)+2,mtidst(j)+2) + bimHHEX(i,j)*nnz(j)

                                !t1•∂(exH2)/∂n
                                bimAA(mtidst(i)  ,mtidst(j)+3) &
                            &=  bimAA(mtidst(i)  ,mtidst(j)+3) - bimGGEX(i,j)*t1x(j)
                                bimAA(mtidst(i)+1,mtidst(j)+3) &
                            &=  bimAA(mtidst(i)+1,mtidst(j)+3) - bimGGEX(i,j)*t1y(j)
                                bimAA(mtidst(i)+2,mtidst(j)+3) &
                            &=  bimAA(mtidst(i)+2,mtidst(j)+3) - bimGGEX(i,j)*t1z(j)

                                !t2•∂(exH2)/∂n
                                bimAA(mtidst(i)  ,mtidst(j)+4) &
                            &=  bimAA(mtidst(i)  ,mtidst(j)+4) - bimGGEX(i,j)*t2x(j)
                                bimAA(mtidst(i)+1,mtidst(j)+4) &
                            &=  bimAA(mtidst(i)+1,mtidst(j)+4) - bimGGEX(i,j)*t2y(j)
                                bimAA(mtidst(i)+2,mtidst(j)+4) &
                            &=  bimAA(mtidst(i)+2,mtidst(j)+4) - bimGGEX(i,j)*t2z(j)

                                !n•∂(exH2)/∂n
                                bimAA(mtidst(i)  ,mtidst(j)+5) &
                            &=  bimAA(mtidst(i)  ,mtidst(j)+5) - bimGGEX(i,j)*nnx(j)
                                bimAA(mtidst(i)+1,mtidst(j)+5) &
                            &=  bimAA(mtidst(i)+1,mtidst(j)+5) - bimGGEX(i,j)*nny(j)
                                bimAA(mtidst(i)+2,mtidst(j)+5) &
                            &=  bimAA(mtidst(i)+2,mtidst(j)+5) - bimGGEX(i,j)*nnz(j)

                            END IF

                        END DO

                    END IF

                END DO

!---

!---
! Continue the H rows 1--3 for bounded exterior-domain surfaces.
!!!NOTE: on this bounded surface, 'H' is 'inH2'
!!!NOTE: ROW 4 to 6 with GGEX and HHEX in the big Matrix (Table 1) in Die_MS II
!!!NOTE: Eq IDs are still 1 to 3, such as mtidst(i), mtidst(i)+1, mtidst(i)+2

                IF (corelnkshell(ithprtl) > 0) THEN

                    jthprtl = corelnkshell(ithprtl)

                    id_tp = corelnkshell(jthprtl)
                    IF (id_tp == 0) THEN
                        mdm_oi = exmiu_EM/inmiu_EM(jthprtl)
                    ELSE
                        mdm_oi = inmiu_EM(id_tp)/inmiu_EM(jthprtl)
                    END IF

                    DO j = ndstaID(jthprtl), ndendID(jthprtl)

                        IF (BCType_EM(jthprtl) == '2SD') THEN

                            bimAA(mtidst(i)  ,mtidst(j)  ) &
                        &=  bimAA(mtidst(i)  ,mtidst(j)  ) + bimHHEX(i,j)*t1x(j)
                            bimAA(mtidst(i)+1,mtidst(j)  ) &
                        &=  bimAA(mtidst(i)+1,mtidst(j)  ) + bimHHEX(i,j)*t1y(j)
                            bimAA(mtidst(i)+2,mtidst(j)  ) &
                        &=  bimAA(mtidst(i)+2,mtidst(j)  ) + bimHHEX(i,j)*t1z(j)

                            bimAA(mtidst(i)  ,mtidst(j)+1) &
                        &=  bimAA(mtidst(i)  ,mtidst(j)+1) + bimHHEX(i,j)*t2x(j)
                            bimAA(mtidst(i)+1,mtidst(j)+1) &
                        &=  bimAA(mtidst(i)+1,mtidst(j)+1) + bimHHEX(i,j)*t2y(j)
                            bimAA(mtidst(i)+2,mtidst(j)+1) &
                        &=  bimAA(mtidst(i)+2,mtidst(j)+1) + bimHHEX(i,j)*t2z(j)

                            bimAA(mtidst(i)  ,mtidst(j)+2) &
                        &=  bimAA(mtidst(i)  ,mtidst(j)+2) + bimHHEX(i,j)*nnx(j)*mdm_oi
                            bimAA(mtidst(i)+1,mtidst(j)+2) &
                        &=  bimAA(mtidst(i)+1,mtidst(j)+2) + bimHHEX(i,j)*nny(j)*mdm_oi
                            bimAA(mtidst(i)+2,mtidst(j)+2) &
                        &=  bimAA(mtidst(i)+2,mtidst(j)+2) + bimHHEX(i,j)*nnz(j)*mdm_oi

![t1•∂(inH2)/∂n,t2•∂(inH2)/∂n,n•∂(inH2)/∂n] by BC:
![exH2t1,exH1t2,exH2nn,t1•∂(exH2)/∂n,t2•∂(exH2)/∂n]
                            DO kk = 1, mxnmbrndlnknd2ndslf
                                IF (ndlnknd2ndslf(j,kk) /= 0) THEN

                                    NdA = ndlnknd2ndslf(j,kk)
                                    ii = 6*(kk-1)

                                    DO id_tp = 1, 6

                                        bimAA(mtidst(i)  ,mtidst(NdA)+id_tp-1) &
                                    &=  bimAA(mtidst(i)  ,mtidst(NdA)+id_tp-1) &
                                    &       -bimGGEX(i,j)*( t1x(j)*TanAAt1(j,ii+id_tp) &
                                    &                      +t2x(j)*TanAAt2(j,ii+id_tp) &
                                    &                      +nnx(j)*DivFrAA(j,ii+id_tp) )

                                        bimAA(mtidst(i)+1,mtidst(NdA)+id_tp-1) &
                                    &=  bimAA(mtidst(i)+1,mtidst(NdA)+id_tp-1) &
                                    &       -bimGGEX(i,j)*( t1y(j)*TanAAt1(j,ii+id_tp) &
                                    &                      +t2y(j)*TanAAt2(j,ii+id_tp) &
                                    &                      +nny(j)*DivFrAA(j,ii+id_tp) )

                                        bimAA(mtidst(i)+2,mtidst(NdA)+id_tp-1) &
                                    &=  bimAA(mtidst(i)+2,mtidst(NdA)+id_tp-1) &
                                    &       -bimGGEX(i,j)*( t1z(j)*TanAAt1(j,ii+id_tp) &
                                    &                      +t2z(j)*TanAAt2(j,ii+id_tp) &
                                    &                      +nnz(j)*DivFrAA(j,ii+id_tp) )

                                    END DO

                                END IF
                            END DO

                            ztp1 = t1x(j)*(exH1x_EM(j)-inH1x_EM(j)) &
                            &     +t1y(j)*(exH1y_EM(j)-inH1y_EM(j)) &
                            &     +t1z(j)*(exH1z_EM(j)-inH1z_EM(j))
                            ztp2 = t2x(j)*(exH1x_EM(j)-inH1x_EM(j)) &
                            &     +t2y(j)*(exH1y_EM(j)-inH1y_EM(j)) &
                            &     +t2z(j)*(exH1z_EM(j)-inH1z_EM(j))
                            ztp3 = nnx(j)* exH1x_EM(j) &
                            &     +nny(j)* exH1y_EM(j) &
                            &     +nnz(j)* exH1z_EM(j)
                            ztp3 = ztp3 * mdm_oi
                            ztp  = nnx(j)* inH1x_EM(j) &
                            &     +nny(j)* inH1y_EM(j) &
                            &     +nnz(j)* inH1z_EM(j)

                            bimBB(mtidst(i)  ) = bimBB(mtidst(i)  ) &
                            &       -bimHHEX(i,j) * ( t1x(j)*ztp1+t2x(j)*ztp2 &
                            &                        +nnx(j)*ztp3-nnx(j)*ztp  )
                            bimBB(mtidst(i)+1) = bimBB(mtidst(i)+1) &
                            &       -bimHHEX(i,j) * ( t1y(j)*ztp1+t2y(j)*ztp2 &
                            &                        +nny(j)*ztp3-nny(j)*ztp  )
                            bimBB(mtidst(i)+2) = bimBB(mtidst(i)+2) &
                            &       -bimHHEX(i,j) * ( t1z(j)*ztp1+t2z(j)*ztp2 &
                            &                        +nnz(j)*ztp3-nnz(j)*ztp  )

![t1•∂(inH2)/∂n,t2•∂(inH2)/∂n,n•∂(inH2)/∂n] by BC:
![exH2t1,exH1t2,exH2nn,t1•∂(exH2)/∂n,t2•∂(exH2)/∂n]
                            ztp1 = t1x(j)*TanBBt1(j) &
                            &     +t2x(j)*TanBBt2(j) &
                            &     +nnx(j)*DivFrBB(j)
                            ztp2 = t1y(j)*TanBBt1(j) &
                            &     +t2y(j)*TanBBt2(j) &
                            &     +nny(j)*DivFrBB(j)
                            ztp3 = t1z(j)*TanBBt1(j) &
                            &     +t2z(j)*TanBBt2(j) &
                            &     +nnz(j)*DivFrBB(j)
                            bimBB(mtidst(i)  ) = bimBB(mtidst(i)  ) + bimGGEX(i,j) * ztp1
                            bimBB(mtidst(i)+1) = bimBB(mtidst(i)+1) + bimGGEX(i,j) * ztp2
                            bimBB(mtidst(i)+2) = bimBB(mtidst(i)+2) + bimGGEX(i,j) * ztp3

                        END IF

                    END DO

                END IF

!---

            END IF

!End External field
!***

!***
! Interior-field equations.
            IF (BCType_EM(ithprtl) == '2SD') THEN

!---
! Helmholtz equations for H, rows 4--6 of the assembled system.

                DO jthprtl = 1, nmbrprtl

                    IF (jthprtl == ithprtl) THEN

                        id_tp = corelnkshell(jthprtl)
                        IF (id_tp == 0) THEN
                            mdm_oi = exmiu_EM/inmiu_EM(jthprtl)
                        ELSE
                            mdm_oi = inmiu_EM(id_tp)/inmiu_EM(jthprtl)
                        END IF

                        DO j = ndstaID(jthprtl), ndendID(jthprtl)

                            bimAA(mtidst(i)+3,mtidst(j)  ) &
                        &=  bimAA(mtidst(i)+3,mtidst(j)  ) + bimHHIN(i,j)*t1x(j)
                            bimAA(mtidst(i)+4,mtidst(j)  ) &
                        &=  bimAA(mtidst(i)+4,mtidst(j)  ) + bimHHIN(i,j)*t1y(j)
                            bimAA(mtidst(i)+5,mtidst(j)  ) &
                        &=  bimAA(mtidst(i)+5,mtidst(j)  ) + bimHHIN(i,j)*t1z(j)

                            bimAA(mtidst(i)+3,mtidst(j)+1) &
                        &=  bimAA(mtidst(i)+3,mtidst(j)+1) + bimHHIN(i,j)*t2x(j)
                            bimAA(mtidst(i)+4,mtidst(j)+1) &
                        &=  bimAA(mtidst(i)+4,mtidst(j)+1) + bimHHIN(i,j)*t2y(j)
                            bimAA(mtidst(i)+5,mtidst(j)+1) &
                        &=  bimAA(mtidst(i)+5,mtidst(j)+1) + bimHHIN(i,j)*t2z(j)

                            bimAA(mtidst(i)+3,mtidst(j)+2) &
                        &=  bimAA(mtidst(i)+3,mtidst(j)+2) + bimHHIN(i,j)*nnx(j)*mdm_oi
                            bimAA(mtidst(i)+4,mtidst(j)+2) &
                        &=  bimAA(mtidst(i)+4,mtidst(j)+2) + bimHHIN(i,j)*nny(j)*mdm_oi
                            bimAA(mtidst(i)+5,mtidst(j)+2) &
                        &=  bimAA(mtidst(i)+5,mtidst(j)+2) + bimHHIN(i,j)*nnz(j)*mdm_oi

![t1•∂(inH2)/∂n,t2•∂(inH2)/∂n,n•∂(inH2)/∂n] by BC:
![exH2t1,exH1t2,exH2nn,t1•∂(exH2)/∂n,t2•∂(exH2)/∂n]
                            DO kk = 1, mxnmbrndlnknd2ndslf
                                IF (ndlnknd2ndslf(j,kk) /= 0) THEN

                                    NdA = ndlnknd2ndslf(j,kk)
                                    ii = 6*(kk-1)

                                    DO id_tp = 1, 6

                                        bimAA(mtidst(i)+3,mtidst(NdA)+id_tp-1) &
                                    &=  bimAA(mtidst(i)+3,mtidst(NdA)+id_tp-1) &
                                    &       -bimGGIN(i,j)*( t1x(j)*TanAAt1(j,ii+id_tp) &
                                    &                      +t2x(j)*TanAAt2(j,ii+id_tp) &
                                    &                      +nnx(j)*DivFrAA(j,ii+id_tp) )

                                        bimAA(mtidst(i)+4,mtidst(NdA)+id_tp-1) &
                                    &=  bimAA(mtidst(i)+4,mtidst(NdA)+id_tp-1) &
                                    &       -bimGGIN(i,j)*( t1y(j)*TanAAt1(j,ii+id_tp) &
                                    &                      +t2y(j)*TanAAt2(j,ii+id_tp) &
                                    &                      +nny(j)*DivFrAA(j,ii+id_tp) )

                                        bimAA(mtidst(i)+5,mtidst(NdA)+id_tp-1) &
                                    &=  bimAA(mtidst(i)+5,mtidst(NdA)+id_tp-1) &
                                    &       -bimGGIN(i,j)*( t1z(j)*TanAAt1(j,ii+id_tp) &
                                    &                      +t2z(j)*TanAAt2(j,ii+id_tp) &
                                    &                      +nnz(j)*DivFrAA(j,ii+id_tp) )

                                    END DO

                                END IF
                            END DO

                            ztp1 = t1x(j)*(exH1x_EM(j)-inH1x_EM(j)) &
                            &     +t1y(j)*(exH1y_EM(j)-inH1y_EM(j)) &
                            &     +t1z(j)*(exH1z_EM(j)-inH1z_EM(j))
                            ztp2 = t2x(j)*(exH1x_EM(j)-inH1x_EM(j)) &
                            &     +t2y(j)*(exH1y_EM(j)-inH1y_EM(j)) &
                            &     +t2z(j)*(exH1z_EM(j)-inH1z_EM(j))
                            ztp3 = nnx(j)* exH1x_EM(j) &
                            &     +nny(j)* exH1y_EM(j) &
                            &     +nnz(j)* exH1z_EM(j)
                            ztp3 = ztp3 * mdm_oi
                            ztp  = nnx(j)* inH1x_EM(j) &
                            &     +nny(j)* inH1y_EM(j) &
                            &     +nnz(j)* inH1z_EM(j)

                            bimBB(mtidst(i)+3) = bimBB(mtidst(i)+3) &
                            &       -bimHHIN(i,j) * ( t1x(j)*ztp1+t2x(j)*ztp2 &
                            &                        +nnx(j)*ztp3-nnx(j)*ztp  )
                            bimBB(mtidst(i)+4) = bimBB(mtidst(i)+4) &
                            &       -bimHHIN(i,j) * ( t1y(j)*ztp1+t2y(j)*ztp2 &
                            &                        +nny(j)*ztp3-nny(j)*ztp  )
                            bimBB(mtidst(i)+5) = bimBB(mtidst(i)+5) &
                            &       -bimHHIN(i,j) * ( t1z(j)*ztp1+t2z(j)*ztp2 &
                            &                        +nnz(j)*ztp3-nnz(j)*ztp  )

![t1•∂(inH2)/∂n,t2•∂(inH2)/∂n,n•∂(inH2)/∂n] by BC:
![exH2t1,exH1t2,exH2nn,t1•∂(exH2)/∂n,t2•∂(exH2)/∂n]
                            ztp1 = t1x(j)*TanBBt1(j) &
                            &     +t2x(j)*TanBBt2(j) &
                            &     +nnx(j)*DivFrBB(j)
                            ztp2 = t1y(j)*TanBBt1(j) &
                            &     +t2y(j)*TanBBt2(j) &
                            &     +nny(j)*DivFrBB(j)
                            ztp3 = t1z(j)*TanBBt1(j) &
                            &     +t2z(j)*TanBBt2(j) &
                            &     +nnz(j)*DivFrBB(j)
                            bimBB(mtidst(i)+3) = bimBB(mtidst(i)+3) + bimGGIN(i,j) * ztp1
                            bimBB(mtidst(i)+4) = bimBB(mtidst(i)+4) + bimGGIN(i,j) * ztp2
                            bimBB(mtidst(i)+5) = bimBB(mtidst(i)+5) + bimGGIN(i,j) * ztp3
                        END DO

                    END IF

!---

!---
!Part2: the surfaces that enclosed by the surface where x0 locates
!(just one level down needed)
!!NOTE: 'H' is 'exH2'
!!NOTE: ROW 1 to 3 with GGIN and HHIN in the big Matrix (Table 1) in Die_MS II
!!NOTE: Eq IDs are still 4 to 6, such as mtidst(i)+3, mtidst(i)+4, mtidst(i)+5

                    IF (corelnkshell(jthprtl) == ithprtl) THEN

                        DO j = ndstaID(jthprtl), ndendID(jthprtl)

                            IF (BCType_EM(jthprtl) == 'PEC') THEN

                                !exH2t1
                                bimAA(mtidst(i)+3,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)+3,mtidst(j)  ) + bimHHIN(i,j)*t1x(j)
                                bimAA(mtidst(i)+4,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)+4,mtidst(j)  ) + bimHHIN(i,j)*t1y(j)
                                bimAA(mtidst(i)+5,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)+5,mtidst(j)  ) + bimHHIN(i,j)*t1z(j)

                                !exH2t2
                                bimAA(mtidst(i)+3,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)+3,mtidst(j)+1) + bimHHIN(i,j)*t2x(j)
                                bimAA(mtidst(i)+4,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)+4,mtidst(j)+1) + bimHHIN(i,j)*t2y(j)
                                bimAA(mtidst(i)+5,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)+5,mtidst(j)+1) + bimHHIN(i,j)*t2z(j)

                                !n •∂(exH2)/∂n
                                bimAA(mtidst(i)+3,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)+3,mtidst(j)+2) - bimGGIN(i,j)*nnx(j)
                                bimAA(mtidst(i)+4,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)+4,mtidst(j)+2) - bimGGIN(i,j)*nny(j)
                                bimAA(mtidst(i)+5,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)+5,mtidst(j)+2) - bimGGIN(i,j)*nnz(j)

                            ![t1•∂(exH2)/∂n,t2•∂(exH2)/∂n] by BC: [exH2t1,exH2t2]
                                DO kk = 1, mxnmbrndlnknd2ndslf
                                    IF (ndlnknd2ndslf(j,kk) /= 0) THEN

                                        NdA = ndlnknd2ndslf(j,kk)
                                        ii = 6*(kk-1)

                                        bimAA(mtidst(i)+3,mtidst(NdA)  ) &
                                    &=  bimAA(mtidst(i)+3,mtidst(NdA)  ) &
                                    &       -bimGGIN(i,j)*t1x(j)*TanAAt1(j,ii+1)
                                        bimAA(mtidst(i)+4,mtidst(NdA)  ) &
                                    &=  bimAA(mtidst(i)+4,mtidst(NdA)  ) &
                                    &       -bimGGIN(i,j)*t1y(j)*TanAAt1(j,ii+1)
                                        bimAA(mtidst(i)+5,mtidst(NdA)  ) &
                                    &=  bimAA(mtidst(i)+5,mtidst(NdA)  ) &
                                    &       -bimGGIN(i,j)*t1z(j)*TanAAt1(j,ii+1)

                                        bimAA(mtidst(i)+3,mtidst(NdA)+1) &
                                    &=  bimAA(mtidst(i)+3,mtidst(NdA)+1) &
                                    &       -bimGGIN(i,j)*t2x(j)*TanAAt2(j,ii+2)
                                        bimAA(mtidst(i)+4,mtidst(NdA)+1) &
                                    &=  bimAA(mtidst(i)+4,mtidst(NdA)+1) &
                                    &       -bimGGIN(i,j)*t2y(j)*TanAAt2(j,ii+2)
                                        bimAA(mtidst(i)+5,mtidst(NdA)+1) &
                                    &=  bimAA(mtidst(i)+5,mtidst(NdA)+1) &
                                    &       -bimGGIN(i,j)*t2z(j)*TanAAt2(j,ii+2)

                                    END IF
                                END DO

                                ![t1•∂(exH2)/∂n,t2•∂(exH2)/∂n] by BC: [exH1t1,exH1t2,exH1nn]
                                bimBB(mtidst(i)+3) = bimBB(mtidst(i)+3) &
                                &                   +bimGGIN(i,j)*( t1x(j)*TanBBt1(j) &
                                &                                  +t2x(j)*TanBBt2(j) )
                                bimBB(mtidst(i)+4) = bimBB(mtidst(i)+4) &
                                &                   +bimGGIN(i,j)*( t1y(j)*TanBBt1(j) &
                                &                                  +t2y(j)*TanBBt2(j) )
                                bimBB(mtidst(i)+5) = bimBB(mtidst(i)+5) &
                                &                   +bimGGIN(i,j)*( t1z(j)*TanBBt1(j) &
                                &                                  +t2z(j)*TanBBt2(j) )

                                ztp1 = exH1x_EM(j)*nnx(j)&
                                &     +exH1y_EM(j)*nny(j)&
                                &     +exH1z_EM(j)*nnz(j)

                                bimBB(mtidst(i)+3) = bimBB(mtidst(i)+3) &
                                &       +bimHHIN(i,j) * ztp1*nnx(j)

                                bimBB(mtidst(i)+4) = bimBB(mtidst(i)+4) &
                                &       +bimHHIN(i,j) * ztp1*nny(j)

                                bimBB(mtidst(i)+5) = bimBB(mtidst(i)+5) &
                                &       +bimHHIN(i,j) * ztp1*nnz(j)

                            END IF

                            IF (BCType_EM(jthprtl) == '2SD') THEN

                                !exH2t1
                                bimAA(mtidst(i)+3,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)+3,mtidst(j)  ) + bimHHIN(i,j)*t1x(j)
                                bimAA(mtidst(i)+4,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)+4,mtidst(j)  ) + bimHHIN(i,j)*t1y(j)
                                bimAA(mtidst(i)+5,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)+5,mtidst(j)  ) + bimHHIN(i,j)*t1z(j)

                                !exH2t2
                                bimAA(mtidst(i)+3,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)+3,mtidst(j)+1) + bimHHIN(i,j)*t2x(j)
                                bimAA(mtidst(i)+4,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)+4,mtidst(j)+1) + bimHHIN(i,j)*t2y(j)
                                bimAA(mtidst(i)+5,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)+5,mtidst(j)+1) + bimHHIN(i,j)*t2z(j)

                                !exH2nn
                                bimAA(mtidst(i)+3,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)+3,mtidst(j)+2) + bimHHIN(i,j)*nnx(j)
                                bimAA(mtidst(i)+4,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)+4,mtidst(j)+2) + bimHHIN(i,j)*nny(j)
                                bimAA(mtidst(i)+5,mtidst(j)+2) &
                            &=  bimAA(mtidst(i)+5,mtidst(j)+2) + bimHHIN(i,j)*nnz(j)

                                !t1•∂(exH2)/∂n
                                bimAA(mtidst(i)+3,mtidst(j)+3) &
                            &=  bimAA(mtidst(i)+3,mtidst(j)+3) - bimGGIN(i,j)*t1x(j)
                                bimAA(mtidst(i)+4,mtidst(j)+3) &
                            &=  bimAA(mtidst(i)+4,mtidst(j)+3) - bimGGIN(i,j)*t1y(j)
                                bimAA(mtidst(i)+5,mtidst(j)+3) &
                            &=  bimAA(mtidst(i)+5,mtidst(j)+3) - bimGGIN(i,j)*t1z(j)

                                !t2•∂(exH2)/∂n
                                bimAA(mtidst(i)+3,mtidst(j)+4) &
                            &=  bimAA(mtidst(i)+3,mtidst(j)+4) - bimGGIN(i,j)*t2x(j)
                                bimAA(mtidst(i)+4,mtidst(j)+4) &
                            &=  bimAA(mtidst(i)+4,mtidst(j)+4) - bimGGIN(i,j)*t2y(j)
                                bimAA(mtidst(i)+5,mtidst(j)+4) &
                            &=  bimAA(mtidst(i)+5,mtidst(j)+4) - bimGGIN(i,j)*t2z(j)

                                !n•∂(exH2)/∂n
                                bimAA(mtidst(i)+3,mtidst(j)+5) &
                            &=  bimAA(mtidst(i)+3,mtidst(j)+5) - bimGGIN(i,j)*nnx(j)
                                bimAA(mtidst(i)+4,mtidst(j)+5) &
                            &=  bimAA(mtidst(i)+4,mtidst(j)+5) - bimGGIN(i,j)*nny(j)
                                bimAA(mtidst(i)+5,mtidst(j)+5) &
                            &=  bimAA(mtidst(i)+5,mtidst(j)+5) - bimGGIN(i,j)*nnz(j)

                            END IF

                        END DO

                    END IF

                END DO

!---

            END IF

!End Internal field
!***

        END DO
!$OMP END DO
!$OMP END PARALLEL

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!second argument = 0 -> LU decomposition
        CALL Slv_LinearMatrixCOMP_LU(MatSize, 0)

!$OMP PARALLEL PRIVATE (i,j,ii,jj,kk,NdA) &
!$OMP & PRIVATE (ithprtl,jthprtl,id_tp) &
!$OMP & PRIVATE (tp,tp1,tp2,tp3,ztp,ztp1,ztp2,ztp3,mdm_oi,mdm_io)
!$OMP DO

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

            IF (BCType_EM(ithprtl) == 'PEC') THEN

                ztp1 = exH1x_EM(i)*nnx(i)&
                &     +exH1y_EM(i)*nny(i)&
                &     +exH1z_EM(i)*nnz(i)
                ztp1 =-ztp1

                exH2x_EM(i) =              ztp1*nnx(i) &
                &            +bimX(mtidst(i)  )*t1x(i) &
                &            +bimX(mtidst(i)+1)*t2x(i)
                exH2y_EM(i) =              ztp1*nny(i) &
                &            +bimX(mtidst(i)  )*t1y(i) &
                &            +bimX(mtidst(i)+1)*t2y(i)
                exH2z_EM(i) =              ztp1*nnz(i) &
                &            +bimX(mtidst(i)  )*t1z(i) &
                &            +bimX(mtidst(i)+1)*t2z(i)

                exH3x_EM(i) = exH1x_EM(i) + exH2x_EM(i)
                exH3y_EM(i) = exH1y_EM(i) + exH2y_EM(i)
                exH3z_EM(i) = exH1z_EM(i) + exH2z_EM(i)

                ztp1 = ztpzero
                ztp2 = ztpzero
                DO kk = 1, mxnmbrndlnknd2ndslf
                    IF (ndlnknd2ndslf(i,kk) /= 0) THEN

                        NdA = ndlnknd2ndslf(i,kk)
                        ii = bcscl*(kk-1)

                        ztp1 = ztp1 + bimX(mtidst(NdA)  ) * TanAAt1(i,ii+1) &
                        &             + bimX(mtidst(NdA)+1) * TanAAt1(i,ii+2)
                        ztp2 = ztp2 + bimX(mtidst(NdA)  ) * TanAAt2(i,ii+1) &
                        &             + bimX(mtidst(NdA)+1) * TanAAt2(i,ii+2)

                    END IF
                END DO

                ztp1 = ztp1 + TanBBt1(i)
                ztp2 = ztp2 + TanBBt2(i)

                exH2xdnn_EM(i) =  bimX(mtidst(i)+2)*nnx(i) &
                &                +             ztp1*t1x(i) &
                &                +             ztp2*t2x(i)
                exH2ydnn_EM(i) =  bimX(mtidst(i)+2)*nny(i) &
                &                +             ztp1*t1y(i) &
                &                +             ztp2*t2y(i)
                exH2zdnn_EM(i) =  bimX(mtidst(i)+2)*nnz(i) &
                &                +             ztp1*t1z(i) &
                &                +             ztp2*t2z(i)

                exH3xdnn_EM(i) = exH1xdnn_EM(i) + exH2xdnn_EM(i)
                exH3ydnn_EM(i) = exH1ydnn_EM(i) + exH2ydnn_EM(i)
                exH3zdnn_EM(i) = exH1zdnn_EM(i) + exH2zdnn_EM(i)

            END IF

            IF (BCType_EM(ithprtl) == '2SD') THEN

                id_tp = corelnkshell(ithprtl)
                IF (id_tp == 0) THEN
                    mdm_oi = exmiu_EM/inmiu_EM(ithprtl)
                ELSE
                    mdm_oi = inmiu_EM(id_tp)/inmiu_EM(ithprtl)
                END IF

                exH2x_EM(i) = bimX(mtidst(i)  )*t1x(i) &
                &            +bimX(mtidst(i)+1)*t2x(i) &
                &            +bimX(mtidst(i)+2)*nnx(i)
                exH2y_EM(i) = bimX(mtidst(i)  )*t1y(i) &
                &            +bimX(mtidst(i)+1)*t2y(i) &
                &            +bimX(mtidst(i)+2)*nny(i)
                exH2z_EM(i) = bimX(mtidst(i)  )*t1z(i) &
                &            +bimX(mtidst(i)+1)*t2z(i) &
                &            +bimX(mtidst(i)+2)*nnz(i)

                exH3x_EM(i) = exH1x_EM(i) + exH2x_EM(i)
                exH3y_EM(i) = exH1y_EM(i) + exH2y_EM(i)
                exH3z_EM(i) = exH1z_EM(i) + exH2z_EM(i)

                exH2xdnn_EM(i) = bimX(mtidst(i)+3)*t1x(i) &
                &               +bimX(mtidst(i)+4)*t2x(i) &
                &               +bimX(mtidst(i)+5)*nnx(i)
                exH2ydnn_EM(i) = bimX(mtidst(i)+3)*t1y(i) &
                &               +bimX(mtidst(i)+4)*t2y(i) &
                &               +bimX(mtidst(i)+5)*nny(i)
                exH2zdnn_EM(i) = bimX(mtidst(i)+3)*t1z(i) &
                &               +bimX(mtidst(i)+4)*t2z(i) &
                &               +bimX(mtidst(i)+5)*nnz(i)

                exH3xdnn_EM(i) = exH1xdnn_EM(i) + exH2xdnn_EM(i)
                exH3ydnn_EM(i) = exH1ydnn_EM(i) + exH2ydnn_EM(i)
                exH3zdnn_EM(i) = exH1zdnn_EM(i) + exH2zdnn_EM(i)

                ztp1 = (exH3x_EM(i)       -inH1x_EM(i))*t1x(i)&
                &     +(exH3y_EM(i)       -inH1y_EM(i))*t1y(i)&
                &     +(exH3z_EM(i)       -inH1z_EM(i))*t1z(i)
                ztp2 = (exH3x_EM(i)       -inH1x_EM(i))*t2x(i)&
                &     +(exH3y_EM(i)       -inH1y_EM(i))*t2y(i)&
                &     +(exH3z_EM(i)       -inH1z_EM(i))*t2z(i)
                ztp3 = (exH3x_EM(i)*mdm_oi-inH1x_EM(i))*nnx(i)&
                &     +(exH3y_EM(i)*mdm_oi-inH1y_EM(i))*nny(i)&
                &     +(exH3z_EM(i)*mdm_oi-inH1z_EM(i))*nnz(i)

                inH2x_EM(i) = ztp1*t1x(i) &
                &            +ztp2*t2x(i) &
                &            +ztp3*nnx(i)
                inH2y_EM(i) = ztp1*t1y(i) &
                &            +ztp2*t2y(i) &
                &            +ztp3*nny(i)
                inH2z_EM(i) = ztp1*t1z(i) &
                &            +ztp2*t2z(i) &
                &            +ztp3*nnz(i)

                inH3x_EM(i) = inH1x_EM(i) + inH2x_EM(i)
                inH3y_EM(i) = inH1y_EM(i) + inH2y_EM(i)
                inH3z_EM(i) = inH1z_EM(i) + inH2z_EM(i)

                ztp1 = ztpzero
                ztp2 = ztpzero
                ztp3 = ztpzero
                DO kk = 1, mxnmbrndlnknd2ndslf
                    IF (ndlnknd2ndslf(i,kk) /= 0) THEN

                        NdA = ndlnknd2ndslf(i,kk)
                        ii = 6*(kk-1)

                        DO id_tp = 1,  6

                            ztp1 = ztp1 + bimX(mtidst(NdA)+id_tp-1) * TanAAt1(i,ii+id_tp)
                            ztp2 = ztp2 + bimX(mtidst(NdA)+id_tp-1) * TanAAt2(i,ii+id_tp)
                            ztp3 = ztp3 + bimX(mtidst(NdA)+id_tp-1) * DivFrAA(i,ii+id_tp)

                        END DO

                    END IF
                END DO

                ztp1 = ztp1 + TanBBt1(i)
                ztp2 = ztp2 + TanBBt2(i)
                ztp3 = ztp3 + DivFrBB(i)

                inH2xdnn_EM(i) = ztp1*t1x(i) &
                &               +ztp2*t2x(i) &
                &               +ztp3*nnx(i)
                inH2ydnn_EM(i) = ztp1*t1y(i) &
                &               +ztp2*t2y(i) &
                &               +ztp3*nny(i)
                inH2zdnn_EM(i) = ztp1*t1z(i) &
                &               +ztp2*t2z(i) &
                &               +ztp3*nnz(i)

                inH3xdnn_EM(i) = inH1xdnn_EM(i) + inH2xdnn_EM(i)
                inH3ydnn_EM(i) = inH1ydnn_EM(i) + inH2ydnn_EM(i)
                inH3zdnn_EM(i) = inH1zdnn_EM(i) + inH2zdnn_EM(i)

            END IF

        END DO

!$OMP END DO
!$OMP END PARALLEL

        DEALLOCATE (bimX,TanAAt1,TanAAt2,TanBBt1,TanBBt2,DivFrAA,DivFrBB,mtidst)
        PRINT *, 'H OK by Hsc'

        IF (PostProcdn_on == 1) THEN

            DO icmpt = 1, 3

                Ttlbim = ttlnmbrnd
                ALLOCATE (mtidst(ttlnmbrnd))

                MatSize = 0
                DO ithprtl = 1, nmbrprtl
                    IF (BCType_EM(ithprtl) == 'PEC') THEN
                        id_tp = 1
                        DO i = ndstaID(ithprtl), ndendID(ithprtl)
                            mtidst(i) = MatSize + id_tp
                            id_tp = id_tp + 1
                        END DO
                        MatSize = MatSize + 1*nmbrnd(ithprtl)
                    END IF
                    IF (BCType_EM(ithprtl) == '2SD') THEN
                        id_tp = 1
                        DO i = ndstaID(ithprtl), ndendID(ithprtl)
                            mtidst(i) = MatSize + id_tp
                            id_tp = id_tp + 2
                        END DO
                        MatSize = MatSize + 2*nmbrnd(ithprtl)
                    END IF
                END DO

                ALLOCATE(bimAA(MatSize,MatSize))
                ALLOCATE(bimBB(MatSize))
                ALLOCATE(bimX(MatSize))

!$OMP PARALLEL PRIVATE (i,j)
!$OMP DO
                DO i = 1, MatSize
                    DO j = 1, MatSize
                        bimAA(i,j) = ztpzero
                    END DO
                    bimBB(i) = ztpzero
                    bimX(i) = ztpzero
                END DO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE (i,j,ii,jj,kk,NdA) &
!$OMP & PRIVATE (ithprtl,jthprtl,id_tp) &
!$OMP & PRIVATE (tp,tp1,tp2,tp3,ztp,ztp1,ztp2,ztp3,mdm_oi,mdm_io)
!$OMP DO
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

!***
! Exterior-field equations.
                    IF (     BCType_EM(ithprtl) == 'PEC' &
                    &   .OR. BCType_EM(ithprtl) == '2SD') THEN

!---
! Helmholtz equations for E, rows 1--3 of the assembled system.
!Part 1: surfaces share the same 'corelnkshell'
!!!NOTE: 'E' is 'exE2'

                        DO jthprtl = 1, nmbrprtl

                            id_tp = corelnkshell(jthprtl)

                            IF (corelnkshell(jthprtl) == corelnkshell(ithprtl)) THEN

                                DO j = ndstaID(jthprtl), ndendID(jthprtl)

                                    IF (     BCType_EM(jthprtl) == 'PEC' &
                                    &   .OR. BCType_EM(jthprtl) == '2SD') THEN

                                        !∂(exH2x)/∂n,∂(exH2y)/∂n,∂(exH2z)/∂n
                                        bimAA(mtidst(i)  ,mtidst(j)  ) &
                                    &=  bimAA(mtidst(i)  ,mtidst(j)  ) + bimGGEX(i,j)

                                        IF (icmpt == 1) &
                                        &   bimBB(mtidst(i)  ) = bimBB(mtidst(i)  ) + bimHHEX(i,j)*exH2x_EM(j)
                                        IF (icmpt == 2) &
                                        &   bimBB(mtidst(i)  ) = bimBB(mtidst(i)  ) + bimHHEX(i,j)*exH2y_EM(j)
                                        IF (icmpt == 3) &
                                        &   bimBB(mtidst(i)  ) = bimBB(mtidst(i)  ) + bimHHEX(i,j)*exH2z_EM(j)

                                    END IF



                                END DO

                            END IF

                        END DO

!---

!---
! Continue the E rows 1--3 for bounded exterior-domain surfaces.
!!!NOTE: on this bounded surface, 'E' is 'inE2'
!!!NOTE: ROW 4 to 6 with GGEX and HHEX in the big Matrix (Table 1) in Die_MS II
!!!NOTE: Eq IDs are still 1 to 3, such as mtidst(i), mtidst(i)+1, mtidst(i)+2

                        IF (corelnkshell(ithprtl) > 0) THEN

                            jthprtl = corelnkshell(ithprtl)

                            id_tp = corelnkshell(jthprtl)
                            IF (id_tp == 0) THEN
                                mdm_oi = exeps_EM/ineps_EM(jthprtl)
                            ELSE
                                mdm_oi = ineps_EM(id_tp)/ineps_EM(jthprtl)
                            END IF

                            DO j = ndstaID(jthprtl), ndendID(jthprtl)

                                IF (BCType_EM(jthprtl) == '2SD') THEN

                                    !∂(exH2x)/∂n,∂(exH2y)/∂n,∂(exH2z)/∂n
                                    bimAA(mtidst(i)  ,mtidst(j)  ) &
                                &=  bimAA(mtidst(i)  ,mtidst(j)  ) + bimGGEX(i,j)

                                    IF (icmpt == 1) &
                                    &   bimBB(mtidst(i)  ) = bimBB(mtidst(i)  ) + bimHHEX(i,j)*inH2x_EM(j)
                                    IF (icmpt == 2) &
                                    &   bimBB(mtidst(i)  ) = bimBB(mtidst(i)  ) + bimHHEX(i,j)*inH2y_EM(j)
                                    IF (icmpt == 3) &
                                    &   bimBB(mtidst(i)  ) = bimBB(mtidst(i)  ) + bimHHEX(i,j)*inH2z_EM(j)

                                END IF

                            END DO

                        END IF

!---

                    END IF

!End External field
!***

!***
! Interior-field equations.
                    IF (BCType_EM(ithprtl) == '2SD') THEN

!---
! Helmholtz equations for E, rows 4--6 of the assembled system.

                        DO jthprtl = 1, nmbrprtl

                            IF (jthprtl == ithprtl) THEN

                                id_tp = corelnkshell(jthprtl)
                                IF (id_tp == 0) THEN
                                    mdm_oi = exeps_EM/ineps_EM(jthprtl)
                                ELSE
                                    mdm_oi = ineps_EM(id_tp)/ineps_EM(jthprtl)
                                END IF

                                DO j = ndstaID(jthprtl), ndendID(jthprtl)

                                    !∂(inH2x)/∂n,∂(inH2y)/∂n,∂(inH2z)/∂n
                                    bimAA(mtidst(i)+1,mtidst(j)+1) &
                                &=  bimAA(mtidst(i)+1,mtidst(j)+1) + bimGGIN(i,j)

                                    IF (icmpt == 1) &
                                    &   bimBB(mtidst(i)+1) = bimBB(mtidst(i)+1) + bimHHIN(i,j)*inH2x_EM(j)
                                    IF (icmpt == 2) &
                                    &   bimBB(mtidst(i)+1) = bimBB(mtidst(i)+1) + bimHHIN(i,j)*inH2y_EM(j)
                                    IF (icmpt == 3) &
                                    &   bimBB(mtidst(i)+1) = bimBB(mtidst(i)+1) + bimHHIN(i,j)*inH2z_EM(j)

                                END DO

                            END IF

!---

!---
!Part2: the surfaces that enclosed by the surface where x0 locates
!(just one level down needed)
!!NOTE: 'E' is 'exE2'
!!NOTE: ROW 1 to 3 with GGIN and HHIN in the big Matrix (Table 1) in Die_MS II
!!NOTE: Eq IDs are still 4 to 6, such as mtidst(i)+3, mtidst(i)+4, mtidst(i)+5

                            IF (corelnkshell(jthprtl) == ithprtl) THEN

                                DO j = ndstaID(jthprtl), ndendID(jthprtl)

                                    IF (     BCType_EM(jthprtl) == 'PEC' &
                                    &   .OR. BCType_EM(jthprtl) == '2SD') THEN

                                        !∂(inH2x)/∂n,∂(inH2y)/∂n,∂(inH2z)/∂n
                                        bimAA(mtidst(i)+1,mtidst(j)+1) &
                                    &=  bimAA(mtidst(i)+1,mtidst(j)+1) + bimGGIN(i,j)

                                        IF (icmpt == 1) &
                                        &   bimBB(mtidst(i)+1) = bimBB(mtidst(i)+1) + bimHHIN(i,j)*exH2x_EM(j)
                                        IF (icmpt == 2) &
                                        &   bimBB(mtidst(i)+1) = bimBB(mtidst(i)+1) + bimHHIN(i,j)*exH2y_EM(j)
                                        IF (icmpt == 3) &
                                        &   bimBB(mtidst(i)+1) = bimBB(mtidst(i)+1) + bimHHIN(i,j)*exH2z_EM(j)

                                    END IF

                                END DO

                            END IF

                        END DO

!---

                    END IF

!End Internal field
!***

                END DO
!$OMP END DO
!$OMP END PARALLEL

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!second argument = 0 -> LU decomposition
                CALL Slv_LinearMatrixCOMP_LU(MatSize, 0)

!$OMP PARALLEL PRIVATE (i,j,ii,jj,kk,NdA) &
!$OMP & PRIVATE (ithprtl,jthprtl,id_tp) &
!$OMP & PRIVATE (tp,tp1,tp2,tp3,ztp,ztp1,ztp2,ztp3,mdm_oi,mdm_io)
!$OMP DO

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

                    IF (     BCType_EM(ithprtl) == 'PEC' &
                    &   .OR. BCType_EM(ithprtl) == '2SD') THEN

                        IF (icmpt == 1) THEN
                            exH2xdnn_EM(i) = bimX(mtidst(i)  )
                            exH3xdnn_EM(i) = exH1xdnn_EM(i) + exH2xdnn_EM(i)
                        END IF
                        IF (icmpt == 2) THEN
                            exH2ydnn_EM(i) = bimX(mtidst(i)  )
                            exH3ydnn_EM(i) = exH1ydnn_EM(i) + exH2ydnn_EM(i)
                        END IF
                        IF (icmpt == 3) THEN
                            exH2zdnn_EM(i) = bimX(mtidst(i)  )
                            exH3zdnn_EM(i) = exH1zdnn_EM(i) + exH2zdnn_EM(i)
                        END IF

                    END IF

                    IF (BCType_EM(ithprtl) == '2SD') THEN

                        IF (icmpt == 1) THEN
                            inH2xdnn_EM(i) = bimX(mtidst(i)+1)
                            inH3xdnn_EM(i) = inH1xdnn_EM(i) + inH2xdnn_EM(i)
                        END IF
                        IF (icmpt == 2) THEN
                            inH2ydnn_EM(i) = bimX(mtidst(i)+1)
                            inH3ydnn_EM(i) = inH1ydnn_EM(i) + inH2ydnn_EM(i)
                        END IF
                        IF (icmpt == 3) THEN
                            inH2zdnn_EM(i) = bimX(mtidst(i)+1)
                            inH3zdnn_EM(i) = inH1zdnn_EM(i) + inH2zdnn_EM(i)
                        END IF

                    END IF

                END DO

!$OMP END DO
!$OMP END PARALLEL

                DEALLOCATE (bimX,mtidst)
                PRINT *, 'Hdn OK by post-processing'

            END DO

        END IF

        DEALLOCATE (bimGGEX, bimHHEX, bimGGIN, bimHHIN)

    END SUBROUTINE





    ! Dispatch the electric and magnetic boundary solves required by the
    ! selected excitation type. Released plane-wave cases execute both.

    SUBROUTINE SlvPrblm_EM

        IF (excitetype_EM == 'spe') THEN
            CALL SlvPrblm_EM_Esc
        ELSE IF (excitetype_EM == 'sph') THEN
            CALL SlvPrblm_EM_Hsc
        ELSE
            CALL SlvPrblm_EM_Esc
            CALL SlvPrblm_EM_Hsc
        END IF

    END SUBROUTINE




    ! Solve bimAA*bimX=bimBB with LAPACK ZGESV.
    ! TtlbimA is the dense system order. CalMtrxCNDet_On is retained for
    ! call-site compatibility and is not used by the release direct-solve path.

    SUBROUTINE Slv_LinearMatrixCOMP_LU(TtlbimA,CalMtrxCNDet_On)

        INTEGER, INTENT(IN) :: TtlbimA,CalMtrxCNDet_On
        INTEGER :: i, slv_info
        INTEGER, ALLOCATABLE :: slv_ipiv(:)

        ! The manuscript examples use the LAPACK direct-solve path.
        ! Legacy custom LU, iterative BiCGStab and preconditioner branches
        ! were unreachable in these drivers and are intentionally excluded.
        ALLOCATE(slv_ipiv(TtlbimA))
        CALL zgesv(TtlbimA,1,bimAA,TtlbimA,slv_ipiv,bimBB,TtlbimA,slv_info)
        IF (slv_info /= 0) THEN
            WRITE(*,*) "LAPACK zgesv failed with INFO = ",slv_info
            STOP 91
        END IF
        DO i=1,TtlbimA
            bimX(i)=bimBB(i)
        END DO
        DEALLOCATE(slv_ipiv,bimAA,bimBB)

    END SUBROUTINE Slv_LinearMatrixCOMP_LU

END MODULE
