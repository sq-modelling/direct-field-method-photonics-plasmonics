! SPDX-FileCopyrightText: 2026 Qiang Sun
! SPDX-License-Identifier: BSD-3-Clause
!
! Surface integration of time-averaged electromagnetic force and torque.

MODULE EM_SurfCal_ForceTorque

    USE omp_lib

    USE Pre_Constants
    USE Pre_csvformat

    USE Geom_GlobalData

    USE EM_SurfCal_GlobalData

    IMPLICIT NONE

    CONTAINS

    ! Contract: mesh coordinates, normals, quadrature weights, centres of mass,
    ! material parameters and solved exterior incident/scattered/total phasors
    ! must be initialised in the imported global modules.  For every outermost
    ! particle (corelnkshell=0), integrate traction on linear or Q6 elements and
    ! populate the direct, electric/magnetic and incident/scattered/interference
    ! force and torque arrays.  Embedded interfaces are returned as zero.
    !
    ! Time averaging uses <a b> = 0.5*Re(a*conjg(b)).  Define D_r=epsilon_r*E
    ! and B_r=mu_r*H (not the SI D and B); the implemented symmetrised stress
    ! uses D_r,E and B_r,H with epsilon_0 and mu_0 prefactors.  NrmlInOut=1
    ! points toward the supplied reference point.  For the released closed
    ! particles that point is interior, so the traction sign is reversed for
    ! force on the body.  Torque is integrated as (r-r_cm) cross traction.
    ! This routine applies no output normalisation; the calling driver selects
    ! any reporting scale.
    !
    ! Background: P. W. Milonni and R. W. Boyd, "Momentum of light in a
    ! dielectric medium," Adv. Opt. Photonics 2, 519 (2010); D. Gao et al.,
    ! "Optical manipulation from the microscale to the nanoscale," Light Sci.
    ! Appl. 6, e17039.
    SUBROUTINE GetFrcTrq_EM

        DOUBLE PRECISION :: tp, tp1, tp2, tp3, tp4, tp5, tp6, tprx, tpry, tprz, tpkappa0

        COMPLEX(KIND=KIND(1.0D0)) :: ztp, ztp1, ztp2, ztp3, ztp4, ztp5, ztp6

        DOUBLE PRECISION :: tpnx, tpny, tpnz

        INTEGER :: ithprtl, GLQi, k, i, icnt

        INTEGER :: ndoffset1, elmntoffset1

        DOUBLE PRECISION, DIMENSION (ttlnmbrnd) :: tpndFrcx_EM, tpndFrcy_EM, tpndFrcz_EM

        DOUBLE PRECISION, DIMENSION (ttlnmbrnd) :: tpndTrqx_EM, tpndTrqy_EM, tpndTrqz_EM

        ! Direct integral of the complete exterior total field.
        DO ithprtl = 1, nmbrprtl

            Frcx_EM(ithprtl) = 0.0d0
            Frcy_EM(ithprtl) = 0.0d0
            Frcz_EM(ithprtl) = 0.0d0

            Trqx_EM(ithprtl) = 0.0d0
            Trqy_EM(ithprtl) = 0.0d0
            Trqz_EM(ithprtl) = 0.0d0

            IF (corelnkshell(ithprtl) == 0) THEN

                DO i = ndstaID(ithprtl), ndendID(ithprtl)

                    tpnx = nnx(i)
                    tpny = nny(i)
                    tpnz = nnz(i)

                    tp4 = 0.0d0
                    tp5 = 0.0d0
                    tp6 = 0.0d0

                    tp = 0.5d0*REAL( exeps_EM*exE3x_EM(i) * DCONJG(exE3x_EM(i)) &
                    &               +exeps_EM*exE3y_EM(i) * DCONJG(exE3y_EM(i)) &
                    &               +exeps_EM*exE3z_EM(i) * DCONJG(exE3z_EM(i))  )
                    tp4 = tp4 - tp*tpnx
                    tp5 = tp5 - tp*tpny
                    tp6 = tp6 - tp*tpnz

                    ztp1 = exE3x_EM(i)*tpnx &
                    &     +exE3y_EM(i)*tpny &
                    &     +exE3z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(exeps_EM*ztp1*DCONJG(exE3x_EM(i)))
                    tp2 = 0.5d0*REAL(exeps_EM*ztp1*DCONJG(exE3y_EM(i)))
                    tp3 = 0.5d0*REAL(exeps_EM*ztp1*DCONJG(exE3z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3

                    ztp1 = exE3x_EM(i)*tpnx &
                    &     +exE3y_EM(i)*tpny &
                    &     +exE3z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(ztp1*DCONJG(exeps_EM*exE3x_EM(i)))
                    tp2 = 0.5d0*REAL(ztp1*DCONJG(exeps_EM*exE3y_EM(i)))
                    tp3 = 0.5d0*REAL(ztp1*DCONJG(exeps_EM*exE3z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3

                    tpndFrcx_EM(i) = 0.5d0*vcm_eps0*tp4
                    tpndFrcy_EM(i) = 0.5d0*vcm_eps0*tp5
                    tpndFrcz_EM(i) = 0.5d0*vcm_eps0*tp6

                    tp4 = 0.0d0
                    tp5 = 0.0d0
                    tp6 = 0.0d0

                    tp = 0.5d0*REAL( exmiu_EM*exH3x_EM(i) * DCONJG(exH3x_EM(i)) &
                    &               +exmiu_EM*exH3y_EM(i) * DCONJG(exH3y_EM(i)) &
                    &               +exmiu_EM*exH3z_EM(i) * DCONJG(exH3z_EM(i))  )
                    tp4 = tp4 - tp*tpnx
                    tp5 = tp5 - tp*tpny
                    tp6 = tp6 - tp*tpnz

                    ztp1 = exH3x_EM(i)*tpnx &
                    &     +exH3y_EM(i)*tpny &
                    &     +exH3z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(exmiu_EM*ztp1*DCONJG(exH3x_EM(i)))
                    tp2 = 0.5d0*REAL(exmiu_EM*ztp1*DCONJG(exH3y_EM(i)))
                    tp3 = 0.5d0*REAL(exmiu_EM*ztp1*DCONJG(exH3z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3

                    ztp1 = exH3x_EM(i)*tpnx &
                    &     +exH3y_EM(i)*tpny &
                    &     +exH3z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(ztp1*DCONJG(exmiu_EM*exH3x_EM(i)))
                    tp2 = 0.5d0*REAL(ztp1*DCONJG(exmiu_EM*exH3y_EM(i)))
                    tp3 = 0.5d0*REAL(ztp1*DCONJG(exmiu_EM*exH3z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3

                    tpndFrcx_EM(i) = tpndFrcx_EM(i) + 0.5d0*vcm_mu0*tp4
                    tpndFrcy_EM(i) = tpndFrcy_EM(i) + 0.5d0*vcm_mu0*tp5
                    tpndFrcz_EM(i) = tpndFrcz_EM(i) + 0.5d0*vcm_mu0*tp6

                    ! Convert an inward code normal to the body-traction sign.
                    IF (NrmlInOut(ithprtl) == 1) THEN
                        tpndFrcx_EM(i) =-tpndFrcx_EM(i)
                        tpndFrcy_EM(i) =-tpndFrcy_EM(i)
                        tpndFrcz_EM(i) =-tpndFrcz_EM(i)
                    END IF

                    ! Form torque density about the particle centre of mass.
                    tprx = xnd(i) - xmssctr(ithprtl)
                    tpry = ynd(i) - ymssctr(ithprtl)
                    tprz = znd(i) - zmssctr(ithprtl)
                    tpndTrqx_EM(i) = tpry * tpndFrcz_EM(i) - tprz * tpndFrcy_EM(i)
                    tpndTrqy_EM(i) = tprz * tpndFrcx_EM(i) - tprx * tpndFrcz_EM(i)
                    tpndTrqz_EM(i) = tprx * tpndFrcy_EM(i) - tpry * tpndFrcx_EM(i)

                END DO

                IF (MeshType == 'L') THEN

                    DO k = elstaID(ithprtl), elendID(ithprtl)

                        DO GLQi = 1, n_glqtr2d

                            icnt = n_glqtr2d*(k-1)+GLQi

                            Frcx_EM(ithprtl) = Frcx_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcx_EM(elmntlnknd(k,3))

                            Frcy_EM(ithprtl) = Frcy_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcy_EM(elmntlnknd(k,3))

                            Frcz_EM(ithprtl) = Frcz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcz_EM(elmntlnknd(k,3))

                            Trqx_EM(ithprtl) = Trqx_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqx_EM(elmntlnknd(k,3))

                            Trqy_EM(ithprtl) = Trqy_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqy_EM(elmntlnknd(k,3))

                            Trqz_EM(ithprtl) = Trqz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqz_EM(elmntlnknd(k,3))

                        END DO

                    END DO

                END IF

                IF (MeshType == 'Q') THEN

                    DO k = elstaID(ithprtl), elendID(ithprtl)

                        DO GLQi = 1, n_glqtr2d

                            icnt = n_glqtr2d*(k-1)+GLQi

                            Frcx_EM(ithprtl) = Frcx_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcx_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndFrcx_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndFrcx_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndFrcx_EM(elmntlnknd(k,6))

                            Frcy_EM(ithprtl) = Frcy_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcy_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndFrcy_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndFrcy_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndFrcy_EM(elmntlnknd(k,6))

                            Frcz_EM(ithprtl) = Frcz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcz_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndFrcz_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndFrcz_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndFrcz_EM(elmntlnknd(k,6))

                            Trqx_EM(ithprtl) = Trqx_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqx_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndTrqx_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndTrqx_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndTrqx_EM(elmntlnknd(k,6))

                            Trqy_EM(ithprtl) = Trqy_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqy_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndTrqy_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndTrqy_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndTrqy_EM(elmntlnknd(k,6))

                            Trqz_EM(ithprtl) = Trqz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqz_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndTrqz_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndTrqz_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndTrqz_EM(elmntlnknd(k,6))


                        END DO

                    END DO


                END IF

            END IF

        END DO



        ! Electric contribution to the complete total-field integral.
        DO ithprtl = 1, nmbrprtl

            FrcElcx_EM(ithprtl) = 0.0d0
            FrcElcy_EM(ithprtl) = 0.0d0
            FrcElcz_EM(ithprtl) = 0.0d0

            TrqElcx_EM(ithprtl) = 0.0d0
            TrqElcy_EM(ithprtl) = 0.0d0
            TrqElcz_EM(ithprtl) = 0.0d0

            IF (corelnkshell(ithprtl) == 0) THEN

                DO i = ndstaID(ithprtl), ndendID(ithprtl)

                    tpnx = nnx(i)
                    tpny = nny(i)
                    tpnz = nnz(i)

                    tp4 = 0.0d0
                    tp5 = 0.0d0
                    tp6 = 0.0d0

                    tp = 0.5d0*REAL( exeps_EM*exE3x_EM(i) * DCONJG(exE3x_EM(i)) &
                    &               +exeps_EM*exE3y_EM(i) * DCONJG(exE3y_EM(i)) &
                    &               +exeps_EM*exE3z_EM(i) * DCONJG(exE3z_EM(i))  )
                    tp4 = tp4 - tp*tpnx
                    tp5 = tp5 - tp*tpny
                    tp6 = tp6 - tp*tpnz

                    ztp1 = exE3x_EM(i)*tpnx &
                    &     +exE3y_EM(i)*tpny &
                    &     +exE3z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(exeps_EM*ztp1*DCONJG(exE3x_EM(i)))
                    tp2 = 0.5d0*REAL(exeps_EM*ztp1*DCONJG(exE3y_EM(i)))
                    tp3 = 0.5d0*REAL(exeps_EM*ztp1*DCONJG(exE3z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3

                    ztp1 = exE3x_EM(i)*tpnx &
                    &     +exE3y_EM(i)*tpny &
                    &     +exE3z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(ztp1*DCONJG(exeps_EM*exE3x_EM(i)))
                    tp2 = 0.5d0*REAL(ztp1*DCONJG(exeps_EM*exE3y_EM(i)))
                    tp3 = 0.5d0*REAL(ztp1*DCONJG(exeps_EM*exE3z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3

                    tpndFrcx_EM(i) = 0.5d0*vcm_eps0*tp4
                    tpndFrcy_EM(i) = 0.5d0*vcm_eps0*tp5
                    tpndFrcz_EM(i) = 0.5d0*vcm_eps0*tp6

                    IF (NrmlInOut(ithprtl) == 1) THEN
                        tpndFrcx_EM(i) =-tpndFrcx_EM(i)
                        tpndFrcy_EM(i) =-tpndFrcy_EM(i)
                        tpndFrcz_EM(i) =-tpndFrcz_EM(i)
                    END IF

                    tprx = xnd(i) - xmssctr(ithprtl)
                    tpry = ynd(i) - ymssctr(ithprtl)
                    tprz = znd(i) - zmssctr(ithprtl)
                    tpndTrqx_EM(i) = tpry * tpndFrcz_EM(i) - tprz * tpndFrcy_EM(i)
                    tpndTrqy_EM(i) = tprz * tpndFrcx_EM(i) - tprx * tpndFrcz_EM(i)
                    tpndTrqz_EM(i) = tprx * tpndFrcy_EM(i) - tpry * tpndFrcx_EM(i)

                END DO

                IF (MeshType == 'L') THEN

                    DO k = elstaID(ithprtl), elendID(ithprtl)

                        DO GLQi = 1, n_glqtr2d

                            icnt = n_glqtr2d*(k-1)+GLQi

                            FrcElcx_EM(ithprtl) = FrcElcx_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcx_EM(elmntlnknd(k,3))

                            FrcElcy_EM(ithprtl) = FrcElcy_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcy_EM(elmntlnknd(k,3))

                            FrcElcz_EM(ithprtl) = FrcElcz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcz_EM(elmntlnknd(k,3))

                            TrqElcx_EM(ithprtl) = TrqElcx_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqx_EM(elmntlnknd(k,3))

                            TrqElcy_EM(ithprtl) = TrqElcy_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqy_EM(elmntlnknd(k,3))

                            TrqElcz_EM(ithprtl) = TrqElcz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqz_EM(elmntlnknd(k,3))

                        END DO

                    END DO

                END IF

                IF (MeshType == 'Q') THEN

                    DO k = elstaID(ithprtl), elendID(ithprtl)

                        DO GLQi = 1, n_glqtr2d

                            icnt = n_glqtr2d*(k-1)+GLQi

                            FrcElcx_EM(ithprtl) = FrcElcx_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcx_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndFrcx_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndFrcx_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndFrcx_EM(elmntlnknd(k,6))

                            FrcElcy_EM(ithprtl) = FrcElcy_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcy_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndFrcy_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndFrcy_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndFrcy_EM(elmntlnknd(k,6))

                            FrcElcz_EM(ithprtl) = FrcElcz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcz_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndFrcz_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndFrcz_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndFrcz_EM(elmntlnknd(k,6))

                            TrqElcx_EM(ithprtl) = TrqElcx_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqx_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndTrqx_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndTrqx_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndTrqx_EM(elmntlnknd(k,6))

                            TrqElcy_EM(ithprtl) = TrqElcy_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqy_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndTrqy_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndTrqy_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndTrqy_EM(elmntlnknd(k,6))

                            TrqElcz_EM(ithprtl) = TrqElcz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqz_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndTrqz_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndTrqz_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndTrqz_EM(elmntlnknd(k,6))


                        END DO

                    END DO


                END IF

            END IF

        END DO


        ! Magnetic contribution to the complete total-field integral.
        DO ithprtl = 1, nmbrprtl

            FrcMagx_EM(ithprtl) = 0.0d0
            FrcMagy_EM(ithprtl) = 0.0d0
            FrcMagz_EM(ithprtl) = 0.0d0

            TrqMagx_EM(ithprtl) = 0.0d0
            TrqMagy_EM(ithprtl) = 0.0d0
            TrqMagz_EM(ithprtl) = 0.0d0

            IF (corelnkshell(ithprtl) == 0) THEN

                DO i = ndstaID(ithprtl), ndendID(ithprtl)

                    tpnx = nnx(i)
                    tpny = nny(i)
                    tpnz = nnz(i)

                    tp4 = 0.0d0
                    tp5 = 0.0d0
                    tp6 = 0.0d0

                    tp = 0.5d0*REAL( exmiu_EM*exH3x_EM(i) * DCONJG(exH3x_EM(i)) &
                    &               +exmiu_EM*exH3y_EM(i) * DCONJG(exH3y_EM(i)) &
                    &               +exmiu_EM*exH3z_EM(i) * DCONJG(exH3z_EM(i))  )
                    tp4 = tp4 - tp*tpnx
                    tp5 = tp5 - tp*tpny
                    tp6 = tp6 - tp*tpnz

                    ztp1 = exH3x_EM(i)*tpnx &
                    &     +exH3y_EM(i)*tpny &
                    &     +exH3z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(exmiu_EM*ztp1*DCONJG(exH3x_EM(i)))
                    tp2 = 0.5d0*REAL(exmiu_EM*ztp1*DCONJG(exH3y_EM(i)))
                    tp3 = 0.5d0*REAL(exmiu_EM*ztp1*DCONJG(exH3z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3

                    ztp1 = exH3x_EM(i)*tpnx &
                    &     +exH3y_EM(i)*tpny &
                    &     +exH3z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(ztp1*DCONJG(exmiu_EM*exH3x_EM(i)))
                    tp2 = 0.5d0*REAL(ztp1*DCONJG(exmiu_EM*exH3y_EM(i)))
                    tp3 = 0.5d0*REAL(ztp1*DCONJG(exmiu_EM*exH3z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3

                    tpndFrcx_EM(i) = 0.5d0*vcm_mu0*tp4
                    tpndFrcy_EM(i) = 0.5d0*vcm_mu0*tp5
                    tpndFrcz_EM(i) = 0.5d0*vcm_mu0*tp6

                    IF (NrmlInOut(ithprtl) == 1) THEN
                        tpndFrcx_EM(i) =-tpndFrcx_EM(i)
                        tpndFrcy_EM(i) =-tpndFrcy_EM(i)
                        tpndFrcz_EM(i) =-tpndFrcz_EM(i)
                    END IF

                    tprx = xnd(i) - xmssctr(ithprtl)
                    tpry = ynd(i) - ymssctr(ithprtl)
                    tprz = znd(i) - zmssctr(ithprtl)
                    tpndTrqx_EM(i) = tpry * tpndFrcz_EM(i) - tprz * tpndFrcy_EM(i)
                    tpndTrqy_EM(i) = tprz * tpndFrcx_EM(i) - tprx * tpndFrcz_EM(i)
                    tpndTrqz_EM(i) = tprx * tpndFrcy_EM(i) - tpry * tpndFrcx_EM(i)

                END DO

                IF (MeshType == 'L') THEN

                    DO k = elstaID(ithprtl), elendID(ithprtl)

                        DO GLQi = 1, n_glqtr2d

                            icnt = n_glqtr2d*(k-1)+GLQi

                            FrcMagx_EM(ithprtl) = FrcMagx_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcx_EM(elmntlnknd(k,3))

                            FrcMagy_EM(ithprtl) = FrcMagy_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcy_EM(elmntlnknd(k,3))

                            FrcMagz_EM(ithprtl) = FrcMagz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcz_EM(elmntlnknd(k,3))

                            TrqMagx_EM(ithprtl) = TrqMagx_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqx_EM(elmntlnknd(k,3))

                            TrqMagy_EM(ithprtl) = TrqMagy_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqy_EM(elmntlnknd(k,3))

                            TrqMagz_EM(ithprtl) = TrqMagz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqz_EM(elmntlnknd(k,3))

                        END DO

                    END DO

                END IF

                IF (MeshType == 'Q') THEN

                    DO k = elstaID(ithprtl), elendID(ithprtl)

                        DO GLQi = 1, n_glqtr2d

                            icnt = n_glqtr2d*(k-1)+GLQi

                            FrcMagx_EM(ithprtl) = FrcMagx_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcx_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndFrcx_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndFrcx_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndFrcx_EM(elmntlnknd(k,6))

                            FrcMagy_EM(ithprtl) = FrcMagy_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcy_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndFrcy_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndFrcy_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndFrcy_EM(elmntlnknd(k,6))

                            FrcMagz_EM(ithprtl) = FrcMagz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcz_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndFrcz_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndFrcz_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndFrcz_EM(elmntlnknd(k,6))

                            TrqMagx_EM(ithprtl) = TrqMagx_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqx_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndTrqx_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndTrqx_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndTrqx_EM(elmntlnknd(k,6))

                            TrqMagy_EM(ithprtl) = TrqMagy_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqy_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndTrqy_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndTrqy_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndTrqy_EM(elmntlnknd(k,6))

                            TrqMagz_EM(ithprtl) = TrqMagz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqz_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndTrqz_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndTrqz_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndTrqz_EM(elmntlnknd(k,6))


                        END DO

                    END DO


                END IF

            END IF

        END DO



        ! Incident-field contribution.
        DO ithprtl = 1, nmbrprtl

            FrcIncx_EM(ithprtl) = 0.0d0
            FrcIncy_EM(ithprtl) = 0.0d0
            FrcIncz_EM(ithprtl) = 0.0d0

            TrqIncx_EM(ithprtl) = 0.0d0
            TrqIncy_EM(ithprtl) = 0.0d0
            TrqIncz_EM(ithprtl) = 0.0d0

            IF (corelnkshell(ithprtl) == 0) THEN

                DO i = ndstaID(ithprtl), ndendID(ithprtl)

                    tpnx = nnx(i)
                    tpny = nny(i)
                    tpnz = nnz(i)

                    tp4 = 0.0d0
                    tp5 = 0.0d0
                    tp6 = 0.0d0

                    tp = 0.5d0*REAL( exeps_EM*exE1x_EM(i) * DCONJG(exE1x_EM(i)) &
                    &               +exeps_EM*exE1y_EM(i) * DCONJG(exE1y_EM(i)) &
                    &               +exeps_EM*exE1z_EM(i) * DCONJG(exE1z_EM(i))  )
                    tp4 = tp4 - tp*tpnx
                    tp5 = tp5 - tp*tpny
                    tp6 = tp6 - tp*tpnz

                    ztp1 = exE1x_EM(i)*tpnx &
                    &     +exE1y_EM(i)*tpny &
                    &     +exE1z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(exeps_EM*ztp1*DCONJG(exE1x_EM(i)))
                    tp2 = 0.5d0*REAL(exeps_EM*ztp1*DCONJG(exE1y_EM(i)))
                    tp3 = 0.5d0*REAL(exeps_EM*ztp1*DCONJG(exE1z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3

                    ztp1 = exE1x_EM(i)*tpnx &
                    &     +exE1y_EM(i)*tpny &
                    &     +exE1z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(ztp1*DCONJG(exeps_EM*exE1x_EM(i)))
                    tp2 = 0.5d0*REAL(ztp1*DCONJG(exeps_EM*exE1y_EM(i)))
                    tp3 = 0.5d0*REAL(ztp1*DCONJG(exeps_EM*exE1z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3

                    tpndFrcx_EM(i) = 0.5d0*vcm_eps0*tp4
                    tpndFrcy_EM(i) = 0.5d0*vcm_eps0*tp5
                    tpndFrcz_EM(i) = 0.5d0*vcm_eps0*tp6

                    tp4 = 0.0d0
                    tp5 = 0.0d0
                    tp6 = 0.0d0

                    tp = 0.5d0*REAL( exmiu_EM*exH1x_EM(i) * DCONJG(exH1x_EM(i)) &
                    &               +exmiu_EM*exH1y_EM(i) * DCONJG(exH1y_EM(i)) &
                    &               +exmiu_EM*exH1z_EM(i) * DCONJG(exH1z_EM(i))  )
                    tp4 = tp4 - tp*tpnx
                    tp5 = tp5 - tp*tpny
                    tp6 = tp6 - tp*tpnz

                    ztp1 = exH1x_EM(i)*tpnx &
                    &     +exH1y_EM(i)*tpny &
                    &     +exH1z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(exmiu_EM*ztp1*DCONJG(exH1x_EM(i)))
                    tp2 = 0.5d0*REAL(exmiu_EM*ztp1*DCONJG(exH1y_EM(i)))
                    tp3 = 0.5d0*REAL(exmiu_EM*ztp1*DCONJG(exH1z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3

                    ztp1 = exH1x_EM(i)*tpnx &
                    &     +exH1y_EM(i)*tpny &
                    &     +exH1z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(ztp1*DCONJG(exmiu_EM*exH1x_EM(i)))
                    tp2 = 0.5d0*REAL(ztp1*DCONJG(exmiu_EM*exH1y_EM(i)))
                    tp3 = 0.5d0*REAL(ztp1*DCONJG(exmiu_EM*exH1z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3

                    tpndFrcx_EM(i) = tpndFrcx_EM(i) + 0.5d0*vcm_mu0*tp4
                    tpndFrcy_EM(i) = tpndFrcy_EM(i) + 0.5d0*vcm_mu0*tp5
                    tpndFrcz_EM(i) = tpndFrcz_EM(i) + 0.5d0*vcm_mu0*tp6

                    IF (NrmlInOut(ithprtl) == 1) THEN
                        tpndFrcx_EM(i) =-tpndFrcx_EM(i)
                        tpndFrcy_EM(i) =-tpndFrcy_EM(i)
                        tpndFrcz_EM(i) =-tpndFrcz_EM(i)
                    END IF

                    tprx = xnd(i) - xmssctr(ithprtl)
                    tpry = ynd(i) - ymssctr(ithprtl)
                    tprz = znd(i) - zmssctr(ithprtl)
                    tpndTrqx_EM(i) = tpry * tpndFrcz_EM(i) - tprz * tpndFrcy_EM(i)
                    tpndTrqy_EM(i) = tprz * tpndFrcx_EM(i) - tprx * tpndFrcz_EM(i)
                    tpndTrqz_EM(i) = tprx * tpndFrcy_EM(i) - tpry * tpndFrcx_EM(i)

                END DO

                IF (MeshType == 'L') THEN

                    DO k = elstaID(ithprtl), elendID(ithprtl)

                        DO GLQi = 1, n_glqtr2d

                            icnt = n_glqtr2d*(k-1)+GLQi

                            FrcIncx_EM(ithprtl) = FrcIncx_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcx_EM(elmntlnknd(k,3))

                            FrcIncy_EM(ithprtl) = FrcIncy_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcy_EM(elmntlnknd(k,3))

                            FrcIncz_EM(ithprtl) = FrcIncz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcz_EM(elmntlnknd(k,3))

                            TrqIncx_EM(ithprtl) = TrqIncx_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqx_EM(elmntlnknd(k,3))

                            TrqIncy_EM(ithprtl) = TrqIncy_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqy_EM(elmntlnknd(k,3))

                            TrqIncz_EM(ithprtl) = TrqIncz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqz_EM(elmntlnknd(k,3))

                        END DO

                    END DO

                END IF

                IF (MeshType == 'Q') THEN

                    DO k = elstaID(ithprtl), elendID(ithprtl)

                        DO GLQi = 1, n_glqtr2d

                            icnt = n_glqtr2d*(k-1)+GLQi

                            FrcIncx_EM(ithprtl) = FrcIncx_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcx_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndFrcx_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndFrcx_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndFrcx_EM(elmntlnknd(k,6))

                            FrcIncy_EM(ithprtl) = FrcIncy_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcy_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndFrcy_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndFrcy_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndFrcy_EM(elmntlnknd(k,6))

                            FrcIncz_EM(ithprtl) = FrcIncz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcz_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndFrcz_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndFrcz_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndFrcz_EM(elmntlnknd(k,6))

                            TrqIncx_EM(ithprtl) = TrqIncx_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqx_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndTrqx_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndTrqx_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndTrqx_EM(elmntlnknd(k,6))

                            TrqIncy_EM(ithprtl) = TrqIncy_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqy_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndTrqy_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndTrqy_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndTrqy_EM(elmntlnknd(k,6))

                            TrqIncz_EM(ithprtl) = TrqIncz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqz_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndTrqz_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndTrqz_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndTrqz_EM(elmntlnknd(k,6))


                        END DO

                    END DO


                END IF

            END IF

        END DO




        ! Scattered-field contribution.
        DO ithprtl = 1, nmbrprtl

            FrcScax_EM(ithprtl) = 0.0d0
            FrcScay_EM(ithprtl) = 0.0d0
            FrcScaz_EM(ithprtl) = 0.0d0

            TrqScax_EM(ithprtl) = 0.0d0
            TrqScay_EM(ithprtl) = 0.0d0
            TrqScaz_EM(ithprtl) = 0.0d0

            IF (corelnkshell(ithprtl) == 0) THEN

                DO i = ndstaID(ithprtl), ndendID(ithprtl)

                    tpnx = nnx(i)
                    tpny = nny(i)
                    tpnz = nnz(i)

                    tp4 = 0.0d0
                    tp5 = 0.0d0
                    tp6 = 0.0d0

                    tp = 0.5d0*REAL( exeps_EM*exE2x_EM(i) * DCONJG(exE2x_EM(i)) &
                    &               +exeps_EM*exE2y_EM(i) * DCONJG(exE2y_EM(i)) &
                    &               +exeps_EM*exE2z_EM(i) * DCONJG(exE2z_EM(i))  )
                    tp4 = tp4 - tp*tpnx
                    tp5 = tp5 - tp*tpny
                    tp6 = tp6 - tp*tpnz

                    ztp1 = exE2x_EM(i)*tpnx &
                    &     +exE2y_EM(i)*tpny &
                    &     +exE2z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(exeps_EM*ztp1*DCONJG(exE2x_EM(i)))
                    tp2 = 0.5d0*REAL(exeps_EM*ztp1*DCONJG(exE2y_EM(i)))
                    tp3 = 0.5d0*REAL(exeps_EM*ztp1*DCONJG(exE2z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3

                    ztp1 = exE2x_EM(i)*tpnx &
                    &     +exE2y_EM(i)*tpny &
                    &     +exE2z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(ztp1*DCONJG(exeps_EM*exE2x_EM(i)))
                    tp2 = 0.5d0*REAL(ztp1*DCONJG(exeps_EM*exE2y_EM(i)))
                    tp3 = 0.5d0*REAL(ztp1*DCONJG(exeps_EM*exE2z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3

                    tpndFrcx_EM(i) = 0.5d0*vcm_eps0*tp4
                    tpndFrcy_EM(i) = 0.5d0*vcm_eps0*tp5
                    tpndFrcz_EM(i) = 0.5d0*vcm_eps0*tp6

                    tp4 = 0.0d0
                    tp5 = 0.0d0
                    tp6 = 0.0d0

                    tp = 0.5d0*REAL( exmiu_EM*exH2x_EM(i) * DCONJG(exH2x_EM(i)) &
                    &               +exmiu_EM*exH2y_EM(i) * DCONJG(exH2y_EM(i)) &
                    &               +exmiu_EM*exH2z_EM(i) * DCONJG(exH2z_EM(i))  )
                    tp4 = tp4 - tp*tpnx
                    tp5 = tp5 - tp*tpny
                    tp6 = tp6 - tp*tpnz

                    ztp1 = exH2x_EM(i)*tpnx &
                    &     +exH2y_EM(i)*tpny &
                    &     +exH2z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(exmiu_EM*ztp1*DCONJG(exH2x_EM(i)))
                    tp2 = 0.5d0*REAL(exmiu_EM*ztp1*DCONJG(exH2y_EM(i)))
                    tp3 = 0.5d0*REAL(exmiu_EM*ztp1*DCONJG(exH2z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3

                    ztp1 = exH2x_EM(i)*tpnx &
                    &     +exH2y_EM(i)*tpny &
                    &     +exH2z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(ztp1*DCONJG(exmiu_EM*exH2x_EM(i)))
                    tp2 = 0.5d0*REAL(ztp1*DCONJG(exmiu_EM*exH2y_EM(i)))
                    tp3 = 0.5d0*REAL(ztp1*DCONJG(exmiu_EM*exH2z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3

                    tpndFrcx_EM(i) = tpndFrcx_EM(i) + 0.5d0*vcm_mu0*tp4
                    tpndFrcy_EM(i) = tpndFrcy_EM(i) + 0.5d0*vcm_mu0*tp5
                    tpndFrcz_EM(i) = tpndFrcz_EM(i) + 0.5d0*vcm_mu0*tp6

                    IF (NrmlInOut(ithprtl) == 1) THEN
                        tpndFrcx_EM(i) =-tpndFrcx_EM(i)
                        tpndFrcy_EM(i) =-tpndFrcy_EM(i)
                        tpndFrcz_EM(i) =-tpndFrcz_EM(i)
                    END IF

                    tprx = xnd(i) - xmssctr(ithprtl)
                    tpry = ynd(i) - ymssctr(ithprtl)
                    tprz = znd(i) - zmssctr(ithprtl)
                    tpndTrqx_EM(i) = tpry * tpndFrcz_EM(i) - tprz * tpndFrcy_EM(i)
                    tpndTrqy_EM(i) = tprz * tpndFrcx_EM(i) - tprx * tpndFrcz_EM(i)
                    tpndTrqz_EM(i) = tprx * tpndFrcy_EM(i) - tpry * tpndFrcx_EM(i)

                END DO

                IF (MeshType == 'L') THEN

                    DO k = elstaID(ithprtl), elendID(ithprtl)

                        DO GLQi = 1, n_glqtr2d

                            icnt = n_glqtr2d*(k-1)+GLQi

                            FrcScax_EM(ithprtl) = FrcScax_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcx_EM(elmntlnknd(k,3))

                            FrcScay_EM(ithprtl) = FrcScay_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcy_EM(elmntlnknd(k,3))

                            FrcScaz_EM(ithprtl) = FrcScaz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcz_EM(elmntlnknd(k,3))

                            TrqScax_EM(ithprtl) = TrqScax_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqx_EM(elmntlnknd(k,3))

                            TrqScay_EM(ithprtl) = TrqScay_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqy_EM(elmntlnknd(k,3))

                            TrqScaz_EM(ithprtl) = TrqScaz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqz_EM(elmntlnknd(k,3))

                        END DO

                    END DO

                END IF

                IF (MeshType == 'Q') THEN

                    DO k = elstaID(ithprtl), elendID(ithprtl)

                        DO GLQi = 1, n_glqtr2d

                            icnt = n_glqtr2d*(k-1)+GLQi

                            FrcScax_EM(ithprtl) = FrcScax_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcx_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndFrcx_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndFrcx_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndFrcx_EM(elmntlnknd(k,6))

                            FrcScay_EM(ithprtl) = FrcScay_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcy_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndFrcy_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndFrcy_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndFrcy_EM(elmntlnknd(k,6))

                            FrcScaz_EM(ithprtl) = FrcScaz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcz_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndFrcz_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndFrcz_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndFrcz_EM(elmntlnknd(k,6))

                            TrqScax_EM(ithprtl) = TrqScax_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqx_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndTrqx_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndTrqx_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndTrqx_EM(elmntlnknd(k,6))

                            TrqScay_EM(ithprtl) = TrqScay_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqy_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndTrqy_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndTrqy_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndTrqy_EM(elmntlnknd(k,6))

                            TrqScaz_EM(ithprtl) = TrqScaz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqz_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndTrqz_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndTrqz_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndTrqz_EM(elmntlnknd(k,6))


                        END DO

                    END DO


                END IF

            END IF

        END DO



        ! Incident--scattered interference (extinction) contribution.
        DO ithprtl = 1, nmbrprtl

            FrcExtx_EM(ithprtl) = 0.0d0
            FrcExty_EM(ithprtl) = 0.0d0
            FrcExtz_EM(ithprtl) = 0.0d0

            TrqExtx_EM(ithprtl) = 0.0d0
            TrqExty_EM(ithprtl) = 0.0d0
            TrqExtz_EM(ithprtl) = 0.0d0

            IF (corelnkshell(ithprtl) == 0) THEN

                DO i = ndstaID(ithprtl), ndendID(ithprtl)

                    tpnx = nnx(i)
                    tpny = nny(i)
                    tpnz = nnz(i)

                    tp4 = 0.0d0
                    tp5 = 0.0d0
                    tp6 = 0.0d0

                    tp = 0.5d0*REAL( exeps_EM*exE1x_EM(i) * DCONJG(exE2x_EM(i)) &
                    &               +exeps_EM*exE1y_EM(i) * DCONJG(exE2y_EM(i)) &
                    &               +exeps_EM*exE1z_EM(i) * DCONJG(exE2z_EM(i))  ) &
                    &   +0.5d0*REAL( exeps_EM*exE2x_EM(i) * DCONJG(exE1x_EM(i)) &
                    &               +exeps_EM*exE2y_EM(i) * DCONJG(exE1y_EM(i)) &
                    &               +exeps_EM*exE2z_EM(i) * DCONJG(exE1z_EM(i))  )
                    tp4 = tp4 - tp*tpnx
                    tp5 = tp5 - tp*tpny
                    tp6 = tp6 - tp*tpnz

                    ztp1 = exE1x_EM(i)*tpnx &
                    &     +exE1y_EM(i)*tpny &
                    &     +exE1z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(exeps_EM*ztp1*DCONJG(exE2x_EM(i)))
                    tp2 = 0.5d0*REAL(exeps_EM*ztp1*DCONJG(exE2y_EM(i)))
                    tp3 = 0.5d0*REAL(exeps_EM*ztp1*DCONJG(exE2z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3
                    ztp1 = exE2x_EM(i)*tpnx &
                    &     +exE2y_EM(i)*tpny &
                    &     +exE2z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(exeps_EM*ztp1*DCONJG(exE1x_EM(i)))
                    tp2 = 0.5d0*REAL(exeps_EM*ztp1*DCONJG(exE1y_EM(i)))
                    tp3 = 0.5d0*REAL(exeps_EM*ztp1*DCONJG(exE1z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3

                    ztp1 = exE1x_EM(i)*tpnx &
                    &     +exE1y_EM(i)*tpny &
                    &     +exE1z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(ztp1*DCONJG(exeps_EM*exE2x_EM(i)))
                    tp2 = 0.5d0*REAL(ztp1*DCONJG(exeps_EM*exE2y_EM(i)))
                    tp3 = 0.5d0*REAL(ztp1*DCONJG(exeps_EM*exE2z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3
                    ztp1 = exE2x_EM(i)*tpnx &
                    &     +exE2y_EM(i)*tpny &
                    &     +exE2z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(ztp1*DCONJG(exeps_EM*exE1x_EM(i)))
                    tp2 = 0.5d0*REAL(ztp1*DCONJG(exeps_EM*exE1y_EM(i)))
                    tp3 = 0.5d0*REAL(ztp1*DCONJG(exeps_EM*exE1z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3

                    tpndFrcx_EM(i) = 0.5d0*vcm_eps0*tp4
                    tpndFrcy_EM(i) = 0.5d0*vcm_eps0*tp5
                    tpndFrcz_EM(i) = 0.5d0*vcm_eps0*tp6

                    tp4 = 0.0d0
                    tp5 = 0.0d0
                    tp6 = 0.0d0

                    tp = 0.5d0*REAL( exmiu_EM*exH1x_EM(i) * DCONJG(exH2x_EM(i)) &
                    &               +exmiu_EM*exH1y_EM(i) * DCONJG(exH2y_EM(i)) &
                    &               +exmiu_EM*exH1z_EM(i) * DCONJG(exH2z_EM(i))  ) &
                    &   +0.5d0*REAL( exmiu_EM*exH2x_EM(i) * DCONJG(exH1x_EM(i)) &
                    &               +exmiu_EM*exH2y_EM(i) * DCONJG(exH1y_EM(i)) &
                    &               +exmiu_EM*exH2z_EM(i) * DCONJG(exH1z_EM(i))  )
                    tp4 = tp4 - tp*tpnx
                    tp5 = tp5 - tp*tpny
                    tp6 = tp6 - tp*tpnz

                    ztp1 = exH1x_EM(i)*tpnx &
                    &     +exH1y_EM(i)*tpny &
                    &     +exH1z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(exmiu_EM*ztp1*DCONJG(exH2x_EM(i)))
                    tp2 = 0.5d0*REAL(exmiu_EM*ztp1*DCONJG(exH2y_EM(i)))
                    tp3 = 0.5d0*REAL(exmiu_EM*ztp1*DCONJG(exH2z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3
                    ztp1 = exH2x_EM(i)*tpnx &
                    &     +exH2y_EM(i)*tpny &
                    &     +exH2z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(exmiu_EM*ztp1*DCONJG(exH1x_EM(i)))
                    tp2 = 0.5d0*REAL(exmiu_EM*ztp1*DCONJG(exH1y_EM(i)))
                    tp3 = 0.5d0*REAL(exmiu_EM*ztp1*DCONJG(exH1z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3

                    ztp1 = exH1x_EM(i)*tpnx &
                    &     +exH1y_EM(i)*tpny &
                    &     +exH1z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(ztp1*DCONJG(exmiu_EM*exH2x_EM(i)))
                    tp2 = 0.5d0*REAL(ztp1*DCONJG(exmiu_EM*exH2y_EM(i)))
                    tp3 = 0.5d0*REAL(ztp1*DCONJG(exmiu_EM*exH2z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3
                    ztp1 = exH2x_EM(i)*tpnx &
                    &     +exH2y_EM(i)*tpny &
                    &     +exH2z_EM(i)*tpnz
                    tp1 = 0.5d0*REAL(ztp1*DCONJG(exmiu_EM*exH1x_EM(i)))
                    tp2 = 0.5d0*REAL(ztp1*DCONJG(exmiu_EM*exH1y_EM(i)))
                    tp3 = 0.5d0*REAL(ztp1*DCONJG(exmiu_EM*exH1z_EM(i)))
                    tp4 = tp4 + tp1
                    tp5 = tp5 + tp2
                    tp6 = tp6 + tp3

                    tpndFrcx_EM(i) = tpndFrcx_EM(i) + 0.5d0*vcm_mu0*tp4
                    tpndFrcy_EM(i) = tpndFrcy_EM(i) + 0.5d0*vcm_mu0*tp5
                    tpndFrcz_EM(i) = tpndFrcz_EM(i) + 0.5d0*vcm_mu0*tp6

                    IF (NrmlInOut(ithprtl) == 1) THEN
                        tpndFrcx_EM(i) =-tpndFrcx_EM(i)
                        tpndFrcy_EM(i) =-tpndFrcy_EM(i)
                        tpndFrcz_EM(i) =-tpndFrcz_EM(i)
                    END IF

                    tprx = xnd(i) - xmssctr(ithprtl)
                    tpry = ynd(i) - ymssctr(ithprtl)
                    tprz = znd(i) - zmssctr(ithprtl)
                    tpndTrqx_EM(i) = tpry * tpndFrcz_EM(i) - tprz * tpndFrcy_EM(i)
                    tpndTrqy_EM(i) = tprz * tpndFrcx_EM(i) - tprx * tpndFrcz_EM(i)
                    tpndTrqz_EM(i) = tprx * tpndFrcy_EM(i) - tpry * tpndFrcx_EM(i)

                END DO

                IF (MeshType == 'L') THEN

                    DO k = elstaID(ithprtl), elendID(ithprtl)

                        DO GLQi = 1, n_glqtr2d

                            icnt = n_glqtr2d*(k-1)+GLQi

                            FrcExtx_EM(ithprtl) = FrcExtx_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcx_EM(elmntlnknd(k,3))

                            FrcExty_EM(ithprtl) = FrcExty_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcy_EM(elmntlnknd(k,3))

                            FrcExtz_EM(ithprtl) = FrcExtz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcz_EM(elmntlnknd(k,3))

                            TrqExtx_EM(ithprtl) = TrqExtx_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqx_EM(elmntlnknd(k,3))

                            TrqExty_EM(ithprtl) = TrqExty_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqy_EM(elmntlnknd(k,3))

                            TrqExtz_EM(ithprtl) = TrqExtz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqz_EM(elmntlnknd(k,3))

                        END DO

                    END DO

                END IF

                IF (MeshType == 'Q') THEN

                    DO k = elstaID(ithprtl), elendID(ithprtl)

                        DO GLQi = 1, n_glqtr2d

                            icnt = n_glqtr2d*(k-1)+GLQi

                            FrcExtx_EM(ithprtl) = FrcExtx_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcx_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndFrcx_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndFrcx_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndFrcx_EM(elmntlnknd(k,6))

                            FrcExty_EM(ithprtl) = FrcExty_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcy_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndFrcy_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndFrcy_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndFrcy_EM(elmntlnknd(k,6))

                            FrcExtz_EM(ithprtl) = FrcExtz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndFrcz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndFrcz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndFrcz_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndFrcz_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndFrcz_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndFrcz_EM(elmntlnknd(k,6))

                            TrqExtx_EM(ithprtl) = TrqExtx_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqx_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqx_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqx_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndTrqx_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndTrqx_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndTrqx_EM(elmntlnknd(k,6))

                            TrqExty_EM(ithprtl) = TrqExty_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqy_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqy_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqy_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndTrqy_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndTrqy_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndTrqy_EM(elmntlnknd(k,6))

                            TrqExtz_EM(ithprtl) = TrqExtz_EM(ithprtl) &
                            &   + srcfmm_wtnd(1,icnt) * tpndTrqz_EM(elmntlnknd(k,1)) &
                            &   + srcfmm_wtnd(2,icnt) * tpndTrqz_EM(elmntlnknd(k,2)) &
                            &   + srcfmm_wtnd(3,icnt) * tpndTrqz_EM(elmntlnknd(k,3)) &
                            &   + srcfmm_wtnd(4,icnt) * tpndTrqz_EM(elmntlnknd(k,4)) &
                            &   + srcfmm_wtnd(5,icnt) * tpndTrqz_EM(elmntlnknd(k,5)) &
                            &   + srcfmm_wtnd(6,icnt) * tpndTrqz_EM(elmntlnknd(k,6))


                        END DO

                    END DO


                END IF

            END IF

        END DO

        ! Recombine the incident, scattered and interference decomposition.
        DO ithprtl = 1, nmbrprtl

            FrcTotx_EM(ithprtl) = FrcIncx_EM(ithprtl) + FrcScax_EM(ithprtl) + FrcExtx_EM(ithprtl)
            FrcToty_EM(ithprtl) = FrcIncy_EM(ithprtl) + FrcScay_EM(ithprtl) + FrcExty_EM(ithprtl)
            FrcTotz_EM(ithprtl) = FrcIncz_EM(ithprtl) + FrcScaz_EM(ithprtl) + FrcExtz_EM(ithprtl)

            TrqTotx_EM(ithprtl) = TrqIncx_EM(ithprtl) + TrqScax_EM(ithprtl) + TrqExtx_EM(ithprtl)
            TrqToty_EM(ithprtl) = TrqIncy_EM(ithprtl) + TrqScay_EM(ithprtl) + TrqExty_EM(ithprtl)
            TrqTotz_EM(ithprtl) = TrqIncz_EM(ithprtl) + TrqScaz_EM(ithprtl) + TrqExtz_EM(ithprtl)

        END DO

        ! Recombine the electric and magnetic decomposition.
        DO ithprtl = 1, nmbrprtl

            FrcTtlx_EM(ithprtl) = FrcElcx_EM(ithprtl) + FrcMagx_EM(ithprtl)
            FrcTtly_EM(ithprtl) = FrcElcy_EM(ithprtl) + FrcMagy_EM(ithprtl)
            FrcTtlz_EM(ithprtl) = FrcElcz_EM(ithprtl) + FrcMagz_EM(ithprtl)

            TrqTtlx_EM(ithprtl) = TrqElcx_EM(ithprtl) + TrqMagx_EM(ithprtl)
            TrqTtly_EM(ithprtl) = TrqElcy_EM(ithprtl) + TrqMagy_EM(ithprtl)
            TrqTtlz_EM(ithprtl) = TrqElcz_EM(ithprtl) + TrqMagz_EM(ithprtl)

        END DO

    END SUBROUTINE

END MODULE
