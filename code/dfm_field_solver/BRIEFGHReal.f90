
! SPDX-FileCopyrightText: 2026 Qiang Sun
! SPDX-License-Identifier: BSD-3-Clause

MODULE BRIEFGHReal

    ! Element-level quadrature kernels for the real modified-Helmholtz BRIEF
    ! formulation.  For an integration point x, collocation point x0, and
    ! R = |x-x0|, this module uses Gk = exp(-k*R)/R; the conventional factor
    ! 1/(4*pi) is omitted.  Hk denotes dGk/dn(x), with n(x) taken from
    ! srcfmm_nrm.  The sign is fixed by tpdp = n(x).(x0-x), so Hk is not a
    ! derivative along the collocation normal n(x0).
    !
    ! The G and H arrays contain integrals of each local shape function times
    ! Gk and Hk, respectively.  The two scalar NSBIM outputs are the element
    ! contributions to the k-independent Laplace subtraction used to cancel
    ! the coincident-point singularities in the assembled boundary equation:
    !   gnsbim = integral[(n0.(x-x0))*dG0/dn(x) - (n0.n(x))*G0] dS,
    !   hnsbim = integral[-dG0/dn(x)] dS,             G0 = 1/R.

    USE Pre_Constants

    USE Geom_GlobalData

    IMPLICIT NONE

    CONTAINS

    ! Integrate Gk and Hk against the three linear-triangle shape functions.
    ! Outputs a:b:c follow local element nodes 1:2:3; the scalar outputs are
    ! the regularising Laplace integrals defined in the module header.
    SUBROUTINE CalGHLnrBRIEFLnrREAL(dmKwn, dmelmnid, &
    &                               dmp0x, dmp0y, dmp0z, dmp0nnx, dmp0nny, dmp0nnz, &
    &                               dmga, dmgb, dmgc, &
    &                               dmha, dmhb, dmhc, &
    &                               dmgnsbim, dmhnsbim )

        DOUBLE PRECISION, INTENT(IN) ::  dmKwn
        INTEGER, INTENT(IN) ::  dmelmnid
        DOUBLE PRECISION, INTENT(IN) ::  dmp0x, dmp0y, dmp0z, dmp0nnx, dmp0nny, dmp0nnz

        DOUBLE PRECISION, INTENT(OUT) :: dmga, dmgb, dmgc

        DOUBLE PRECISION, INTENT(OUT) :: dmha, dmhb, dmhc

        DOUBLE PRECISION, INTENT(OUT) :: dmgnsbim, dmhnsbim

        DOUBLE PRECISION :: tpga, tpgb, tpgc, tpha, tphb, tphc, tpgnsbim, tphnsbim
        DOUBLE PRECISION :: tpGrnK, tpdGKdn, tpGrn0, tpdG0dn
        DOUBLE PRECISION :: tpkwni, tpExpkwnr
        DOUBLE PRECISION :: tpn0xx0, tpn0En

        DOUBLE PRECISION :: tpEnx, tpEny, tpEnz

        DOUBLE PRECISION :: tprx, tpry, tprz
        DOUBLE PRECISION :: tprr0x, tprr0y, tprr0z, mdl_r0r, &
        &                   over_r0r1, over_r0r2, over_r0r3, tpdp

        INTEGER :: GLQi, icnt

        tpkwni = (-1.0d0)*dmKwn

        tpga = 0.0d0
        tpgb = 0.0d0
        tpgc = 0.0d0

        tpha = 0.0d0
        tphb = 0.0d0
        tphc = 0.0d0

        tpgnsbim = 0.0d0
        tphnsbim = 0.0d0

        DO GLQi = 1, n_glqtr2d

            icnt = n_glqtr2d*(dmelmnid-1)+GLQi

            tprx = srcfmm_vec(1, icnt)
            tpry = srcfmm_vec(2, icnt)
            tprz = srcfmm_vec(3, icnt)

            tpEnx = srcfmm_nrm(1, icnt)
            tpEny = srcfmm_nrm(2, icnt)
            tpEnz = srcfmm_nrm(3, icnt)

            tprr0x = tprx - dmp0x
            tprr0y = tpry - dmp0y
            tprr0z = tprz - dmp0z

            mdl_r0r = DSQRT(tprr0x**2 + tprr0y**2 + tprr0z**2)
            over_r0r1 = 1.0d0/mdl_r0r
            over_r0r2 = over_r0r1*over_r0r1
            over_r0r3 = over_r0r2*over_r0r1
            tpdp = -(tpEnx*tprr0x + tpEny*tprr0y + tpEnz*tprr0z)

            tpExpkwnr = DEXP(tpkwni*mdl_r0r)

            tpGrnK = over_r0r1*tpExpkwnr
            tpdGKdn = tpdp*(over_r0r1-tpkwni)*tpExpkwnr*over_r0r2
            tpGrn0 = over_r0r1
            tpdG0dn = tpdp*over_r0r3

            tpga = tpga + srcfmm_wtnd(1,icnt) * tpGrnK
            tpgb = tpgb + srcfmm_wtnd(2,icnt) * tpGrnK
            tpgc = tpgc + srcfmm_wtnd(3,icnt) * tpGrnK

            tpn0xx0 = dmp0nnx*tprr0x+dmp0nny*tprr0y+dmp0nnz*tprr0z
            tpn0En  = dmp0nnx*tpEnx+dmp0nny*tpEny+dmp0nnz*tpEnz

            tpgnsbim = tpgnsbim + srcfmm_wght(icnt) * ( tpn0xx0*tpdG0dn -tpn0En*tpGrn0 )

            tpha = tpha + srcfmm_wtnd(1,icnt) * tpdGKdn
            tphb = tphb + srcfmm_wtnd(2,icnt) * tpdGKdn
            tphc = tphc + srcfmm_wtnd(3,icnt) * tpdGKdn

            tphnsbim = tphnsbim + srcfmm_wght(icnt) * (-tpdG0dn)


        END DO

        dmga = tpga
        dmgb = tpgb
        dmgc = tpgc

        dmha = tpha
        dmhb = tphb
        dmhc = tphc

        dmgnsbim = tpgnsbim
        dmhnsbim = tphnsbim

    END SUBROUTINE

    ! Integrate Gk and Hk against the six quadratic-triangle (Q6) shape
    ! functions.  Outputs a:b:c:d:e:f follow local nodes 1:2:3:4:5:6;
    ! the scalar outputs are the same shape-function-independent BRIEF terms.
    SUBROUTINE CalGHQdrBRIEFLnrREAL(dmKwn, dmelmnid, &
    &                               dmp0x, dmp0y, dmp0z, dmp0nnx, dmp0nny, dmp0nnz,  &
    &                               dmga, dmgb, dmgc, dmgd, dmge, dmgf, &
    &                               dmha, dmhb, dmhc, dmhd, dmhe, dmhf, &
    &                               dmgnsbim, dmhnsbim )

        DOUBLE PRECISION, INTENT(IN) ::  dmKwn
        INTEGER, INTENT(IN) ::  dmelmnid
        DOUBLE PRECISION, INTENT(IN) ::  dmp0x, dmp0y, dmp0z, dmp0nnx, dmp0nny, dmp0nnz

        DOUBLE PRECISION, INTENT(OUT) :: dmga, dmgb, dmgc, dmgd, dmge, dmgf

        DOUBLE PRECISION, INTENT(OUT) :: dmha, dmhb, dmhc, dmhd, dmhe, dmhf

        DOUBLE PRECISION, INTENT(OUT) :: dmgnsbim, dmhnsbim

        DOUBLE PRECISION :: tpga, tpgb, tpgc, tpgd, tpge, tpgf
        DOUBLE PRECISION :: tpha, tphb, tphc, tphd, tphe, tphf
        DOUBLE PRECISION :: tpgnsbim, tphnsbim
        DOUBLE PRECISION :: tpGrnK, tpdGKdn, tpGrn0, tpdG0dn
        DOUBLE PRECISION :: tpkwni, tpExpkwnr
        DOUBLE PRECISION :: tpn0xx0, tpn0En

        DOUBLE PRECISION :: tprx, tpry, tprz
        DOUBLE PRECISION :: tprr0x, tprr0y, tprr0z, mdl_r0r, &
        &                   over_r0r1, over_r0r2, over_r0r3, tpdp

        DOUBLE PRECISION :: drx_deps,drx_dyet,dry_deps,dry_dyet,drz_deps,drz_dyet
        DOUBLE PRECISION :: tpEnx, tpEny, tpEnz, JcbDtmn
        DOUBLE PRECISION :: tpxieta

        INTEGER :: GLQi, icnt

        tpkwni = (-1.0d0)*dmKwn

        tpga = 0.0d0
        tpgb = 0.0d0
        tpgc = 0.0d0
        tpgd = 0.0d0
        tpge = 0.0d0
        tpgf = 0.0d0

        tpha = 0.0d0
        tphb = 0.0d0
        tphc = 0.0d0
        tphd = 0.0d0
        tphe = 0.0d0
        tphf = 0.0d0

        tpgnsbim = 0.0d0
        tphnsbim = 0.0d0

        DO GLQi = 1, n_glqtr2d

            icnt = n_glqtr2d*(dmelmnid-1)+GLQi

            tprx = srcfmm_vec(1, icnt)
            tpry = srcfmm_vec(2, icnt)
            tprz = srcfmm_vec(3, icnt)

            tpEnx = srcfmm_nrm(1, icnt)
            tpEny = srcfmm_nrm(2, icnt)
            tpEnz = srcfmm_nrm(3, icnt)

            tprr0x = tprx - dmp0x
            tprr0y = tpry - dmp0y
            tprr0z = tprz - dmp0z

            mdl_r0r = DSQRT(tprr0x**2 + tprr0y**2 + tprr0z**2)
            over_r0r1 = 1.0d0/mdl_r0r
            over_r0r2 = over_r0r1*over_r0r1
            over_r0r3 = over_r0r1*over_r0r2
            tpdp = -(tpEnx*tprr0x + tpEny*tprr0y + tpEnz*tprr0z)

            tpExpkwnr = DEXP(tpkwni*mdl_r0r)

            tpGrnK = over_r0r1*tpExpkwnr
            tpdGKdn = tpdp*(over_r0r1-tpkwni)*tpExpkwnr*over_r0r2
            tpGrn0 = over_r0r1
            tpdG0dn = tpdp*over_r0r3

            tpn0xx0 = dmp0nnx*tprr0x+dmp0nny*tprr0y+dmp0nnz*tprr0z
            tpn0En  = dmp0nnx*tpEnx+dmp0nny*tpEny+dmp0nnz*tpEnz


            tpga = tpga + srcfmm_wtnd(1,icnt) * tpGrnK
            tpgb = tpgb + srcfmm_wtnd(2,icnt) * tpGrnK
            tpgc = tpgc + srcfmm_wtnd(3,icnt) * tpGrnK
            tpgd = tpgd + srcfmm_wtnd(4,icnt) * tpGrnK
            tpge = tpge + srcfmm_wtnd(5,icnt) * tpGrnK
            tpgf = tpgf + srcfmm_wtnd(6,icnt) * tpGrnK

            tpha = tpha + srcfmm_wtnd(1,icnt) * tpdGKdn
            tphb = tphb + srcfmm_wtnd(2,icnt) * tpdGKdn
            tphc = tphc + srcfmm_wtnd(3,icnt) * tpdGKdn
            tphd = tphd + srcfmm_wtnd(4,icnt) * tpdGKdn
            tphe = tphe + srcfmm_wtnd(5,icnt) * tpdGKdn
            tphf = tphf + srcfmm_wtnd(6,icnt) * tpdGKdn

            tpgnsbim = tpgnsbim + srcfmm_wght(icnt) * ( tpn0xx0*tpdG0dn -tpn0En*tpGrn0 )

            tphnsbim = tphnsbim + srcfmm_wght(icnt) * (-tpdG0dn)


        END DO

        dmga = tpga
        dmgb = tpgb
        dmgc = tpgc
        dmgd = tpgd
        dmge = tpge
        dmgf = tpgf

        dmha = tpha
        dmhb = tphb
        dmhc = tphc
        dmhd = tphd
        dmhe = tphe
        dmhf = tphf

        dmgnsbim = tpgnsbim
        dmhnsbim = tphnsbim

    END SUBROUTINE
END MODULE
