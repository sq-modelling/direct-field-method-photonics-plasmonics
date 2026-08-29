
! SPDX-FileCopyrightText: 2026 Qiang Sun
! SPDX-License-Identifier: BSD-3-Clause

MODULE Pre_Constants

    IMPLICIT NONE


    DOUBLE PRECISION, PARAMETER :: pai = DATAN(1.0d0)*4.0d0

    COMPLEX(KIND=KIND(1.0D0)), PARAMETER :: ztpzero = DCMPLX(0.0d0, 0.0d0), &
    &                                       ztponei = DCMPLX(0.0d0, 1.0d0), &
    &                                       ztponej = DCMPLX(0.0d0,-1.0d0), &
    &                                       ztpone  = DCMPLX(1.0d0, 0.0d0)


    ! The field solver uses normalised vacuum constants: eps0 = mu0 = c = 1.
    ! A separate SI epsilon0 value is retained for downstream compatibility.
    DOUBLE PRECISION, PARAMETER :: vcm_eps0 = 1.0d0, vcm_mu0=1.0d0
    DOUBLE PRECISION, PARAMETER :: vcm_lgtspd = 1.0d0
    DOUBLE PRECISION, PARAMETER :: vcm_sqrteps0mu0 = DSQRT(vcm_eps0/vcm_mu0)
    DOUBLE PRECISION, PARAMETER :: vcm_sqrtmu0eps0 = DSQRT(vcm_mu0/vcm_eps0)

    DOUBLE PRECISION, PARAMETER :: epsilon0 = 8.854187817620E-12    !in C^2/(m^2*N)

    DOUBLE PRECISION, PARAMETER :: NAvog = 6.02214085774E23

    DOUBLE PRECISION, PARAMETER :: BoltzmannC = 1.3805E-23          !in J/K

    DOUBLE PRECISION, PARAMETER :: elcchrg = 1.602176565E-19        !in C
    DOUBLE PRECISION, PARAMETER :: elcmass = 9.10938356E-31    !in kg

    DOUBLE PRECISION, PARAMETER :: Planckh = 6.62607004E-34    !in J*s
    DOUBLE PRECISION, PARAMETER :: Planckhbar = Planckh/(2.0d0*pai)

END MODULE
