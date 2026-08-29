
! SPDX-FileCopyrightText: 2026 Qiang Sun
! SPDX-License-Identifier: BSD-3-Clause

MODULE EM_SurfCal_Input

! Read the fixed-layout electromagnetic input file and initialise shared material,
! excitation, boundary, source, and post-processing data.
!
! Input_Phys_EM.dat is record-oriented: no-list READ statements deliberately consume
! description, blank, and marker records.  Keep those records in the documented order.
! All coordinates, extents, radii, and wavelengths use the mesh length unit; a supplied
! wavenumber uses its inverse.  A supplied angular frequency must use the consistent
! unit implied by vcm_eps0 and vcm_mu0.  The code performs no unit conversion.  Material
! n,k pairs and direction vectors are dimensionless.  phase_EM is supplied in degrees.

    USE omp_lib

    USE Pre_Constants

    USE Geom_GlobalData

    USE EM_SurfCal_GlobalData

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE GetPhysInputInt_EM

! Purpose: read the 18-record global header of Input_Phys_EM.dat and derive the
!          normalised propagation direction, free-space wavelength/wavenumber,
!          angular frequency, complex exterior material, and phase increment.
! Output : module data in EM_SurfCal_GlobalData and slice controls in Geom_GlobalData.
! Failure: stop if Input_Phys_EM.dat cannot be opened.

        INTEGER :: IOS, ithprtl

        DOUBLE PRECISION :: tp,tp1,tp2,tp3,tp4,tp5,tp6,tp7,tp8,tp9,tpread
        COMPLEX(KIND=KIND(1.0D0)) :: ztp,ztp1,ztp2,ztp3,ztp4,ztp5,ztp6,ztp7,ztp8,ztp9

        OPEN (71, FILE = "Input_Phys_EM.dat", STATUS = "OLD", IOSTAT = IOS)
        IF (IOS /= 0) THEN
            PRINT*, "'Input_Phys_EM.dat' does not exist! Please check!"
            STOP
        END IF

! Records 1-2: description, then slice flag, plane selector, centre, full extents,
! and grid count.  The selector is xy, yz, zx, or 3D; PostEgnd_RMpln is the count
! along the first axis of a selected plane and the second count preserves aspect ratio.
        READ (71, *)
        READ (71, *) on_fldcal,Postxyz_RMpln,PostOfstx_RMpln,PostOfsty_RMpln,PostOfstz_RMpln,&
                    &PostEgSZx_RMpln,PostEgSZy_RMpln,PostEgSZz_RMpln,PostEgnd_RMpln
        READ (71, *)
! Records 4-5: description, then enhancement flag, cross-section flag, sampling-sphere
! radius, and centre.  The radius and centre use the mesh length unit.
        READ (71, *)
        READ (71, *) on_enhance, on_xs,size_xs,xcen_xs,ycen_xs,zcen_xs
        READ (71, *)
! Records 7-8: description, then optional pulse/FFT controls in the following READ order.
        READ (71, *)
        READ (71, *) FFTpulse_EM,FFTpulsetype_EM,FFTnmbrFrq_EM,FFTWidth_EM,&
        &            FFTnmbrBIM_EM,FFTnmbrOsPulse_EM,FFTalfa_EM,FFTphsPulse_EM,FFTBIMCal_EM
        READ (71, *)
! Records 10-11: description, then excitation code, polarisation code, component
! receiving the relative phase, relative phase/feature value, spectral selector and
! value, incident-field vector, propagation vector, and output phase step in degrees.
! wlorwn_EM = l, k, or w selects wavelength, wavenumber, or angular frequency.
! incFeature_EM is in radians when it is used as a component-to-component phase.
! For pwe/pwh, incOrder_EM = 1, 2, or 3 applies that phase to y, z, or x, respectively.
        READ (71, *)
        READ (71, *) excitetype_EM, poltype_EM, incOrder_EM, incFeature_EM, &
        &            wlorwn_EM, tpread, &
        &            incFieldx_EM, incFieldy_EM, incFieldz_EM, inckx_EM, incky_EM, inckz_EM, &
        &            phase_EM
        READ (71, *)
! Records 13-14: description, then electric n,k, magnetic n,k, and prism n,k for
! the parent/exterior medium.  Relative eps and mu are formed as (n + i*kappa)**2.
        READ (71, *)
        READ (71, *) exepsn_EM, exepsk_EM, exmiun_EM, exmiuk_EM, &
        &            prismn_EM, prismk_EM
        READ (71, *)
! Records 16-17: description, then the number of sources in domain 0.
        READ (71, *)
        READ (71, *) exnmbrsrc_EM
        READ (71, *)

        CLOSE (71)

        IF (excitetype_EM == 'spe' .OR. excitetype_EM == 'sph') THEN
            AngFrqnc_EM = 0.0d0
            vcmwn_EM = 0.0d0
            vcmwl_EM = 1E18
        END IF

        IF (     excitetype_EM == 'pwe' .OR. excitetype_EM == 'pwh' &
        &   .OR. excitetype_EM == 'swe' .OR. excitetype_EM == 'swh' &
        &   .OR. excitetype_EM == 'eva' .OR. excitetype_EM == 'gau' &
        &   .OR. excitetype_EM == 'gb5' .OR. excitetype_EM == 'bsl' &
        &   .OR. excitetype_EM == 'mch' ) THEN
            IF (wlorwn_EM == "l") THEN
                vcmwl_EM = tpread
                vcmwn_EM = 2.0d0*pai/vcmwl_EM
                AngFrqnc_EM = vcmwn_EM/DSQRT(vcm_eps0*vcm_mu0)
            END IF
            IF (wlorwn_EM == "k") THEN
                vcmwn_EM = tpread
                vcmwl_EM = 2.0d0*pai/vcmwn_EM
                AngFrqnc_EM = vcmwn_EM/DSQRT(vcm_eps0*vcm_mu0)
            END IF
            IF (wlorwn_EM == "w") THEN
                AngFrqnc_EM = tpread
                vcmwn_EM = AngFrqnc_EM*DSQRT(vcm_eps0*vcm_mu0)
                vcmwl_EM = 2.0d0*pai/vcmwn_EM
            END IF
        END IF

        incFieldmdl_EM = DSQRT(incFieldx_EM**2 + incFieldy_EM**2 + incFieldz_EM**2)

        IF (excitetype_EM == 'pwe' .OR. excitetype_EM == 'pwh') THEN
            IF (incOrder_EM /= 1 .AND. incOrder_EM /= 2 .AND. incOrder_EM /= 3) incOrder_EM = 1
            IF (incOrder_EM == 1) THEN
                tp1 =  DABS(incFieldx_EM)
                tp2 = CDABS(incFieldy_EM*CDEXP(ztponei*incFeature_EM))
                tp3 =  DABS(incFieldz_EM)
                tp = DSQRT(tp1**2+tp2**2+tp3**2)
                incFieldmdl_EM = tp
                tp = 1.0d0/tp
                incFieldx_EM = incFieldx_EM*tp
                incFieldy_EM = incFieldy_EM*tp
                incFieldz_EM = incFieldz_EM*tp
            END IF
            IF (incOrder_EM == 2) THEN
                tp1 =  DABS(incFieldx_EM)
                tp2 =  DABS(incFieldy_EM)
                tp3 = CDABS(incFieldz_EM*CDEXP(ztponei*incFeature_EM))
                tp = DSQRT(tp1**2+tp2**2+tp3**2)
                incFieldmdl_EM = tp
                tp = 1.0d0/tp
                incFieldx_EM = incFieldx_EM*tp
                incFieldy_EM = incFieldy_EM*tp
                incFieldz_EM = incFieldz_EM*tp
            END IF
            IF (incOrder_EM == 3) THEN
                tp1 = CDABS(incFieldx_EM*CDEXP(ztponei*incFeature_EM))
                tp2 =  DABS(incFieldy_EM)
                tp3 =  DABS(incFieldz_EM)
                tp = DSQRT(tp1**2+tp2**2+tp3**2)
                incFieldmdl_EM = tp
                tp = 1.0d0/tp
                incFieldx_EM = incFieldx_EM*tp
                incFieldy_EM = incFieldy_EM*tp
                incFieldz_EM = incFieldz_EM*tp
            END IF
        END IF

        ztp = DCMPLX(exmiun_EM, exmiuk_EM)
        exmiu_EM = ztp*ztp

        ztp = DCMPLX(exepsn_EM, exepsk_EM)
        exeps_EM = ztp*ztp

        ztp = exeps_EM*exmiu_EM
        exk_EM = vcmwn_EM*CDSQRT(ztp)

        IF (ABS(phase_EM) < 0.1d0) phase_EM = 5.0d0
        phase_EM = phase_EM*pai/180.0d0

        tp = DSQRT(inckx_EM**2+incky_EM**2+inckz_EM**2)
        inckx_EM = inckx_EM/tp
        incky_EM = incky_EM/tp
        inckz_EM = inckz_EM/tp

    END SUBROUTINE


    SUBROUTINE GetPhysInput_EM

! Purpose: allocate all per-surface EM arrays and read one fixed 10-record block for
!          each of the nmbrprtl surfaces after the 18-record global header.
! Input  : nmbrprtl and the global spectral data must already be initialised.
! Output : boundary records, owned-domain material properties and wavenumbers, source
!          counts/ranges, and storage for later force and enhancement data.
! Contract: each surface block contains start marker; boundary description/data;
!           separator; material description/data; separator; source description/count;
!           and an end marker.  Marker text is consumed but not parsed.

        INTEGER :: i, j, k, IOS, ithprtl

        DOUBLE PRECISION :: tp,tp1,tp2,tp3,tp4,tp5,tp6,tp7,tp8,tp9,tpread
        COMPLEX(KIND=KIND(1.0D0)) :: ztp,ztp1,ztp2,ztp3,ztp4,ztp5,ztp6,ztp7,ztp8,ztp9

        ALLOCATE (BCType_EM(nmbrprtl))
        ALLOCATE (BCValue_EM(nmbrprtl))
        ALLOCATE (BCRead_EM(nmbrprtl))

        ALLOCATE (exsurfQ1_EM(nmbrprtl))
        ALLOCATE (exsurfQ2_EM(nmbrprtl))
        ALLOCATE (exsurfQ3_EM(nmbrprtl))
        ALLOCATE (insurfQ1_EM(nmbrprtl))
        ALLOCATE (insurfQ2_EM(nmbrprtl))
        ALLOCATE (insurfQ3_EM(nmbrprtl))

        ALLOCATE (Frcx_EM(nmbrprtl))
        ALLOCATE (Frcy_EM(nmbrprtl))
        ALLOCATE (Frcz_EM(nmbrprtl))
        ALLOCATE (Trqx_EM(nmbrprtl))
        ALLOCATE (Trqy_EM(nmbrprtl))
        ALLOCATE (Trqz_EM(nmbrprtl))
        ALLOCATE (FrcElcx_EM(nmbrprtl))
        ALLOCATE (FrcElcy_EM(nmbrprtl))
        ALLOCATE (FrcElcz_EM(nmbrprtl))
        ALLOCATE (TrqElcx_EM(nmbrprtl))
        ALLOCATE (TrqElcy_EM(nmbrprtl))
        ALLOCATE (TrqElcz_EM(nmbrprtl))
        ALLOCATE (FrcMagx_EM(nmbrprtl))
        ALLOCATE (FrcMagy_EM(nmbrprtl))
        ALLOCATE (FrcMagz_EM(nmbrprtl))
        ALLOCATE (TrqMagx_EM(nmbrprtl))
        ALLOCATE (TrqMagy_EM(nmbrprtl))
        ALLOCATE (TrqMagz_EM(nmbrprtl))
        ALLOCATE (FrcTtlx_EM(nmbrprtl))
        ALLOCATE (FrcTtly_EM(nmbrprtl))
        ALLOCATE (FrcTtlz_EM(nmbrprtl))
        ALLOCATE (TrqTtlx_EM(nmbrprtl))
        ALLOCATE (TrqTtly_EM(nmbrprtl))
        ALLOCATE (TrqTtlz_EM(nmbrprtl))
        ALLOCATE (FrcIncx_EM(nmbrprtl))
        ALLOCATE (FrcIncy_EM(nmbrprtl))
        ALLOCATE (FrcIncz_EM(nmbrprtl))
        ALLOCATE (TrqIncx_EM(nmbrprtl))
        ALLOCATE (TrqIncy_EM(nmbrprtl))
        ALLOCATE (TrqIncz_EM(nmbrprtl))
        ALLOCATE (FrcScax_EM(nmbrprtl))
        ALLOCATE (FrcScay_EM(nmbrprtl))
        ALLOCATE (FrcScaz_EM(nmbrprtl))
        ALLOCATE (TrqScax_EM(nmbrprtl))
        ALLOCATE (TrqScay_EM(nmbrprtl))
        ALLOCATE (TrqScaz_EM(nmbrprtl))
        ALLOCATE (FrcExtx_EM(nmbrprtl))
        ALLOCATE (FrcExty_EM(nmbrprtl))
        ALLOCATE (FrcExtz_EM(nmbrprtl))
        ALLOCATE (TrqExtx_EM(nmbrprtl))
        ALLOCATE (TrqExty_EM(nmbrprtl))
        ALLOCATE (TrqExtz_EM(nmbrprtl))
        ALLOCATE (FrcTotx_EM(nmbrprtl))
        ALLOCATE (FrcToty_EM(nmbrprtl))
        ALLOCATE (FrcTotz_EM(nmbrprtl))
        ALLOCATE (TrqTotx_EM(nmbrprtl))
        ALLOCATE (TrqToty_EM(nmbrprtl))
        ALLOCATE (TrqTotz_EM(nmbrprtl))

        ALLOCATE (inepsn_EM(nmbrprtl))
        ALLOCATE (inepsk_EM(nmbrprtl))
        ALLOCATE (inmiun_EM(nmbrprtl))
        ALLOCATE (inmiuk_EM(nmbrprtl))
        ALLOCATE (ink_EM(nmbrprtl))
        ALLOCATE (ineps_EM(nmbrprtl))
        ALLOCATE (inmiu_EM(nmbrprtl))

        ALLOCATE(L1SurfIntg(nmbrprtl))
        ALLOCATE(L2SurfIntg(nmbrprtl))
        ALLOCATE(L4SurfIntg(nmbrprtl))
        ALLOCATE(L1tip(nmbrprtl))
        ALLOCATE(L2tip(nmbrprtl))
        ALLOCATE(L4tip(nmbrprtl))

        ALLOCATE (nmbrsrc_EM(0:nmbrprtl))
        ALLOCATE (srcstaID_EM(0:nmbrprtl))
        ALLOCATE (srcendID_EM(0:nmbrprtl))

        DO ithprtl = 0, nmbrprtl
            nmbrsrc_EM(ithprtl) = 0
            srcstaID_EM(ithprtl) = 0
            srcendID_EM(ithprtl) = 0
        END DO

        nmbrsrc_EM(0) = exnmbrsrc_EM

        IF (nmbrsrc_EM(0) > 0) THEN
            srcstaID_EM(0) = 1
            srcendID_EM(0) = nmbrsrc_EM(0)
        END IF

        ttlnmbrsrc_EM = nmbrsrc_EM(0)

        OPEN (71, FILE = "Input_Phys_EM.dat", STATUS = "OLD", IOSTAT = IOS)
        IF (IOS /= 0) THEN
            PRINT*, "'Input_Phys_EM.dat' does not exist! Please check!"
            STOP
        END IF

! Skip the global header read by GetPhysInputInt_EM; its length is part of the format.
        DO i = 1, 18
            READ (71, *)
        END DO

        DO ithprtl = 1, nmbrprtl

            READ (71, *)
! Boundary record: type, constant parameter, and read flag.  Common types are PEC
! (perfect electric conductor) and 2SD (two-sided dielectric interface).
            READ (71, *)
            READ (71, *) BCType_EM(ithprtl), BCValue_EM(ithprtl), BCRead_EM(ithprtl)
            READ (71, *)
! Material owned by this surface: electric n,k followed by magnetic n,k; all are
! dimensionless and relative eps and mu are formed as (n + i*kappa)**2.
            READ (71, *)
            READ (71, *) inepsn_EM(ithprtl), inepsk_EM(ithprtl), &
            &            inmiun_EM(ithprtl), inmiuk_EM(ithprtl)
            READ (71, *)

            ztp = DCMPLX(inmiun_EM(ithprtl), inmiuk_EM(ithprtl))
            inmiu_EM(ithprtl) = ztp*ztp

            ztp = DCMPLX(inepsn_EM(ithprtl), inepsk_EM(ithprtl))
            ineps_EM(ithprtl) = ztp*ztp

            ztp = ineps_EM(ithprtl)*inmiu_EM(ithprtl)
            ink_EM(ithprtl) = vcmwn_EM*CDSQRT(ztp)

! Number of prescribed sources in the domain owned by this surface.
            READ (71, *)
            READ (71, *) nmbrsrc_EM(ithprtl)
            READ (71, *)

            IF (nmbrsrc_EM(ithprtl) > 0) THEN
                srcstaID_EM(ithprtl) = ttlnmbrsrc_EM + 1
                srcendID_EM(ithprtl) = ttlnmbrsrc_EM + nmbrsrc_EM(ithprtl)
            END IF

            ttlnmbrsrc_EM = ttlnmbrsrc_EM + nmbrsrc_EM(ithprtl)

        END DO

        CLOSE (71)

    END SUBROUTINE

END MODULE
