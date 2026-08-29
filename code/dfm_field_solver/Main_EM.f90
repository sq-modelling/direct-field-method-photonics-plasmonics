! SPDX-FileCopyrightText: 2026 Qiang Sun
! SPDX-License-Identifier: BSD-3-Clause

! DFM electromagnetic examples, release 0.1.0.
! Primary software developer: Qiang Sun.
!
! Purpose: run the general direct-field boundary-element solver used for
! manuscript Figs. 1--3 and export surface or observation-domain fields.
! Case scope: PEC or dielectric particles represented by generated or imported
! linear/Q6 surface meshes; optional planar samples and angular RCS cuts are
! selected by the supplied case files.  The specialised PEC-H reconstruction
! for the lowest-frequency Fig. 2 panels is not included in this release.
!
! Inputs (current working directory): Input_Geom.dat, Input_Phys_EM.dat and,
! when impressed sources are declared, Input_Source_EM.dat.  Imported cases
! also require their Prtl_*.msh files beside these inputs.
! Principal outputs: Rslt_SurfCal{1,2,3}_EM.dat (incident, scattered and total
! surface fields), optional Rslt_Dmn*_plot.dat planar fields, optional
! Rslt_RCS_{xy,yz,zx}.dat cuts, and the Rslt_PrcdSmm.dat run record.
! Existing field/RCS data files may be appended; run each case in a fresh copy.
!
! Conventions: fields are complex phasors with exp(-i*omega*t).  This driver
! preserves the input length and field-amplitude units.  The RCS writer uses
! 4*pi*r**2*|E_sca|**2/|E_inc|**2 and reports 10*log10 of the same numeric
! value in its decibel column; no additional driver-level normalisation occurs.
!
! Main workflow: read sweep-wide controls; for each case build the geometry and
! normals, initialise materials and boundary data, solve the surface fields,
! write requested products, and release all case-local allocations.

INCLUDE './Pre_csvformat.f90'
INCLUDE './Pre_GQ_FDM.f90'
INCLUDE './Pre_Constants.f90'

INCLUDE './Geom_GlobalData.f90'
INCLUDE './Geom_Input.f90'
INCLUDE './Geom_MeshSphereCircle.f90'
INCLUDE './Geom_Mesh.f90'
INCLUDE './Geom_NormVec.f90'

INCLUDE './BRIEFGHReal.f90'
INCLUDE './BRIEFGHComp.f90'

INCLUDE './EM_SurfCal_GlobalData.f90'
INCLUDE './EM_SurfCal_Input.f90'
INCLUDE './EM_SurfCal_VariableInt.f90'
INCLUDE './EM_SurfCal_Output.f90'
INCLUDE './EM_SurfCal_PhysBC.f90'
INCLUDE './EM_SurfCal_Solver.f90'
INCLUDE './EM_DmnCal.f90'
INCLUDE './EM_DmnCal_XSec.f90'
INCLUDE './EM_DmnCal_SliceOutput.f90'

PROGRAM Main_EM

    USE omp_lib

    USE Pre_Constants
    USE Pre_csvformat

    USE Geom_GlobalData
    USE Geom_Input
    USE Geom_Mesh
    USE Geom_NormVec

    USE EM_SurfCal_GlobalData
    USE EM_SurfCal_Input
    USE EM_SurfCal_VariableInt
    USE EM_SurfCal_Output
    USE EM_SurfCal_PhysBC
    USE EM_SurfCal_Solver
    USE EM_DmnCal
    USE EM_DmnCal_XSec
    USE EM_DmnCal_SliceOutput


    IMPLICIT NONE


    LOGICAL :: filexist

    INTEGER :: ksmcs
    DOUBLE PRECISION :: tempcalcs

    INTEGER :: ithprtl,jthprtl,kthprtl,i,j,k,ii,jj,kk
    DOUBLE PRECISION :: tp,tp1,tp2,tp3,tp4,tp5,tp6,tp7,tp8,tp9, &
    &                   vcm_Wvlngth_read, tp_refractiveindx, tp_extinctioncff
    COMPLEX(KIND=KIND(1.0D0)) :: ztp,ztp1,ztp2,ztp3,ztp4,ztp5,ztp6,ztp7,ztp8,ztp9

    INTEGER, DIMENSION (8) :: LclDtTm
    CHARACTER (LEN = 12), DIMENSION (3) :: READ_REAL_CLOCK

    ! Capture the start time for the run record.
    CALL DATE_AND_TIME ( READ_REAL_CLOCK (1), READ_REAL_CLOCK (2), &
                        &READ_REAL_CLOCK (3), LclDtTm)

    OPEN (9001, FILE = "Rslt_PrcdSmm.dat", STATUS = "REPLACE")

    WRITE(9001,*) "**************************************************************"
    WRITE(9001,*)
    WRITE(9001,*) "**** This is a solver on EM problems by scalar Helm & NSBIM ****"
    WRITE(9001,*) "**** Primary developer: Qiang Sun; release 0.1.0 ****"
    WRITE(9001,*)
    WRITE(9001,*) "**************************************************************"
    WRITE(9001,*)
    WRITE(9001,*)
    WRITE(9001,*) "**************************************************************"
    WRITE(9001,*)
    WRITE(9001, "(A38,I2,A1,I2,A1,I2)") "The simulation started at (HH-MM-SS): ",&
    &                   LclDtTm(5),"-",LclDtTm(6),"-",LclDtTm(7)
    WRITE(9001, "(A17,I4,A1,I2,A1,I2,A12)") "on (YYYY-MM-DD): ",&
    &                   LclDtTm(1),"-",LclDtTm(2),"-",LclDtTm(3), " local time."
    WRITE(9001,*)
    WRITE(9001,*) "**************************************************************"
    WRITE(9001,*)

    CALL GetGeomInput_Int
    print *, 'GetGeomInput_Int OK'

    WRITE(9001,*)
    WRITE(9001,*) "**************************************************************"
    WRITE(9001,*)
    IF (MeshType == "L") WRITE(9001,*) "Linear triangle elements have been used."
    IF (MeshType == "Q") WRITE(9001,*) "Quadratic elements have been used."
    WRITE(9001,*)
    WRITE(9001,*) "**************************************************************"
    WRITE(9001,*)
    WRITE(9001,*)
    WRITE(9001,*) "**************************************************************"
    WRITE(9001,*)
    WRITE(9001,*) "The number of particles: ", nmbrprtl
    WRITE(9001,*)
    CLOSE (9001)

    CALL GetPhysInputInt_EM
    print *, 'GetPhysInputInt_EM OK'

    OPEN (9001, FILE = "Rslt_PrcdSmm.dat", STATUS = "OLD", &
        & POSITION="APPEND", ACTION="WRITE")
    WRITE(9001,*)
    WRITE(9001,*) "**************************************************************"
    WRITE(9001,*)
    WRITE(9001,*) "External relative eps: ", exeps_EM
    WRITE(9001,*) "External relative miu: ", exmiu_EM
    WRITE(9001,*)
    WRITE(9001,*) "**************************************************************"
    WRITE(9001,*)
    CLOSE (9001)

    DO ksmcs = csStartID, csEndID, csStepNumber

        tempcalcs = csStartValue + DBLE(ksmcs)*csStepSize

        CALL GetGeomInput
        print *, 'GetGeomInput OK'

        CALL Meshediting
        print *, 'Meshediting OK'

        IF (MeshType == 'L') CALL Getndnrml
        print *, 'Getndnrml OK'
        IF (MeshType == 'Q') CALL GetndnrmlQdrtcLnr
        print *, 'GetndnrmlQdrtcLnr OK'

        CALL GetPhysInput_EM
        print *, 'GetPhysInput_EM OK'

        CALL GetVariableInt_EM
        print *, 'GetVariableInt_EM OK'

        CALL GetPhysBC_EM
        print *, 'GetPhysBC_EM OK'

        OPEN (9001, FILE = "Rslt_PrcdSmm.dat", STATUS = "OLD", &
            & POSITION="APPEND", ACTION="WRITE")
        WRITE(9001,*)
        WRITE(9001,*) "**************************************************************"
        WRITE(9001,*)
        WRITE(9001,*) "Case ID: ", ksmcs
        WRITE(9001,*)
        WRITE(9001,*) "Total numbers of nodes and elements: ", ttlnmbrnd, ttlnmbrelmnt
        WRITE(9001,*)
        WRITE(9001,*) "External wave number: ", exk_EM
        WRITE(9001,*)
        CLOSE (9001)

        CALL SlvPrblm_EM
        print *, 'SlvPrblm_EM OK'

        CALL RsltSurf_fieldCase(ksmcs)
        print *, 'Surfout OK'

        IF (on_xs == 1) THEN
            CALL GetRCS_EM(ksmcs, 361, size_xs, xcen_xs,ycen_xs,zcen_xs)
            print *, 'GetRCS_EM OK'
        END IF

        IF (on_fldcal == 1) THEN
            IF (Postxyz_RMpln == 'xy') CALL GetDmn_SliceOutput_EM(ksmcs,Postxyz_RMpln)
            IF (Postxyz_RMpln == 'yz') CALL GetDmn_SliceOutput_EM(ksmcs,Postxyz_RMpln)
            IF (Postxyz_RMpln == 'zx') CALL GetDmn_SliceOutput_EM(ksmcs,Postxyz_RMpln)
            IF (Postxyz_RMpln == '3D') THEN
                 CALL GetDmn_SliceOutput_EM(ksmcs,'xy')
                 CALL GetDmn_SliceOutput_EM(ksmcs,'yz')
                 CALL GetDmn_SliceOutput_EM(ksmcs,'zx')
            END IF
            print *, 'Dmn_SliceOutput OK'
        END IF

        OPEN (9001, FILE = "Rslt_PrcdSmm.dat", STATUS = "OLD", &
            & POSITION="APPEND", ACTION="WRITE")
        WRITE(9001,*)
        WRITE(9001,*) "**************************************************************"
        WRITE(9001,*)
        WRITE(9001,*) "Case ID: ", ksmcs
        WRITE(9001,*)
        WRITE(9001,*) "External wave number: ", exk_EM
        WRITE(9001,*)
        DO ithprtl = 1, nmbrprtl
            WRITE(9001,*)
            WRITE(9001,*) "------------------------------------"
            WRITE(9001,*) "Particle ID and its core-link-shell: ", ithprtl, corelnkshell(ithprtl)
            WRITE(9001,*) "number of nodes, node ID sta, node ID end: ", &
            &              nmbrnd(ithprtl), ndstaID(ithprtl), ndendID(ithprtl)
            WRITE(9001,*) "number of elements, element ID sta, element ID end: ", &
            &              nmbrelmnt(ithprtl), elstaID(ithprtl), elendID(ithprtl)
            WRITE(9001,*) "Particle surface area: ", surfarea(ithprtl)
            WRITE(9001,*) "Particle volume: ", volume(ithprtl)
            WRITE(9001,*) "Particle location offset: ", xloctn(ithprtl),yloctn(ithprtl),zloctn(ithprtl)
            WRITE(9001,*) "Particle orientation: ", anglecal_x(ithprtl),anglecal_y(ithprtl),anglecal_z(ithprtl)
            WRITE(9001,*) "Particle centre of mass: ", xmssctr(ithprtl),ymssctr(ithprtl),zmssctr(ithprtl)
            WRITE(9001,*)
            WRITE(9001,*) "Internal relative eps: ", ineps_EM(ithprtl)
            WRITE(9001,*) "Internal relative miu: ", inmiu_EM(ithprtl)
            WRITE(9001,*) "Internal wave number: ", ink_EM(ithprtl)
            WRITE(9001,*) "------------------------------------"
            WRITE(9001,*)
        END DO
        WRITE(9001,*) "**************************************************************"
        WRITE(9001,*)
        CLOSE (9001)

        DEALLOCATE (PrtlType,MeshRead,NrmlInOut,MeshGnrtn,Meshnlvl,MeshRlvlStp, &
        &           xnrmref,ynrmref,znrmref)
        DEALLOCATE (sizezoom,xloctn,yloctn,zloctn)
        DEALLOCATE (corelnkshell)
        DEALLOCATE (bvsa,cvsa,dvsa,bowl_a,bowl_b, &
        &           dfsp_a,dfsp_b,dfsp_c,dfsp_l,dfsp_m,dfsp_n, &
        &           anglecal_x,anglecal_y,anglecal_z,surfarea,volume)
        DEALLOCATE (nmbrnd,nmbrelmnt,ndstaID,ndendID,elstaID,elendID)
        DEALLOCATE (rho_den,xmssctr,ymssctr,zmssctr)

        DEALLOCATE (xnd,ynd,znd,nnx,nny,nnz,t1x,t1y,t1z,t2x,t2y,t2z)
        DEALLOCATE (curvt1,curvt2,curvmn,curvt1th,curvt2th,curvmnth)
        DEALLOCATE (elmntarea,nnxelmnt,nnyelmnt,nnzelmnt)
        DEALLOCATE (srcfmm_vec,srcfmm_nrm,srcfmm_wght,srcfmm_wtnd)
        DEALLOCATE (elmntlnknd,ndlnkelmnt, &
        &           ndlnknd1st,ndlnknd1stslf,ndlnknd2nd,ndlnknd2ndslf)
        IF (MeshType == "Q") DEALLOCATE (elmntlnkndlnr,ndlnkelmntlnr)
        DEALLOCATE (d_dt1,d_dt2,d2_dt1,d2_dt2)

        DEALLOCATE (BCType_EM,BCValue_EM,BCRead_EM)
        DEALLOCATE (inepsn_EM,inepsk_EM,inmiun_EM,inmiuk_EM,ink_EM,ineps_EM,inmiu_EM)
        DEALLOCATE (exsurfQ1_EM,exsurfQ2_EM,exsurfQ3_EM, &
        &           insurfQ1_EM,insurfQ2_EM,insurfQ3_EM  )
        DEALLOCATE (Frcx_EM,Frcy_EM,Frcz_EM,Trqx_EM,Trqy_EM,Trqz_EM, &
        &           FrcIncx_EM, FrcIncy_EM, FrcIncz_EM, TrqIncx_EM, TrqIncy_EM, TrqIncz_EM, &
        &           FrcScax_EM, FrcScay_EM, FrcScaz_EM, TrqScax_EM, TrqScay_EM, TrqScaz_EM, &
        &           FrcExtx_EM, FrcExty_EM, FrcExtz_EM, TrqExtx_EM, TrqExty_EM, TrqExtz_EM, &
        &           FrcTotx_EM, FrcToty_EM, FrcTotz_EM, TrqTotx_EM, TrqToty_EM, TrqTotz_EM, &
        &           FrcElcx_EM, FrcElcy_EM, FrcElcz_EM, TrqElcx_EM, TrqElcy_EM, TrqElcz_EM, &
        &           FrcMagx_EM, FrcMagy_EM, FrcMagz_EM, TrqMagx_EM, TrqMagy_EM, TrqMagz_EM, &
        &           FrcTtlx_EM, FrcTtly_EM, FrcTtlz_EM, TrqTtlx_EM, TrqTtly_EM, TrqTtlz_EM  )
        DEALLOCATE (L1SurfIntg,L2SurfIntg,L4SurfIntg,L1tip,L2tip,L4tip)
        DEALLOCATE (nmbrsrc_EM,srcstaID_EM,srcendID_EM,srcType_EM,srcStrength_EM, &
        &           xsrc_EM,ysrc_EM,zsrc_EM,polxsrc_EM,polysrc_EM,polzsrc_EM )
        DEALLOCATE (&
        &   exE1x_EM, exE1y_EM, exE1z_EM, inE1x_EM, inE1y_EM, inE1z_EM, &
        &   exE1xdnn_EM, exE1ydnn_EM, exE1zdnn_EM, inE1xdnn_EM, inE1ydnn_EM, inE1zdnn_EM, &
        &   exE1xdt1_EM, exE1ydt1_EM, exE1zdt1_EM, inE1xdt1_EM, inE1ydt1_EM, inE1zdt1_EM, &
        &   exE1xdt2_EM, exE1ydt2_EM, exE1zdt2_EM, inE1xdt2_EM, inE1ydt2_EM, inE1zdt2_EM, &
        &   exE2x_EM, exE2y_EM, exE2z_EM, inE2x_EM, inE2y_EM, inE2z_EM, &
        &   exE2xdnn_EM, exE2ydnn_EM, exE2zdnn_EM, inE2xdnn_EM, inE2ydnn_EM, inE2zdnn_EM, &
        &   exE3x_EM, exE3y_EM, exE3z_EM, inE3x_EM, inE3y_EM, inE3z_EM, &
        &   exE3xdnn_EM, exE3ydnn_EM, exE3zdnn_EM, inE3xdnn_EM, inE3ydnn_EM, inE3zdnn_EM, &
        &   &
        &   exH1x_EM, exH1y_EM, exH1z_EM, inH1x_EM, inH1y_EM, inH1z_EM, &
        &   exH1xdnn_EM, exH1ydnn_EM, exH1zdnn_EM, inH1xdnn_EM, inH1ydnn_EM, inH1zdnn_EM, &
        &   exH1xdt1_EM, exH1ydt1_EM, exH1zdt1_EM, inH1xdt1_EM, inH1ydt1_EM, inH1zdt1_EM, &
        &   exH1xdt2_EM, exH1ydt2_EM, exH1zdt2_EM, inH1xdt2_EM, inH1ydt2_EM, inH1zdt2_EM, &
        &   exH2x_EM, exH2y_EM, exH2z_EM, inH2x_EM, inH2y_EM, inH2z_EM, &
        &   exH2xdnn_EM, exH2ydnn_EM, exH2zdnn_EM, inH2xdnn_EM, inH2ydnn_EM, inH2zdnn_EM, &
        &   exH3x_EM, exH3y_EM, exH3z_EM, inH3x_EM, inH3y_EM, inH3z_EM, &
        &   exH3xdnn_EM, exH3ydnn_EM, exH3zdnn_EM, inH3xdnn_EM, inH3ydnn_EM, inH3zdnn_EM  )

    END DO

    DEALLOCATE (wg_glqln1d,xg_glqln1d,wg_glqtr2d,xg_glqtr2d,yg_glqtr2d)

    CALL DATE_AND_TIME (READ_REAL_CLOCK (1), READ_REAL_CLOCK (2), &
    &                   READ_REAL_CLOCK (3), LclDtTm)

    OPEN (9001, FILE = "Rslt_PrcdSmm.dat", STATUS = "OLD", &
            & POSITION="APPEND", ACTION="WRITE")
    WRITE(9001, "(A38,I2,A1,I2,A1,I2)") "The simulation ended at (HH-MM-SS): ",&
    &                   LclDtTm(5),"-",LclDtTm(6),"-",LclDtTm(7)
    WRITE(9001, "(A17,I4,A1,I2,A1,I2,A12)") "on (YYYY-MM-DD): ",&
    &                   LclDtTm(1),"-",LclDtTm(2),"-",LclDtTm(3), " local time."
    WRITE(9001,*)
    WRITE(9001,*) "**************************************************************"
    WRITE(9001,*)
    WRITE(9001,*)
    WRITE(9001,*)
    CLOSE(9001)

END PROGRAM
