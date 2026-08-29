! SPDX-FileCopyrightText: 2026 Qiang Sun
! SPDX-License-Identifier: BSD-3-Clause

! DFM electromagnetic examples, release 0.1.0.
! Primary software developer: Qiang Sun.
!
! Purpose: compute the orientation-dependent optical force and torque used for
! manuscript Fig. 6.
! Case scope: exactly two externally meshed particles, a single-frequency
! steady-state solve, and rotation of particle 1 about the y axis.  The supplied
! inputs generate 73 orientations from 0 to 360 degrees in 5-degree steps.
! Pulse/FFT, field-slice, enhancement, cross-section and RCS branches are not
! part of this release driver.
!
! Inputs (current working directory): Input_Geom.dat, Input_Phys_EM.dat and the
! imported Prtl_0001.msh and Prtl_0002.msh surfaces.  Input_Source_EM.dat is
! read only if impressed sources are declared; the supplied case has none.
! Principal outputs: Rslt_FrcTrq.dat contains angle followed by three force and
! three torque components for each particle; Rslt_PrcdSmm.dat is a run record.
!
! Conventions: fields are complex phasors with exp(-i*omega*t), and theta is in
! degrees.  For the supplied unit-amplitude incident field, output forces are
! F/(epsilon_0*a**2) and torques are N/(epsilon_0*a**3), where
! a=sizezoom(1).  GetFrcTrq_EM itself returns unnormalised surface integrals.
!
! Main workflow: read the fixed controls; for each angle rotate particle 1,
! import and transform both meshes, construct normals and boundary data, solve
! the surface fields, integrate stress-tensor traction, write one row, and
! release all case-local allocations.

INCLUDE './Pre_csvformat.f90'
INCLUDE './Pre_GQ_FDM.f90'
INCLUDE './Pre_Constants.f90'

INCLUDE './Geom_GlobalData.f90'
INCLUDE './Geom_Input.f90'
INCLUDE './Geom_MeshSphereCircle.f90'
INCLUDE './Geom_Mesh.f90'
INCLUDE './Geom_NormVec.f90'

INCLUDE './BRIEFGHComp.f90'

INCLUDE './EM_SurfCal_GlobalData.f90'
INCLUDE './EM_SurfCal_Input.f90'
INCLUDE './EM_SurfCal_VariableInt.f90'
INCLUDE './EM_SurfCal_PhysBC.f90'
INCLUDE './EM_SurfCal_Solver.f90'
INCLUDE './EM_SurfCal_ForceTorque.f90'


PROGRAM Main_EM

    USE Pre_Constants
    USE Pre_csvformat
    USE Geom_GlobalData
    USE Geom_Input
    USE Geom_Mesh
    USE Geom_NormVec
    USE EM_SurfCal_GlobalData
    USE EM_SurfCal_Input
    USE EM_SurfCal_VariableInt
    USE EM_SurfCal_PhysBC
    USE EM_SurfCal_Solver
    USE EM_SurfCal_ForceTorque

    IMPLICIT NONE

    INTEGER :: ksmcs
    DOUBLE PRECISION :: tempcalcs, force_scale, torque_scale

    OPEN(9001,FILE='Rslt_PrcdSmm.dat',STATUS='REPLACE')
    WRITE(9001,'(A)') 'DFM optical-force example'
    WRITE(9001,'(A)') 'Primary developer: Qiang Sun; release 0.1.0'
    CLOSE(9001)

    CALL GetGeomInput_Int
    CALL GetPhysInputInt_EM

    IF (FFTpulse_EM /= 0) THEN
        WRITE(*,'(A)') 'This release driver supports only FFTpulse_EM=0.'
        STOP 2
    END IF
    IF (on_enhance /= 0 .OR. on_xs /= 0 .OR. on_RCS /= 0 .OR. on_fldcal /= 0) THEN
        WRITE(*,'(A)') 'This release driver is restricted to force/torque output.'
        STOP 3
    END IF

    OPEN(21,FILE='Rslt_FrcTrq.dat',STATUS='REPLACE')
    CALL csv_write_char(21,'Variables = "theta", ',.false.,'spc')
    CALL csv_write_char(21,'"Fx_1", ',.false.,'spc')
    CALL csv_write_char(21,'"Fy_1", ',.false.,'spc')
    CALL csv_write_char(21,'"Fz_1", ',.false.,'spc')
    CALL csv_write_char(21,'"Nx_1", ',.false.,'spc')
    CALL csv_write_char(21,'"Ny_1", ',.false.,'spc')
    CALL csv_write_char(21,'"Nz_1", ',.false.,'spc')
    CALL csv_write_char(21,'"Fx_2", ',.false.,'spc')
    CALL csv_write_char(21,'"Fy_2", ',.false.,'spc')
    CALL csv_write_char(21,'"Fz_2", ',.false.,'spc')
    CALL csv_write_char(21,'"Nx_2", ',.false.,'spc')
    CALL csv_write_char(21,'"Ny_2", ',.false.,'spc')
    CALL csv_write_char(21,'"Nz_2"',.true.,'spc')
    CLOSE(21)

    ! Each sweep point is an independent mesh transformation and field solve.
    DO ksmcs = csStartID, csEndID, csStepNumber
        tempcalcs = csStartValue + DBLE(ksmcs)*csStepSize

        CALL GetGeomInput
        IF (nmbrprtl /= 2 .OR. ANY(MeshRead /= 1)) THEN
            WRITE(*,'(A)') 'Fig. 6 requires two externally meshed particles.'
            STOP 4
        END IF

        anglecal_y(1) = tempcalcs/180.0D0*pai
        CALL Meshediting
        IF (MeshType == 'L') CALL Getndnrml
        IF (MeshType == 'Q') CALL GetndnrmlQdrtcLnr

        CALL GetPhysInput_EM
        CALL GetVariableInt_EM
        CALL GetPhysBC_EM
        CALL SlvPrblm_EM
        CALL GetFrcTrq_EM

        force_scale = vcm_eps0*sizezoom(1)**2
        torque_scale = force_scale*sizezoom(1)
        OPEN(21,FILE='Rslt_FrcTrq.dat',STATUS='OLD',POSITION='APPEND',ACTION='WRITE')
        CALL csv_write_dble(21,tempcalcs,.false.,'cmr')
        CALL csv_write_dble(21,Frcx_EM(1)/force_scale,.false.,'cmr')
        CALL csv_write_dble(21,Frcy_EM(1)/force_scale,.false.,'cmr')
        CALL csv_write_dble(21,Frcz_EM(1)/force_scale,.false.,'cmr')
        CALL csv_write_dble(21,Trqx_EM(1)/torque_scale,.false.,'cmr')
        CALL csv_write_dble(21,Trqy_EM(1)/torque_scale,.false.,'cmr')
        CALL csv_write_dble(21,Trqz_EM(1)/torque_scale,.false.,'cmr')
        CALL csv_write_dble(21,Frcx_EM(2)/force_scale,.false.,'cmr')
        CALL csv_write_dble(21,Frcy_EM(2)/force_scale,.false.,'cmr')
        CALL csv_write_dble(21,Frcz_EM(2)/force_scale,.false.,'cmr')
        CALL csv_write_dble(21,Trqx_EM(2)/torque_scale,.false.,'cmr')
        CALL csv_write_dble(21,Trqy_EM(2)/torque_scale,.false.,'cmr')
        CALL csv_write_dble(21,Trqz_EM(2)/torque_scale,.true.,'spc')
        CLOSE(21)

        CALL ReleaseCaseData
    END DO

    DEALLOCATE(wg_glqln1d,xg_glqln1d,wg_glqtr2d,xg_glqtr2d,yg_glqtr2d)

CONTAINS

    ! Contract: after one fully initialised and solved sweep point, deallocate
    ! every case-local geometry and EM array exactly once.  Shared quadrature
    ! tables remain allocated until the driver exits the sweep.
    SUBROUTINE ReleaseCaseData
        DEALLOCATE(PrtlType,MeshRead,NrmlInOut,MeshGnrtn,Meshnlvl,MeshRlvlStp, &
        &          xnrmref,ynrmref,znrmref)
        DEALLOCATE(sizezoom,xloctn,yloctn,zloctn,corelnkshell)
        DEALLOCATE(bvsa,cvsa,dvsa,bowl_a,bowl_b,Cheb_a,Cheb_b,GRSp_a,GRSp_b, &
        &          dfsp_a,dfsp_b,dfsp_c,dfsp_l,dfsp_m,dfsp_n, &
        &          anglecal_x,anglecal_y,anglecal_z,surfarea,volume)
        DEALLOCATE(nmbrnd,nmbrelmnt,ndstaID,ndendID,elstaID,elendID)
        DEALLOCATE(rho_den,xmssctr,ymssctr,zmssctr)

        DEALLOCATE(xnd,ynd,znd,nnx,nny,nnz,t1x,t1y,t1z,t2x,t2y,t2z)
        DEALLOCATE(curvt1,curvt2,curvmn,curvt1th,curvt2th,curvmnth)
        DEALLOCATE(elmntarea,nnxelmnt,nnyelmnt,nnzelmnt)
        DEALLOCATE(srcfmm_vec,srcfmm_nrm,srcfmm_wght,srcfmm_wtnd)
        DEALLOCATE(elmntlnknd,ndlnkelmnt,ndlnknd1st,ndlnknd1stslf, &
        &          ndlnknd2nd,ndlnknd2ndslf)
        IF (MeshType == 'Q') DEALLOCATE(elmntlnkndlnr,ndlnkelmntlnr)
        DEALLOCATE(d_dt1,d_dt2,d2_dt1,d2_dt2)
        DEALLOCATE(dt1dt1_t1,dt1dt2_t1,dt1dt1_t2,dt1dt2_t2,dt1dt1_nn,dt1dt2_nn, &
        &          dt2dt1_t1,dt2dt2_t1,dt2dt1_t2,dt2dt2_t2,dt2dt1_nn,dt2dt2_nn, &
        &          dnndt1_t1,dnndt2_t1,dnndt1_t2,dnndt2_t2,dnndt1_nn,dnndt2_nn)

        DEALLOCATE(BCType_EM,BCValue_EM,BCRead_EM)
        DEALLOCATE(inepsn_EM,inepsk_EM,inmiun_EM,inmiuk_EM,ink_EM,ineps_EM,inmiu_EM)
        DEALLOCATE(exsurfQ1_EM,exsurfQ2_EM,exsurfQ3_EM,insurfQ1_EM,insurfQ2_EM,insurfQ3_EM)
        DEALLOCATE(Frcx_EM,Frcy_EM,Frcz_EM,Trqx_EM,Trqy_EM,Trqz_EM, &
        &          FrcIncx_EM,FrcIncy_EM,FrcIncz_EM,TrqIncx_EM,TrqIncy_EM,TrqIncz_EM, &
        &          FrcScax_EM,FrcScay_EM,FrcScaz_EM,TrqScax_EM,TrqScay_EM,TrqScaz_EM, &
        &          FrcExtx_EM,FrcExty_EM,FrcExtz_EM,TrqExtx_EM,TrqExty_EM,TrqExtz_EM, &
        &          FrcTotx_EM,FrcToty_EM,FrcTotz_EM,TrqTotx_EM,TrqToty_EM,TrqTotz_EM, &
        &          FrcElcx_EM,FrcElcy_EM,FrcElcz_EM,TrqElcx_EM,TrqElcy_EM,TrqElcz_EM, &
        &          FrcMagx_EM,FrcMagy_EM,FrcMagz_EM,TrqMagx_EM,TrqMagy_EM,TrqMagz_EM, &
        &          FrcTtlx_EM,FrcTtly_EM,FrcTtlz_EM,TrqTtlx_EM,TrqTtly_EM,TrqTtlz_EM)
        DEALLOCATE(L1SurfIntg,L2SurfIntg,L4SurfIntg,L1tip,L2tip,L4tip)
        DEALLOCATE(nmbrsrc_EM,srcstaID_EM,srcendID_EM,srcType_EM, &
        &          xsrc_EM,ysrc_EM,zsrc_EM,polxsrc_EM,polysrc_EM,polzsrc_EM)
        DEALLOCATE(exE1x_EM,exE1y_EM,exE1z_EM,inE1x_EM,inE1y_EM,inE1z_EM, &
        &          exE1xdnn_EM,exE1ydnn_EM,exE1zdnn_EM,inE1xdnn_EM,inE1ydnn_EM,inE1zdnn_EM, &
        &          exE1xdt1_EM,exE1ydt1_EM,exE1zdt1_EM,inE1xdt1_EM,inE1ydt1_EM,inE1zdt1_EM, &
        &          exE1xdt2_EM,exE1ydt2_EM,exE1zdt2_EM,inE1xdt2_EM,inE1ydt2_EM,inE1zdt2_EM, &
        &          exE2x_EM,exE2y_EM,exE2z_EM,inE2x_EM,inE2y_EM,inE2z_EM, &
        &          exE2xdnn_EM,exE2ydnn_EM,exE2zdnn_EM,inE2xdnn_EM,inE2ydnn_EM,inE2zdnn_EM, &
        &          exE3x_EM,exE3y_EM,exE3z_EM,inE3x_EM,inE3y_EM,inE3z_EM, &
        &          exE3xdnn_EM,exE3ydnn_EM,exE3zdnn_EM,inE3xdnn_EM,inE3ydnn_EM,inE3zdnn_EM, &
        &          exH1x_EM,exH1y_EM,exH1z_EM,inH1x_EM,inH1y_EM,inH1z_EM, &
        &          exH1xdnn_EM,exH1ydnn_EM,exH1zdnn_EM,inH1xdnn_EM,inH1ydnn_EM,inH1zdnn_EM, &
        &          exH1xdt1_EM,exH1ydt1_EM,exH1zdt1_EM,inH1xdt1_EM,inH1ydt1_EM,inH1zdt1_EM, &
        &          exH1xdt2_EM,exH1ydt2_EM,exH1zdt2_EM,inH1xdt2_EM,inH1ydt2_EM,inH1zdt2_EM, &
        &          exH2x_EM,exH2y_EM,exH2z_EM,inH2x_EM,inH2y_EM,inH2z_EM, &
        &          exH2xdnn_EM,exH2ydnn_EM,exH2zdnn_EM,inH2xdnn_EM,inH2ydnn_EM,inH2zdnn_EM, &
        &          exH3x_EM,exH3y_EM,exH3z_EM,inH3x_EM,inH3y_EM,inH3z_EM, &
        &          exH3xdnn_EM,exH3ydnn_EM,exH3zdnn_EM,inH3xdnn_EM,inH3ydnn_EM,inH3zdnn_EM)
    END SUBROUTINE ReleaseCaseData

END PROGRAM Main_EM
