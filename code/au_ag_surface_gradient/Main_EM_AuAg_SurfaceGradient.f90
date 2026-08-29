! SPDX-FileCopyrightText: 2026 Qiang Sun
! SPDX-License-Identifier: BSD-3-Clause

! DFM electromagnetic examples, release 0.1.0.
! Primary software developer: Qiang Sun.
!
! Purpose: solve the fixed eccentric Au-core/Ag-shell case used for manuscript
! Fig. 5 and export the two-sided total-field traces needed by finite-distance
! surface-gradient postprocessing.
! Case scope: two generated Q6 interfaces at one common mesh level, 480 nm
! vacuum wavelength, a 50 nm Ag outer radius, and a 40 nm Au core displaced by
! 5 nm along +x.  ValidateCase enforces the complete release-case definition.
!
! Inputs (current working directory): Input_Geom.dat and Input_Phys_EM.dat.
! The supplied case declares no impressed sources.  The executable accepts no
! argument or the optional compatibility token "run".
! Principal outputs: surface_traces_both_sides.dat, surface_elements_q6.dat,
! cross_sections.dat, surface_trace_audit_summary.txt and Rslt_PrcdSmm.dat.
!
! Conventions: solver geometry is in um and fields use exp(-i*omega*t).  The
! trace exporter converts coordinates to m, curvature to 1/m, H to A/m, and
! normal derivatives to per-metre SI form; E is in V/m.  Cross sections remain
! in um**2.  Each exported interface normal points toward Ag, and both sides'
! derivatives are reported with respect to that same interface normal.
!
! Main workflow: read and validate the fixed case, generate geometry and
! normals, initialise and solve the two-domain surface system, integrate the
! enclosing-sphere cross sections, then write traces, connectivity and audits.

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
INCLUDE './EM_DmnCal.f90'
INCLUDE './EM_DmnCal_XSec.f90'


PROGRAM Main_EM_AuAg_SurfaceGradient

    USE Pre_Constants
    USE Geom_GlobalData
    USE Geom_Input
    USE Geom_Mesh
    USE Geom_NormVec
    USE EM_SurfCal_GlobalData
    USE EM_SurfCal_Input
    USE EM_SurfCal_VariableInt
    USE EM_SurfCal_PhysBC
    USE EM_SurfCal_Solver
    USE EM_DmnCal_XSec

    IMPLICIT NONE

    INTEGER :: level
    CHARACTER(LEN=32) :: arg
    DOUBLE PRECISION :: sigma_inc,sigma_sca,sigma_ext,sigma_abs
    DOUBLE PRECISION, PARAMETER :: um_to_m=1.0D-6

    CALL GetGeomInput_Int
    CALL GetPhysInputInt_EM

    IF (csStartID /= csEndID) CALL Die('Exactly one geometry case is required.')
    IF (nmbrprtl /= 2) CALL Die('Exactly two interfaces are required.')

    CALL GetGeomInput
    level = Meshnlvl(1)
    IF (Meshnlvl(2) /= level) CALL Die('Both interfaces must use the same Q6 level.')
    CALL GET_COMMAND_ARGUMENT(1,arg)
    IF (LEN_TRIM(arg) > 0) THEN
        IF (TRIM(arg) /= 'run') CALL Die('Usage: au_ag_surface_gradient [run]')
    END IF

    CALL Meshediting
    IF (MeshType /= 'Q') CALL Die('The exporter requires Q6 elements.')
    CALL GetndnrmlQdrtcLnr
    CALL GetPhysInput_EM

    CALL ValidateCase(level)
    CALL GetVariableInt_EM
    CALL GetPhysBC_EM

    WRITE(*,'(A,I0,A,I0)') 'Solving Au-core/Ag-shell, Q6 level ',level,', nodes=',ttlnmbrnd
    CALL SlvPrblm_EM

    OPEN(88,FILE='Rslt_PrcdSmm.dat',STATUS='REPLACE')
    WRITE(88,'(A)') 'Au-core/Ag-shell surface-gradient run'
    CLOSE(88)
    CALL GetCrossSection_EM(size_xs,xcen_xs,ycen_xs,zcen_xs, &
    &                        sigma_inc,sigma_sca,sigma_ext,sigma_abs)

    CALL WriteSurfaceTraces
    CALL WriteElementsQ6
    CALL WriteCrossSections(sigma_sca,sigma_ext,sigma_abs)
    CALL WriteAudit(level,sigma_sca,sigma_ext,sigma_abs)

    WRITE(*,'(A,I0,A,I0)') 'Export complete: nodes=',ttlnmbrnd,', Q6 elements=',ttlnmbrelmnt

CONTAINS

    ! Contract: print MESSAGE and terminate with status 2; this routine does
    ! not return to its caller.
    SUBROUTINE Die(message)
        CHARACTER(LEN=*), INTENT(IN) :: message
        WRITE(*,'(A)') TRIM(message)
        STOP 2
    END SUBROUTINE Die

    ! Contract: after geometry, normals and material data are initialised,
    ! verify every invariant of the fixed Fig. 5 case or terminate via Die.
    SUBROUTINE ValidateCase(mesh_level)
        INTEGER, INTENT(IN) :: mesh_level
        DOUBLE PRECISION :: expected_k

        expected_k = 2.0D0*pai/0.480D0*1.33D0
        IF (mesh_level < 2) CALL Die('Mesh level must be at least 2.')
        IF (ANY(PrtlType /= 'Sphr')) CALL Die('Both interfaces must be generated from spheres.')
        IF (ANY(MeshRead /= 0)) CALL Die('Expected generated meshes, not externally read meshes.')
        IF (ANY(MeshGnrtn /= 'Icshdrl')) CALL Die('Both sphere meshes must use Icshdrl.')
        IF (ANY(ABS(MeshRlvlStp-1.0D0) > 1.0D-12)) CALL Die('Expected uniform Icshdrl mesh spacing.')
        IF (excitetype_EM /= 'pwe') CALL Die('Expected pwe excitation.')
        IF (poltype_EM /= 'p') CALL Die('Expected canonical plane-wave polarization selector p.')
        IF (incOrder_EM /= 1 .OR. ABS(incFeature_EM) > 1.0D-12) &
        &   CALL Die('Expected zero-phase canonical plane-wave component ordering.')
        IF (ABS(vcmwl_EM-0.480D0) > 1.0D-12) CALL Die('Expected lambda0=0.480 um.')
        IF (ABS(exepsn_EM-1.33D0) > 1.0D-12 .OR. ABS(exepsk_EM) > 1.0D-12) &
        &   CALL Die('Expected exterior electric n+ik = 1.33 + i0.')
        IF (ABS(exmiun_EM-1.0D0) > 1.0D-12 .OR. ABS(exmiuk_EM) > 1.0D-12) &
        &   CALL Die('Expected exterior magnetic n+ik = 1 + i0.')
        IF (ABS(inepsn_EM(1)-0.13734195D0) > 1.0D-12 .OR. &
        &   ABS(inepsk_EM(1)-2.63758568D0) > 1.0D-12) &
        &   CALL Die('Expected Ag electric n+ik at 480 nm.')
        IF (ABS(inepsn_EM(2)-1.09553201D0) > 1.0D-12 .OR. &
        &   ABS(inepsk_EM(2)-1.76334636D0) > 1.0D-12) &
        &   CALL Die('Expected Au electric n+ik at 480 nm.')
        IF (ANY(ABS(inmiun_EM-1.0D0) > 1.0D-12) .OR. ANY(ABS(inmiuk_EM) > 1.0D-12)) &
        &   CALL Die('Expected unit magnetic n+ik in Ag and Au.')
        IF (ABS(DBLE(exk_EM)-expected_k) > 1.0D-10 .OR. ABS(DIMAG(exk_EM)) > 1.0D-12) &
        &   CALL Die('Unexpected exterior wave number.')
        IF (MAXVAL(ABS((/inckx_EM,incky_EM,inckz_EM/)-(/0.0D0,0.0D0,1.0D0/))) > 1.0D-12) &
        &   CALL Die('Expected propagation along +z.')
        IF (MAXVAL(ABS((/incFieldx_EM,incFieldy_EM,incFieldz_EM/)-(/1.0D0,0.0D0,0.0D0/))) > 1.0D-12) &
        &   CALL Die('Expected x-polarized unit incident E field.')
        IF (ABS(incFieldmdl_EM-1.0D0) > 1.0D-12) CALL Die('Expected incident E amplitude 1 V/m.')
        IF (on_xs /= 1 .OR. ABS(size_xs-1.0D0) > 1.0D-12 .OR. &
        &   MAXVAL(ABS((/xcen_xs,ycen_xs,zcen_xs/))) > 1.0D-12) &
        &   CALL Die('Expected cross-section integration sphere radius 1 um at the origin.')
        IF (exnmbrsrc_EM /= 0 .OR. ANY(nmbrsrc_EM /= 0)) CALL Die('Expected no impressed sources.')
        IF (ANY(BCType_EM /= '2SD')) CALL Die('Both interfaces must use 2SD.')
        IF (corelnkshell(1) /= 0 .OR. corelnkshell(2) /= 1) CALL Die('Unexpected core-shell topology.')
        IF (NrmlInOut(1) /= 1 .OR. NrmlInOut(2) /= -1) CALL Die('Unexpected code-normal orientation.')
        IF (ABS(sizezoom(1)-0.050D0) > 1.0D-12 .OR. ABS(sizezoom(2)-0.040D0) > 1.0D-12) &
        &   CALL Die('Unexpected radii.')
        IF (ABS(xloctn(1)) > 1.0D-12 .OR. ABS(xloctn(2)-0.005D0) > 1.0D-12 .OR. &
        &   ABS(yloctn(1))+ABS(zloctn(1))+ABS(yloctn(2))+ABS(zloctn(2)) > 1.0D-12) &
        &   CALL Die('Unexpected sphere centres.')
    END SUBROUTINE ValidateCase

    ! Contract: map a one-based global node ID to interface 1 (water/Ag),
    ! interface 2 (Ag/Au), or 0 when the ID lies outside both node ranges.
    INTEGER FUNCTION NodeInterface(node_id)
        INTEGER, INTENT(IN) :: node_id
        IF (node_id >= ndstaID(1) .AND. node_id <= ndendID(1)) THEN
            NodeInterface=1
        ELSE IF (node_id >= ndstaID(2) .AND. node_id <= ndendID(2)) THEN
            NodeInterface=2
        ELSE
            NodeInterface=0
        END IF
    END FUNCTION NodeInterface

    ! Contract: map a one-based global element ID to interface 1 (water/Ag),
    ! interface 2 (Ag/Au), or 0 when the ID lies outside both element ranges.
    INTEGER FUNCTION ElementInterface(element_id)
        INTEGER, INTENT(IN) :: element_id
        IF (element_id >= elstaID(1) .AND. element_id <= elendID(1)) THEN
            ElementInterface=1
        ELSE IF (element_id >= elstaID(2) .AND. element_id <= elendID(2)) THEN
            ElementInterface=2
        ELSE
            ElementInterface=0
        END IF
    END FUNCTION ElementInterface

    ! Contract: append SCALE*VECTOR as three real/imaginary pairs after OFFSET
    ! in VALUES, advancing OFFSET by six; the caller provides sufficient space.
    SUBROUTINE AppendComplex3(values,offset,vector,scale)
        DOUBLE PRECISION, INTENT(INOUT) :: values(63)
        INTEGER, INTENT(INOUT) :: offset
        COMPLEX(KIND=KIND(1.0D0)), INTENT(IN) :: vector(3)
        DOUBLE PRECISION, INTENT(IN) :: scale
        INTEGER :: component
        DO component=1,3
            offset=offset+1
            values(offset)=DBLE(vector(component))*scale
            offset=offset+1
            values(offset)=DIMAG(vector(component))*scale
        END DO
    END SUBROUTINE AppendComplex3

    ! Contract: after SlvPrblm_EM, replace surface_traces_both_sides.dat with
    ! one SI-scaled, two-sided complex trace record per global surface node.
    ! The exported kappa1 and kappa2 columns are directional curvatures along
    ! t1 and t2, not principal curvatures; kappa_sum is their signed trace.
    SUBROUTINE WriteSurfaceTraces
        INTEGER :: i,interface_id,outer_region,inner_region,offset
        DOUBLE PRECISION :: values(63),normal_sign,h_scale,dn_sign
        COMPLEX(KIND=KIND(1.0D0)) :: vector(3)
        DOUBLE PRECISION, PARAMETER :: z0_si=376.730313668D0

        OPEN(91,FILE='surface_traces_both_sides.dat',STATUS='REPLACE')
        WRITE(91,'(A)') '# schema: auag-two-sided-q6-v1; numeric rows only; all fields are total complex phasors.'
        WRITE(91,'(A)') '# normal: interface-specific and always points toward Ag; both d/dn traces use that same normal; SI units.'
        WRITE(91,'(A)') '# columns: node_id interface_id ex_region_id in_region_id ' // &
        & 'x_m y_m z_m nx ny nz t1x t1y t1z t2x t2y t2z ' // &
        & 'kappa1_per_m kappa2_per_m kappa_sum_per_m ' // &
        & 'Ex_ex_re Ex_ex_im Ey_ex_re Ey_ex_im Ez_ex_re Ez_ex_im ' // &
        & 'Hx_ex_re Hx_ex_im Hy_ex_re Hy_ex_im Hz_ex_re Hz_ex_im ' // &
        & 'dEx_ex_dn_re dEx_ex_dn_im dEy_ex_dn_re dEy_ex_dn_im ' // &
        & 'dEz_ex_dn_re dEz_ex_dn_im dHx_ex_dn_re dHx_ex_dn_im ' // &
        & 'dHy_ex_dn_re dHy_ex_dn_im dHz_ex_dn_re dHz_ex_dn_im ' // &
        & 'Ex_in_re Ex_in_im Ey_in_re Ey_in_im Ez_in_re Ez_in_im ' // &
        & 'Hx_in_re Hx_in_im Hy_in_re Hy_in_im Hz_in_re Hz_in_im ' // &
        & 'dEx_in_dn_re dEx_in_dn_im dEy_in_dn_re dEy_in_dn_im ' // &
        & 'dEz_in_dn_re dEz_in_dn_im dHx_in_dn_re dHx_in_dn_im ' // &
        & 'dHy_in_dn_re dHy_in_dn_im dHz_in_dn_re dHz_in_dn_im'

        h_scale=1.0D0/z0_si
        DO i=1,ttlnmbrnd
            interface_id=NodeInterface(i)
            IF (interface_id == 1) THEN
                outer_region=1
                inner_region=2
                normal_sign=1.0D0
            ELSE IF (interface_id == 2) THEN
                outer_region=2
                inner_region=3
                ! The native inner-interface normal already points from Au to Ag.
                normal_sign=1.0D0
            ELSE
                CALL Die('Node is not assigned to an interface.')
            END IF
            dn_sign=normal_sign/um_to_m

            values=0.0D0
            values(1:3)=(/xnd(i),ynd(i),znd(i)/)*um_to_m
            values(4:6)=normal_sign*(/nnx(i),nny(i),nnz(i)/)
            values(7:9)=(/t1x(i),t1y(i),t1z(i)/)
            values(10:12)=normal_sign*(/t2x(i),t2y(i),t2z(i)/)
            values(13:15)=normal_sign*(/curvt1(i),curvt2(i),curvmn(i)/)/um_to_m
            offset=15

            vector=(/exE3x_EM(i),exE3y_EM(i),exE3z_EM(i)/)
            CALL AppendComplex3(values,offset,vector,1.0D0)
            vector=(/exH3x_EM(i),exH3y_EM(i),exH3z_EM(i)/)
            CALL AppendComplex3(values,offset,vector,h_scale)
            vector=(/exE3xdnn_EM(i),exE3ydnn_EM(i),exE3zdnn_EM(i)/)
            CALL AppendComplex3(values,offset,vector,dn_sign)
            vector=(/exH3xdnn_EM(i),exH3ydnn_EM(i),exH3zdnn_EM(i)/)
            CALL AppendComplex3(values,offset,vector,dn_sign*h_scale)
            vector=(/inE3x_EM(i),inE3y_EM(i),inE3z_EM(i)/)
            CALL AppendComplex3(values,offset,vector,1.0D0)
            vector=(/inH3x_EM(i),inH3y_EM(i),inH3z_EM(i)/)
            CALL AppendComplex3(values,offset,vector,h_scale)
            vector=(/inE3xdnn_EM(i),inE3ydnn_EM(i),inE3zdnn_EM(i)/)
            CALL AppendComplex3(values,offset,vector,dn_sign)
            vector=(/inH3xdnn_EM(i),inH3ydnn_EM(i),inH3zdnn_EM(i)/)
            CALL AppendComplex3(values,offset,vector,dn_sign*h_scale)
            IF (offset /= 63) CALL Die('Internal trace-column packing error.')

            WRITE(91,'(4(I0,1X),63(ES24.16E3,1X))') i,interface_id,outer_region,inner_region,values
        END DO
        CLOSE(91)
    END SUBROUTINE WriteSurfaceTraces

    ! Contract: require six-node surface elements and replace
    ! surface_elements_q6.dat with interface-labelled global connectivity.
    SUBROUTINE WriteElementsQ6
        INTEGER :: k
        IF (mxnmbrndlnkelmnt /= 6) CALL Die('Expected six nodes per Q6 element.')
        OPEN(92,FILE='surface_elements_q6.dat',STATUS='REPLACE')
        WRITE(92,'(A)') '# columns: element_id interface_id n1 n2 n3 n4 n5 n6; one-based global node IDs.'
        DO k=1,ttlnmbrelmnt
            WRITE(92,'(8(I0,1X))') k,ElementInterface(k),elmntlnknd(k,1:6)
        END DO
        CLOSE(92)
    END SUBROUTINE WriteElementsQ6

    ! Contract: replace cross_sections.dat with the 480 nm scattering,
    ! extinction and absorption cross sections supplied in um**2.
    SUBROUTINE WriteCrossSections(sca,ext,absorb)
        DOUBLE PRECISION, INTENT(IN) :: sca,ext,absorb
        OPEN(93,FILE='cross_sections.dat',STATUS='REPLACE')
        WRITE(93,'(A)') '# columns: wavelength_nm sigma_sca_um2 sigma_ext_um2 sigma_abs_um2'
        WRITE(93,'(4(ES24.16E3,1X))') 480.0D0,sca,ext,absorb
        CLOSE(93)
    END SUBROUTINE WriteCrossSections

    ! Contract: evaluate boundary-continuity and div(H) residual summaries from
    ! the solved traces, then replace surface_trace_audit_summary.txt.
    SUBROUTINE WriteAudit(mesh_level,sca,ext,absorb)
        INTEGER, INTENT(IN) :: mesh_level
        DOUBLE PRECISION, INTENT(IN) :: sca,ext,absorb
        INTEGER :: i,kk,neighbor,interface_id
        DOUBLE PRECISION :: sum_et2,sum_ht2,sum_dnD2,sum_dnB2,sum_dhex2,sum_dhin2
        DOUBLE PRECISION :: sum_hexdn2,sum_hindn2,sum_hex2,sum_hin2
        DOUBLE PRECISION :: max_et,max_ht,max_dnD,max_dnB
        COMPLEX(KIND=KIND(1.0D0)) :: eps_ex,eps_in,mu_ex,mu_in
        COMPLEX(KIND=KIND(1.0D0)) :: Eex(3),Ein(3),Hex(3),Hin(3),dHex(3),dHin(3)
        COMPLEX(KIND=KIND(1.0D0)) :: dt1Hex(3),dt2Hex(3),dt1Hin(3),dt2Hin(3)
        COMPLEX(KIND=KIND(1.0D0)) :: r1,r2,r3,r4,div_ex,div_in
        DOUBLE PRECISION :: tiny_scale

        sum_et2=0.0D0; sum_ht2=0.0D0; sum_dnD2=0.0D0; sum_dnB2=0.0D0
        sum_dhex2=0.0D0; sum_dhin2=0.0D0; sum_hexdn2=0.0D0; sum_hindn2=0.0D0
        sum_hex2=0.0D0; sum_hin2=0.0D0
        max_et=0.0D0; max_ht=0.0D0; max_dnD=0.0D0; max_dnB=0.0D0

        DO i=1,ttlnmbrnd
            interface_id=NodeInterface(i)
            IF (interface_id == 1) THEN
                eps_ex=exeps_EM; mu_ex=exmiu_EM
                eps_in=ineps_EM(1); mu_in=inmiu_EM(1)
            ELSE
                eps_ex=ineps_EM(1); mu_ex=inmiu_EM(1)
                eps_in=ineps_EM(2); mu_in=inmiu_EM(2)
            END IF
            Eex=(/exE3x_EM(i),exE3y_EM(i),exE3z_EM(i)/)
            Ein=(/inE3x_EM(i),inE3y_EM(i),inE3z_EM(i)/)
            Hex=(/exH3x_EM(i),exH3y_EM(i),exH3z_EM(i)/)
            Hin=(/inH3x_EM(i),inH3y_EM(i),inH3z_EM(i)/)
            dHex=(/exH3xdnn_EM(i),exH3ydnn_EM(i),exH3zdnn_EM(i)/)
            dHin=(/inH3xdnn_EM(i),inH3ydnn_EM(i),inH3zdnn_EM(i)/)

            r1=t1x(i)*(Eex(1)-Ein(1))+t1y(i)*(Eex(2)-Ein(2))+t1z(i)*(Eex(3)-Ein(3))
            r2=t2x(i)*(Eex(1)-Ein(1))+t2y(i)*(Eex(2)-Ein(2))+t2z(i)*(Eex(3)-Ein(3))
            r3=t1x(i)*(Hex(1)-Hin(1))+t1y(i)*(Hex(2)-Hin(2))+t1z(i)*(Hex(3)-Hin(3))
            r4=t2x(i)*(Hex(1)-Hin(1))+t2y(i)*(Hex(2)-Hin(2))+t2z(i)*(Hex(3)-Hin(3))
            sum_et2=sum_et2+ABS(r1)**2+ABS(r2)**2
            sum_ht2=sum_ht2+ABS(r3)**2+ABS(r4)**2
            max_et=MAX(max_et,ABS(r1),ABS(r2)); max_ht=MAX(max_ht,ABS(r3),ABS(r4))

            r1=eps_ex*(nnx(i)*Eex(1)+nny(i)*Eex(2)+nnz(i)*Eex(3)) &
            & -eps_in*(nnx(i)*Ein(1)+nny(i)*Ein(2)+nnz(i)*Ein(3))
            r2=mu_ex*(nnx(i)*Hex(1)+nny(i)*Hex(2)+nnz(i)*Hex(3)) &
            & -mu_in*(nnx(i)*Hin(1)+nny(i)*Hin(2)+nnz(i)*Hin(3))
            sum_dnD2=sum_dnD2+ABS(r1)**2; sum_dnB2=sum_dnB2+ABS(r2)**2
            max_dnD=MAX(max_dnD,ABS(r1)); max_dnB=MAX(max_dnB,ABS(r2))

            dt1Hex=ztpzero; dt2Hex=ztpzero; dt1Hin=ztpzero; dt2Hin=ztpzero
            DO kk=1,mxnmbrndlnknd2ndslf
                neighbor=ndlnknd2ndslf(i,kk)
                IF (neighbor /= 0) THEN
                    dt1Hex=dt1Hex+(/exH3x_EM(neighbor),exH3y_EM(neighbor),exH3z_EM(neighbor)/)*d_dt1(i,kk)
                    dt2Hex=dt2Hex+(/exH3x_EM(neighbor),exH3y_EM(neighbor),exH3z_EM(neighbor)/)*d_dt2(i,kk)
                    dt1Hin=dt1Hin+(/inH3x_EM(neighbor),inH3y_EM(neighbor),inH3z_EM(neighbor)/)*d_dt1(i,kk)
                    dt2Hin=dt2Hin+(/inH3x_EM(neighbor),inH3y_EM(neighbor),inH3z_EM(neighbor)/)*d_dt2(i,kk)
                END IF
            END DO
            div_ex=nnx(i)*dHex(1)+nny(i)*dHex(2)+nnz(i)*dHex(3) &
            & +t1x(i)*dt1Hex(1)+t1y(i)*dt1Hex(2)+t1z(i)*dt1Hex(3) &
            & +t2x(i)*dt2Hex(1)+t2y(i)*dt2Hex(2)+t2z(i)*dt2Hex(3)
            div_in=nnx(i)*dHin(1)+nny(i)*dHin(2)+nnz(i)*dHin(3) &
            & +t1x(i)*dt1Hin(1)+t1y(i)*dt1Hin(2)+t1z(i)*dt1Hin(3) &
            & +t2x(i)*dt2Hin(1)+t2y(i)*dt2Hin(2)+t2z(i)*dt2Hin(3)
            sum_dhex2=sum_dhex2+ABS(div_ex)**2; sum_dhin2=sum_dhin2+ABS(div_in)**2
            sum_hexdn2=sum_hexdn2+SUM(ABS(dHex)**2); sum_hindn2=sum_hindn2+SUM(ABS(dHin)**2)
            sum_hex2=sum_hex2+SUM(ABS(Hex)**2); sum_hin2=sum_hin2+SUM(ABS(Hin)**2)
        END DO

        tiny_scale=TINY(1.0D0)
        OPEN(94,FILE='surface_trace_audit_summary.txt',STATUS='REPLACE')
        WRITE(94,'(A)') 'Eccentric Au-core / Ag-shell two-sided surface-trace audit'
        WRITE(94,'(A,I0)') 'level = ',mesh_level
        WRITE(94,'(A,I0)') 'nodes = ',ttlnmbrnd
        WRITE(94,'(A,I0)') 'Q6 elements = ',ttlnmbrelmnt
        WRITE(94,'(A)') '2SD dH recovery = ex-side derivative is solved directly; in-side recovery contracts all six coefficients'
        WRITE(94,'(A)') 'PEC-only omitted-cross-term defect is outside this all-2SD case; no PEC patch applied'
        WRITE(94,'(A,ES24.16)') 'RMS tangential E jump = ',SQRT(sum_et2/(2.0D0*ttlnmbrnd))
        WRITE(94,'(A,ES24.16)') 'RMS tangential H jump = ',SQRT(sum_ht2/(2.0D0*ttlnmbrnd))
        WRITE(94,'(A,ES24.16)') 'RMS normal D jump = ',SQRT(sum_dnD2/ttlnmbrnd)
        WRITE(94,'(A,ES24.16)') 'RMS normal B jump = ',SQRT(sum_dnB2/ttlnmbrnd)
        WRITE(94,'(A,ES24.16)') 'max tangential E jump = ',max_et
        WRITE(94,'(A,ES24.16)') 'max tangential H jump = ',max_ht
        WRITE(94,'(A,ES24.16)') 'max normal D jump = ',max_dnD
        WRITE(94,'(A,ES24.16)') 'max normal B jump = ',max_dnB
        WRITE(94,'(A,ES24.16)') 'relative exterior divH RMS = ',SQRT(sum_dhex2/MAX(sum_hexdn2,tiny_scale))
        WRITE(94,'(A,ES24.16)') 'relative interior divH RMS = ',SQRT(sum_dhin2/MAX(sum_hindn2,tiny_scale))
        WRITE(94,'(A,3(1X,ES24.16))') 'sigma_sca/ext/abs_um2 = ',sca,ext,absorb
        WRITE(94,'(A,3(1X,ES24.16))') 'archived_L6_sigma_sca/ext/abs_um2 = ', &
        & 0.0190123091586D0,0.0367594792939D0,0.0177471701352D0
        CLOSE(94)
    END SUBROUTINE WriteAudit

END PROGRAM Main_EM_AuAg_SurfaceGradient
