! SPDX-FileCopyrightText: 2026 Qiang Sun
! SPDX-License-Identifier: BSD-3-Clause

! Six-point Gauss-Legendre rule used by all released surface examples.
! The fixed table replaces the unused general root-finding implementation.
    SUBROUTINE GauLegCoeff1D(dmngp,dmxabsc,dmweig)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: dmngp
        DOUBLE PRECISION,INTENT(OUT) :: dmxabsc(dmngp),dmweig(dmngp)

        IF (dmngp /= 6) THEN
            WRITE(*,*) "Released quadrature supports exactly six 1-D points."
            STOP 81
        END IF
        dmxabsc = (/ -0.9324695142031521D0, -0.6612093864662645D0, &
        &             -0.2386191860831969D0,  0.2386191860831969D0, &
        &              0.6612093864662645D0,  0.9324695142031521D0 /)
        dmweig = (/ 0.1713244923791704D0, 0.3607615730481386D0, &
        &            0.4679139345726910D0, 0.4679139345726910D0, &
        &            0.3607615730481386D0, 0.1713244923791704D0 /)
    END SUBROUTINE GauLegCoeff1D

! Sixteen-point degree-eight symmetric quadrature rule on the reference
! triangle. The nodes are the rule published by Zhang, Cui and Liu,
! J. Comput. Math. 27 (2009) 89-96, Section 4.1. Their weights sum to one;
! the weights below are scaled by one half for a reference triangle of area
! one half.
    SUBROUTINE GauLegUniTriCff2D(dmngp,dmxabsc,dmyabsc,dmweig)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: dmngp
        DOUBLE PRECISION,INTENT(OUT) :: dmxabsc(dmngp),dmyabsc(dmngp),dmweig(dmngp)

        IF (dmngp /= 16) THEN
            WRITE(*,*) "Released quadrature supports exactly 16 triangle points."
            STOP 80
        END IF

        dmweig = (/ &
        &   0.07215780383889358412554555524453D0, &
        &   0.04754581713364231239694805219429D0, &
        &   0.04754581713364231239694805219429D0, &
        &   0.04754581713364231239694805219429D0, &
        &   0.05160868526735912514089577514606D0, &
        &   0.05160868526735912514089577514606D0, &
        &   0.05160868526735912514089577514606D0, &
        &   0.01622924881159904015546296417089D0, &
        &   0.01622924881159904015546296417089D0, &
        &   0.01622924881159904015546296417089D0, &
        &   0.01361515708721749713242234503695D0, &
        &   0.01361515708721749713242234503695D0, &
        &   0.01361515708721749713242234503695D0, &
        &   0.01361515708721749713242234503695D0, &
        &   0.01361515708721749713242234503695D0, &
        &   0.01361515708721749713242234503695D0 /)

        dmxabsc = (/ &
        &   0.33333333333333333333333333333333D0, &
        &   0.45929258829272315602881551449417D0, &
        &   0.08141482341455368794236897101166D0, &
        &   0.45929258829272315602881551449417D0, &
        &   0.17056930775176020662229350149146D0, &
        &   0.65886138449647958675541299701707D0, &
        &   0.17056930775176020662229350149146D0, &
        &   0.05054722831703097545842355059660D0, &
        &   0.89890554336593804908315289880680D0, &
        &   0.05054722831703097545842355059660D0, &
        &   0.00839477740995760533721383453929D0, &
        &   0.72849239295540428124100037917606D0, &
        &   0.26311282963463811342178578628464D0, &
        &   0.72849239295540428124100037917606D0, &
        &   0.26311282963463811342178578628464D0, &
        &   0.00839477740995760533721383453929D0 /)

        dmyabsc = (/ &
        &   0.33333333333333333333333333333333D0, &
        &   0.08141482341455368794236897101166D0, &
        &   0.45929258829272315602881551449417D0, &
        &   0.45929258829272315602881551449417D0, &
        &   0.65886138449647958675541299701707D0, &
        &   0.17056930775176020662229350149146D0, &
        &   0.17056930775176020662229350149146D0, &
        &   0.89890554336593804908315289880680D0, &
        &   0.05054722831703097545842355059660D0, &
        &   0.05054722831703097545842355059660D0, &
        &   0.72849239295540428124100037917606D0, &
        &   0.00839477740995760533721383453929D0, &
        &   0.72849239295540428124100037917606D0, &
        &   0.26311282963463811342178578628464D0, &
        &   0.00839477740995760533721383453929D0, &
        &   0.26311282963463811342178578628464D0 /)
    END SUBROUTINE GauLegUniTriCff2D

    subroutine fdcoef(mord,nord,x0,grid,coef)
        implicit none
        save

! Compute finite-difference weights on an arbitrarily spaced 1-D grid.
! The stable recursion follows B. Fornberg, Math. Comp. 51 (1988) 699--706.
! mord is the derivative order; nord is the stencil size; x0 is the
! evaluation point; grid contains the nodes; coef returns the weights.

! Arguments.
        integer          mord,nord
        DOUBLE PRECISION x0,grid(nord),coef(nord)


! Local variables.
        integer          nu,nn,mm,nmmin,mmax,nmax
        parameter        (mmax=8, nmax=20)
        DOUBLE PRECISION weight(mmax,nmax,nmax),c1,c2,c3,c4,pfac


! Initialise the recursion table.
        do nu=1,nord
            do nn=1,nord
                do mm=1,mord
                    weight(mm,nn,nu) = 0.0d0
                enddo
            enddo
        enddo

        weight(1,1,1) = 1.0d0
        c1            = 1.0d0
        nmmin         = min(nord,mord)
        do nn = 2,nord
            c2 = 1.0d0
            do nu=1,nn-1
                c3 = grid(nn) - grid(nu)
                c2 = c2 * c3
                c4 = 1.0d0/c3
                pfac = grid(nn) - x0
                weight(1,nn,nu) = c4 * ( pfac * weight(1,nn-1,nu) )
                do mm=2,nmmin
                    weight(mm,nn,nu) = c4 * ( pfac * weight(mm,nn-1,nu) &
                    &                 - dfloat(mm-1)*weight(mm-1,nn-1,nu) )
                enddo
            enddo
            pfac = (grid(nn-1) - x0)
            weight(1,nn,nn) = c1/c2 * (-pfac*weight(1,nn-1,nn-1))
            c4 = c1/c2
            do mm=2,nmmin
                weight(mm,nn,nn) = c4 * (dfloat(mm-1)*weight(mm-1,nn-1,nn-1) - &
                &                  pfac*weight(mm,nn-1,nn-1))
            enddo
            c1 = c2
        enddo

! Return the requested derivative weights.
        do nu = 1,nord
            coef(nu) = weight(mord,nord,nu)
        enddo
        return

    end
