! $HeadURL$
! $Id$

!> Compute a 3rd-order WENO reconstruction.
!! Given function values \f$u(x_i)\f$ for \f$i\in\left\{0,\ldots,n\right\}\f$
!! compute the reconstruction at \f$u_{r}\left(x_{i+1/2}\right)\f$ following
!! section 3.5 of Liu, Osher, and Chan's 1994 JCP paper.
!!
!! @param ur    Reconstruction \f$u_{r}\left(x_{i+1/2}\right)\f$
!! @param u     Function values \f$u(x_i)\f$
!! @param n     Grid size
!! @param bias  If strictly positive, bias stencil to the left.
!!              Otherwise, bias stencil to the right.
SUBROUTINE reconstruct3 (ur, u, n, bias)
! Equation numbers in the implementation refer to Liu, Osher, and Chan's paper.

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n, bias
  REAL, INTENT(IN) :: u(0:n)
  REAL, INTENT(OUT) :: ur(0:n)
  REAL, PARAMETER :: eps = 1.d-6  ! guarantee nonzero denominator
  REAL :: beta(1:2)               ! smoothness indicators
  REAL :: w(1:2)                  ! nonlinear weights
  REAL :: wt(1:2), wtsumi         ! temporary nonlinear weights
  REAL :: urloc(1:2)              ! the two local reconstructions
  REAL :: a(1:2,1:2)              ! weights in reconstruction
  INTEGER :: i
  REAL :: v(-1:n+2)               ! add on periodic BCs
  REAL :: v0, vp, vm              ! local values

  a(1,1) = -1.d0 / 2.d0
  a(1,2) =  3.d0 / 2.d0
  a(2,1) =  1.d0 / 2.d0
  a(2,2) =  1.d0 / 2.d0

! Add on periodic boundary conditions
! this is wasteful but results in a single loop so the code is easier to read
  v(0:n) = u(0:n)
  v(-1)  = u(n-1)
  v(n+1:n+2) = u(1:2)

 IF (bias > 0) THEN ! Bias to the left, case 1 in section 3.5
   DO i = 0, n, 1
     v0 = v(i)
     vp = v(i+1)
     vm = v(i-1)
! The reconstructed values at x(i+1/2) per p'(j), p'(j+1) from bottom of p205
! Note mistake in the p'j formula, i.e. (x-x).
     urloc(1) = a(1,1) * vm + a(1,2) * v0
     urloc(2) = a(2,1) * v0 + a(2,2) * vp
! Smoothness indicators from p206 just above equation 3.16
     beta(1) = (v0 - vm)**2
     beta(2) = (vp - v0)**2
! Compute nonlinear weights (3.17a)
     wt(1) = 0.5d0 / ((eps + beta(1))**2)
     wt(2) = 1.0d0 / ((eps + beta(2))**2)
     wtsumi = 1.d0 / (wt(1) + wt(2))
     w(1) = wt(1) * wtsumi
     w(2) = wt(2) * wtsumi
! Finally reconstruct, formula (3.16)
     ur(i) = w(1) * urloc(1) + w(2) * urloc(2)
   END DO
  ELSE ! biased to the right, case 2 in section 3.5
    DO i = 1, n+1, 1
      v0 = v(i)
      vp = v(i+1)
      vm = v(i-1)
! The reconstructed values at x(i-1/2) per p'(j), p'(j+1) from bottom of p205
! Note mistake in the p'j formula, i.e. (x-x).
      urloc(1) = a(2,1) * vm + a(2,2) * v0
      urloc(2) = a(1,2) * v0 + a(1,1) * vp
! Smoothness indicators from p206 just above equation 3.16
      beta(1) = (v0 - vm)**2
      beta(2) = (vp - v0)**2
! Compute nonlinear weights (3.17a)
      wt(1) = 1.0d0 / ((eps + beta(1))**2)
      wt(2) = 0.5d0 / ((eps + beta(2))**2)
      wtsumi = 1.d0 / (wt(1) + wt(2))
      w(1) = wt(1) * wtsumi
      w(2) = wt(2) * wtsumi
! Finally reconstruct, formula (3.16)
      ur(i-1) = w(1) * urloc(1) + w(2) * urloc(2)
    END DO
  END IF
END SUBROUTINE reconstruct3
