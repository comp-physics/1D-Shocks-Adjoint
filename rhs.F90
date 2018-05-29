module rhs

    contains

SUBROUTINE rhs_WENO(fout, fin)
  USE prms
  IMPLICIT NONE
  integer :: i
  REAL, INTENT(OUT) :: fout(3,0:n)
  REAL, INTENT(IN) :: fin(3,0:n)

  fout = 0.
  do i = 1,3
    fout(i,6:n-5) = -dxi * (fin(i,6:n-5) - fin(i,5:n-6))
    fout(i,0:5) = 0.; fout(i,n-4:n) = 0.
  end do
END SUBROUTINE rhs_WENO

SUBROUTINE rhs_FV(fout, fin, uin)
  use prms
  IMPLICIT NONE
  integer :: i
  real, INTENT(IN)  :: fin(3,0:n), uin(3,0:n)
  real, INTENT(OUT) :: fout(3,0:n)
  real :: myeps

    myeps = maxwavespeed

    fout = 0.
    do i = 1,3
        fout(i,1:n-1) = -dxi * 0.5 *( fin(i,2:n) - fin(i,0:n-2) - myeps*(uin(i,2:n) -2.*uin(i,1:n-1) + uin(i,0:n-2)) )
    end do
    fout(i,0:5) = 0.; fout(i,n-4:n) = 0.
END SUBROUTINE rhs_FV

SUBROUTINE rhs_FVb(fout, uin)
  use prms
  IMPLICIT NONE
  integer :: i
  real, INTENT(IN)  :: uin(3,0:n)
  real, INTENT(OUT) :: fout(3,0:n)
  real :: myeps

    !myeps = maxwavespeed
    !myeps = 1
    !myeps = minval( (/ maxwavespeed, 1. /) )
    myeps = 0.

    fout = 0.
    do i = 1,3
        !fout(i,1:n-1) = dxi * 0.5 *( uin(i,2:n) - uin(i,0:n-2) + myeps*(uin(i,2:n) -2.*uin(i,1:n-1) + uin(i,0:n-2)) )
        fout(i,5:n-6) = dxi * 0.5 *( uin(i,6:n-5) - uin(i,4:n-7) + myeps*(uin(i,6:n-5) -2.*uin(i,5:n-6) + uin(i,4:n-7)) )
    end do
    fout(i,0:5) = 0.; fout(i,n-4:n) = 0.
END SUBROUTINE rhs_FVb

SUBROUTINE rhs_split(fout, fin, uin)
  use prms
  IMPLICIT NONE
  integer :: i
  real, INTENT(IN)  :: uin(3,0:n), fin(3,0:n)
  real, INTENT(OUT) :: fout(3,0:n)
  real :: myeps

    !myeps = maxwavespeed
    myeps = 0.
    
    fout = 0.
    do i = 1,3
        !fout(i,1:n-1) = dxi * 0.5 *( fin(i,2:n) - fin(i,0:n-2) + myeps*(uin(i,2:n) -2.*uin(i,1:n-1) + uin(i,0:n-2)) )
        fout(i,5:n-6) = dxi * 0.5 *( uin(i,6:n-5) - uin(i,4:n-7) + myeps*(uin(i,6:n-5) -2.*uin(i,5:n-6) + uin(i,4:n-7)) )
    end do
    fout(i,0:5) = 0.; fout(i,n-4:n) = 0.
END SUBROUTINE rhs_split

SUBROUTINE rhs_diag(fout,Lambdas,uin)
  use prms
  IMPLICIT NONE
  integer :: i,j
  real, INTENT(IN) :: Lambdas(3,0:n), uin(3,0:n)
  real, INTENT(OUT) :: fout(3,0:n)
  real :: myeps

    !myeps = maxwavespeed
    !myeps = 1.
    myeps = 0.

    fout = 0.
    do i = 1,3
        fout(i,1:n-1) = dxi * 0.5 *( lambdas(i,1:n-1)*(uin(i,2:n) - uin(i,0:n-2)) + &
                        myeps*(uin(i,2:n) -2.*uin(i,1:n-1) + uin(i,0:n-2)) )
    end do
    fout(:,0:1) = 0.; fout(:,n-1:n) = 0.
END SUBROUTINE rhs_diag

end module rhs
