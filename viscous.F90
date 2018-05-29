SUBROUTINE viscosity(u,fin,flag)
  use prms
  IMPLICIT NONE
  REAL, intent(IN)   :: u(3,0:n)
  REAL, INTENT(INOUT)  :: fin(3,0:n)
  REAL :: uxx(3,0:n)
  REAL               :: c(0:1)
  INTEGER                       :: i
  character(1), INTENT(IN)      :: flag

    c(0) = eps * (-2.)*dxi**2
    c(1) = eps * (+1.)*dxi**2

    uxx = 0.
    do i = 5,n-6
        uxx(:,i) = c(1)*u(:,i-1) + c(0)*u(:,i) + c(1)*u(:,i+1)
    end do

    if (flag.eq.'F') then       !forward viscosity 
        fin = fin + uxx
    else if (flag.eq.'A') then  !adjoint viscosity
        fin = fin + uxx
    end if
END SUBROUTINE viscosity
