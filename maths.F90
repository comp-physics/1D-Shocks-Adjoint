module maths

    contains

SUBROUTINE mv(A, v, matvec)
  USE prms
  IMPLICIT NONE
  real, INTENT(IN)  :: A(3,3,0:n), v(3,0:n)
  real, INTENT(OUT) :: matvec(3,0:n)
  integer :: i

    do i = 1,3
        matvec(i,:) = A(i,1,:)*v(1,:) + A(i,2,:)*v(2,:) + A(i,3,:)*v(3,:)
    end do

END SUBROUTINE mv

SUBROUTINE mymatmul(A,B,C)
  USE prms
  IMPLICIT NONE
  real, INTENT(IN)  :: A(3,3,0:n), B(3,3,0:n)
  real, INTENT(OUT) :: C(3,3,0:n)
  integer :: i,j

    do i = 1,3
        do j = 1,3
            C(i,j,:) = A(i,1,:)*B(1,j,:) + A(i,2,:)*B(2,j,:) + A(i,3,:)*B(3,j,:)
        end do
    end do  
END SUBROUTINE mymatmul

end module maths
