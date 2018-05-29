module matrices

contains
    
SUBROUTINE getA(At, uc)
  !provides the transpose of the 1d inviscid flux jacobian
  use prms
  IMPLICIT NONE
  real, INTENT(IN)  :: uc(3,0:n)
  real, INTENT(OUT) :: At(3,3,0:n)
  integer :: i

    At(1,1,:) = 0.;
    At(2,1,:) = 1.;
    At(3,1,:) = 0.;

    At(1,2,:) = ( (g-3)/2. )*( (uc(2,:)/uc(1,:))**2. )
    At(2,2,:) = (3-g)*uc(2,:)/uc(1,:);
    At(3,2,:) = g-1;

    At(1,3,:) = -1.*g*uc(2,:)*uc(3,:)/(uc(1,:)**2.) + (g-1)*( (uc(2,:)/uc(1,:))**3. );
    At(2,3,:) = g*uc(3,:)/uc(1,:) - 1.5*(g-1)*( (uc(2,:)/uc(1,:))**2. );
    At(3,3,:) = g*uc(2,:)/uc(1,:)
END SUBROUTINE getA

SUBROUTINE getL(Lambda, uc)
  !provides wavespeeds
  use prms
  use conv
  IMPLICIT NONE
  real, INTENT(IN)  :: uc(3,0:n)
  real, INTENT(OUT) :: Lambda(3,0:n)
  real, dimension(0:n) :: c
  real, dimension(3,0:n) :: up
  integer :: i

    call c_to_p(up,uc)

    c(:) = sqrt( g*up(3,:)/up(1,:) )
    Lambda(1,:) = up(2,:)
    Lambda(2,:) = up(2,:)+c(:)
    Lambda(3,:) = up(2,:)-c(:)

END SUBROUTINE getL

SUBROUTINE getB(B, uc, drdt, drdx)
  !B = (Q_t + Lam Q_x )Qinv
  use prms
  IMPLICIT NONE
  real, INTENT(IN)  :: uc(3,0:n), drdt(0:n), drdx(0:n)
  real, INTENT(OUT) :: B(3,3,0:n)
  real, dimension(0:n) :: aa, dd !drdx

    aa(:) = 1./uc(1,:) * ( drdt(:) + drdx(:)*( 1./sqrt(uc(1,:)) + vinf) )
    dd(:) = 1./uc(1,:) * ( drdt(:) + drdx(:)*(-1./sqrt(uc(1,:)) + vinf) )
  
    B(1,1,:) = 0.; B(1,2,:) = 0.; B(1,3,:) = 0.;
    B(2,1,:) = aa(:); B(2,2,:) = -0.75*aa(:); B(2,3,:) = -0.25*aa(:);
    B(3,1,:) = dd(:); B(3,2,:) = -0.25*dd(:); B(3,3,:) = -0.75*dd(:);
END SUBROUTINE getB

SUBROUTINE getBfull(B, uc, dpdt, dpdx)
  !B = (Q_t + Lam Q_x )Qinv
  use prms
  use conv
  use maths
  IMPLICIT NONE
  real, INTENT(IN)  :: uc(3,0:n), dpdt(3,0:n), dpdx(3,0:n)
  real, INTENT(OUT) :: B(3,3,0:n)
  real, dimension(3,3,0:n)  :: dQt, dQx, Qinv, Lambda
  real, dimension(3,3,0:n)  :: LQx, QtLQx
  real, dimension(3,0:n)    :: wavespeeds
    
    call getdQ(dQt,dpdt,uc)
    call getdQ(dQx,dpdx,uc)
    call getL(wavespeeds,uc)
    call getQinv(Qinv,uc)

    Lambda = 0.;
    Lambda(1,1,:) = wavespeeds(1,:)
    Lambda(2,2,:) = wavespeeds(2,:)
    Lambda(3,3,:) = wavespeeds(3,:)

    call mymatmul(Lambda,dQx,LQx)
    QtLQx = dQt + LQx
    call mymatmul(QtLQx,Qinv,B)
    B = -1.*B
END SUBROUTINE getBfull

SUBROUTINE getdAT(Bx, uc, dpdx)
  !(A^T)_x
  use prms
  use conv
  use maths
  IMPLICIT NONE
  real, INTENT(IN)  :: uc(3,0:n), dpdx(3,0:n)
  real, INTENT(OUT) :: Bx(3,3,0:n)
  real, dimension(3,0:n) :: up 
  real :: kappa
    kappa = g/(g-1.)
    call c_to_p(up,uc)

    Bx = 0.
    Bx(1,2,:) = (g-3)*up(2,:)*dpdx(2,:) 
    Bx(1,3,:) = kappa*up(3,:)*up(2,:)*dpdx(1,:)/(up(1,:)*up(1,:)) - 0.5*up(2,:)*(4.-4*g+3*g*up(2,:))*dpdx(2,:)+ &
            & (g*up(2,:)*dpdx(3,:)+g*up(3,:)*dpdx(2,:))/(up(1,:)-g*up(1,:))


    Bx(2,2,:) = (3-g)*dpdx(2,:)
    Bx(2,3,:) = kappa*(up(1,:)*dpdx(3,:) - up(3,:)*dpdx(1,:))/(up(1,:)*up(1,:)) + 0.5*up(2,:)*(2*g-9*(g-1)*up(2,:))*dpdx(2,:) 

    Bx(3,3,:) = g*dpdx(2,:)
END SUBROUTINE getdAT

SUBROUTINE getQ(Q, uc)
  !Q = transpose( P ); P = evecs(A)
  use prms
  use conv
  IMPLICIT NONE
  real, INTENT(IN)  :: uc(3,0:n)
  real, INTENT(OUT) :: Q(3,3,0:n)
  real, dimension(0:n) :: c, H, v2
  real, dimension(3,0:n) :: up 

    call c_to_p(up,uc)

    v2(:) = up(2,:)*up(2,:)

    c(:) = sqrt( g*up(3,:)/up(1,:) )
    H(:) = g/(g-1.) * up(3,:)/up(1,:) + 0.5*v2(:)
  
    Q(1,1,:) = 1.; Q(1,2,:) = up(2,:);         Q(1,3,:) = 0.5*v2(:)
    Q(2,1,:) = 1.; Q(2,2,:) = up(2,:) + c(:);  Q(2,3,:) = H(:) + up(2,:)*c(:)
    Q(3,1,:) = 1.; Q(3,2,:) = up(2,:) - c(:);  Q(3,3,:) = H(:) - up(2,:)*c(:)
END SUBROUTINE getQ

SUBROUTINE getQinv(Q, uc)
  !Qinv = inverse(Q)
  use prms
  use conv
  IMPLICIT NONE
  real, INTENT(IN)  :: uc(3,0:n)
  real, INTENT(OUT) :: Q(3,3,0:n)
  real, dimension(0:n) :: c, H, v2
  real, dimension(3,0:n) :: up

    call c_to_p(up,uc)
    v2(:) = up(2,:)*up(2,:)

    c(:) = sqrt( g*up(3,:)/uc(1,:) )
    H(:) = g/(g-1.) * up(3,:)/uc(1,:) + 0.5*v2(:)
  
    Q(1,1,:) = 2.+ (2.*H(:))/(v2(:)-2.*H(:))
    Q(2,1,:) = -2.*up(2,:)/(v2(:)-2.*H(:));
    Q(3,1,:) = 2./(v2(:)-2.*H(:));

    Q(1,2,:) = up(2,:)/c(:) * ( up(2,:)*(c(:)+up(2,:)) - 2.*H(:) )/(4.*H(:) - 2.*v2(:))
    Q(2,2,:) = 0.5/c(:) + up(2,:)/(v2(:)-2.*H(:))
    Q(3,2,:) = 1./(2.*H(:) - v2(:))

    Q(1,3,:) =  0.5*up(2,:)*(1/c(:) + up(2,:)/(2.*H(:)-v2(:)))
    Q(2,3,:) = -0.5/c(:) + up(2,:)/(v2(:)-2.*H(:))
    Q(3,3,:) = 1./(2.*H(:) - v2(:))
END SUBROUTINE getQinv

SUBROUTINE getdQ(dQ, dup, uc)
  !dQ/dx or dQ/dt
  use prms
  use conv
  IMPLICIT NONE
  real, INTENT(IN)  :: uc(3,0:n), dup(3,0:n)
  real, INTENT(OUT) :: dQ(3,3,0:n)
  real, dimension(0:n) :: c, cx, vcx, Hx
  real, dimension(3,0:n) :: up
  real :: kappa 

    call c_to_p(up,uc)
    
    kappa  = g/(g-1)
    c(:)   = sqrt( g*up(3,:)/uc(1,:) )
    Hx(:)  = kappa*(dup(3,:)/up(1,:) - up(3,:)*dup(1,:)/(up(1,:)*up(1,:)) ) + up(2,:)*dup(2,:)
             !g/(g-1.) * (up(1,:)*dup(3,:) - up(3,:)*dup(1,:)) /&
             ! & ( up(1,:)*up(1,:) ) + up(2,:)*dup(2,:)
    cx(:)  = 0.5*g*(  up(1,:)*dup(3,:) - up(3,:)*dup(1,:) ) /&
             ( sqrt(g*up(3,:)/up(1,:)) * up(1,:)*up(1,:) )  
             !sqrt(g/(up(1,:)*up(3,:))) * &
             !& ( up(1,:)*dup(3,:) - up(3,:)*dup(1,:) ) /&
             !& (2.*up(1,:))
    vcx(:) = sqrt(g/(up(1,:)*up(3,:))) * &
              & ( up(1,:)*(up(2,:)*dup(3,:) + 2*up(3,:)*dup(2,:)) - up(2,:)*up(3,:)*dup(1,:) ) /&
              & (2.*up(1,:))
            !up(2,:)*cx(:) + dup(2,:)*c(:)
    
    dQ(1,1,:) = 0.; dQ(1,2,:) = dup(2,:);          dQ(1,3,:) = up(2,:)*dup(2,:)
    dQ(2,1,:) = 0.; dQ(2,2,:) = dup(2,:) + cx(:);  dQ(2,3,:) = Hx(:) + vcx(:)
    dQ(3,1,:) = 0.; dQ(3,2,:) = dup(2,:) - cx(:);  dQ(3,3,:) = Hx(:) - vcx(:)

!dQ(2,3,:) = dup(2,:)*up(2,:) - ( (dup(1,:)*up(3,:) - dup(3,:)*up(1,:))*(2*up(3,:)*g + up(1,:)*up(2,:)*(g-1)*c(:) ) ) /&
!    (2*up(3,:)*up(1,:)*up(1,:)*(g-1))

!dQ(3,3,:) = dup(2,:)*up(2,:) - ( (dup(1,:)*up(3,:) - dup(3,:)*up(1,:))*(2*up(3,:)*g - up(1,:)*up(2,:)*(g-1)*c(:) ) ) /&
!    (2*up(3,:)*up(1,:)*up(1,:)*(g-1))

END SUBROUTINE getdQ

subroutine p_to_q(p,ucfor)
    use prms
    !use matrices
    use maths
    implicit none
    real, intent(inout) :: p(3,0:n)
    real, intent(in) :: ucfor(3,0:n)
    real :: qdum(3,0:n)
    real :: QQ(3,3,0:n)

    call getQ(QQ,ucfor)
    call mv(QQ,p,qdum)
    p=qdum
end subroutine p_to_q

subroutine q_to_p(p,ucfor)
    use prms
    !use matrices
    use maths
    implicit none
    real, intent(inout) :: p(3,0:n)
    real, intent(in) :: ucfor(3,0:n)
    real :: qdum(3,0:n)
    real :: QQ(3,3,0:n)

    call getQinv(QQ,ucfor)
    call mv(QQ,p,qdum)
    p=qdum
end subroutine q_to_p

end module matrices
