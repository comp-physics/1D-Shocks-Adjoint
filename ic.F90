SUBROUTINE getic(x,uc,up,ucsol)
  USE prms
  USE conv
  USE io
  IMPLICIT NONE
  REAL, INTENT(OUT) :: uc(3,0:n), up(3,0:n), ucsol(0:nsteps,3,0:n)
  REAL :: x(0:n), udum(3,0:n)
    real :: a,b,w
  integer :: i

    FORALL(i=0:n:1) x(i) = i*dx
    x(:) = x(:) - L/2.
    
    if (which_ic.eq.1) then
        !Entropy wave: Gaussian
        do i=0,n
            up(1,i) = amp * ( exp( -(x(i)**2.) ) + 1 )
        end do
        up(2,:) = 1; up(3,:) = 1./g
    else if (which_ic.eq.2) then
        !Entropy wave: Shock
        do i=0,n
            if (x(i) > 0 ) then
                up(1,i) = amp
            else
                up(1,i) = 2.*amp
            end if
        end do
        up(2,:) = 1; up(3,:) = 1./g
    else if (which_ic.eq.3) then
        !Sod shocktube
        do i = 0,n
            if ( x(i) .le. 0 ) then
                up(1,i) = 1.;       up(2,i) = 0.; up(3,i) = 1.
            else
                up(1,i) = 0.125;    up(2,i) = 0.; up(3,i) = 0.1
            end if
        end do
    else if (which_ic.eq.4) then
         !Gaussians
        do i=0,n
            up(1,i) = amp * ( exp( -(x(i)**2.) ) + 1 ) !between amp and 2*amp
            up(2,i) =  0.1*(exp( -(x(i)**2.) )+1) !between 0.1 0.2
            up(3,i) =  0.1*exp( -(x(i)**2.) ) + 5 !between 5 and 5.1
        end do
    else if (which_ic.eq.5) then
         !Square wave
         do i=0,n
            if ( x(i) > -1  .and. x(i) < 1 ) then
                up(1,i) = 2.*amp
            else
                up(1,i) = amp
            end if
        end do
        up(2,:) = 1; up(3,:) = 1./g
    else if (which_ic.eq.6) then
         !Shocks only in density and velocity
         do i=0,n
            if ( x(i) < 0.) then
                up(1,i) = 2.*amp
            else
                up(1,i) = amp
            end if
            if ( x(i) < 1.) then
                up(2,i) = -1
            else
                up(2,i) = 1
            end if
        end do
        up(3,:) = 1./g
    else if (which_ic.eq.7) then
         !Modified Sod shock tube
         do i=0,n
            if ( x(i) < 0.) then
                up(1,i) = 1
            else
                up(1,i) = 0.5
            end if
            if ( x(i) < 1.) then
                up(2,i) = 1.
            else
                up(2,i) = 1.
            end if
            if (x(i) < -0.5 ) then 
                up(3,i) = 0.2
            else
                up(3,i) = 0.1
            end if
        end do
    else if (which_ic.eq.8) then
        !Sod shocktube2
        do i = 0,n
            if ( x(i) .le. 0 ) then
                up(1,i) = 1.;       up(2,i) = 0.1; up(3,i) = 1.
            else
                up(1,i) = 0.125;    up(2,i) = 0.1; up(3,i) = 0.1
            end if
        end do
    else if (which_ic.eq.33) then
        !Sod shocktube with hyperbolic tangent regularization
        do i = 0,n
                w = 100.
                a = 1.; b = 0.125;
                up(1,i) = (a - b)*(( Tanh(-w*x(i)) + 1)/2.) + b
                a = 1.; b = 0.1;
                up(3,i) = (a - b)*(( Tanh(-w*x(i)) + 1)/2.) + b
                
                up(2,i) = 0.
        end do
    end if

    udum = up
    do i = 1,n-1
        up(:,i) = 1./6.*(udum(:,i-1) + 4.*udum(:,i) + udum(:,i+1))
    end do

    call p_to_c(up,uc)
    call output(uc(1,:),x,'rho',0)
    ucsol(0,:,:) = uc(:,:)
END SUBROUTINE getic
