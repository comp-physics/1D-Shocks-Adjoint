PROGRAM main
    USE prms
    USE conv
    USE rhs
    USE maths
    USE matrices
    USE fluxes
    USE io

    IMPLICIT NONE

    real                                :: J0
    real, DIMENSION(:),     ALLOCATABLE :: x, Js, Fd, perts, errs
    real, DIMENSION(:,:),   ALLOCATABLE :: up, uc, frp, frm, fp, fm, f, ff, pctemp, pc, matvec, Bp, Bxp, dJdPhi
    real, dimension(:,:),   ALLOCATABLE :: Lambda, dpdt, dcdt, dpdx, up0, src
    real, DIMENSION(:,:,:), ALLOCATABLE :: urk, prk, ucsol, A, B, Bx
    integer :: i
    character(1) :: flag

#if WENOORDER == 3
#define RECONSTRUCT_FUNCTION reconstruct3
#elif WENOORDER == 5
#define RECONSTRUCT_FUNCTION reconstruct5
#else
  #error "WENOORDER not #defined or unknown"
#endif
    print '(" WENO reconstruction order = ", I2)', WENOORDER

    call getprms

    allocate( x(0:n) )
    allocate( up(3,0:n), uc(3,0:n), frp(3,0:n), frm(3,0:n), f(3,0:n), fp(3,0:n), &
            fm(3,0:n), ff(3,0:n), pc(3,0:n), pctemp(3,0:n), matvec(3,0:n) )
    allocate( urk(3,0:n,1:2), prk(3,0:n,1:2), A(3,3,0:n), ucsol(0:nsteps,3,0:n), up0(3,0:n) )
    allocate( dpdt(3,0:n), dcdt(3,0:n), dpdx(3,0:n), B(3,3,0:n), Lambda(3,0:n) )
    
    allocate( Bp(3,0:n), Bxp(3,0:n), Bx(3,3,0:n) )

    allocate( Js(npert), Fd(npert), perts(npert), src(3,0:n), dJdPhi(3,0:n), errs(npert) )
    
    call getic(x,uc,up,ucsol)
    
    print*, 'Forward simulation'; flag='F'
    do nt = 1,nsteps,1
        if (mod(nt,outfile).eq.0) print*, 'Forward: Timestep = ', nt

        if (FT.eq.'Euler') then
            call c_to_p(up,uc)
            call euler_fluxes(ff,up)
            if (SI.eq.'WENO') then
                call numflux(fp, fm, uc, ff)
                call recon(fp,fm,f)
                call rhs_WENO(fp, f)
            else if (SI.eq.'FV') then
                call rhs_FV(fp,ff,uc)
            end if    
            if (enable_f_viscosity.eq.'True') call viscosity(uc,fp,flag)
            uc = uc + dt * fp
        else if (FT.eq.'RK3') then
            !substep 1
            call c_to_p(up,uc)
            call euler_fluxes(ff,up)
            if (SI.eq.'WENO') then
                call numflux(fp, fm, uc, ff)
                call recon(fp,fm,f)
                call rhs_WENO(fp, f)
            else if (SI.eq.'FV') then
                call rhs_FV(fp,ff,uc)
            end if    
            if (enable_f_viscosity.eq.'True') call viscosity(uc,fp,flag)
            urk(:,:,1)= uc(:,:) + dt * fp(:,:)

            !substep 2
            call c_to_p(up,urk(:,:,1))
            call euler_fluxes(ff,up)
            if (SI.eq.'WENO') then
                call numflux(fp, fm, urk(:,:,1), ff)
                call recon(fp,fm,f)
                call rhs_WENO(fp, f)
            else if (SI.eq.'FV') then
                call rhs_FV(fp,ff,urk(:,:,1))
            end if    
            if (enable_f_viscosity.eq.'True') call viscosity(urk(:,:,1),fp,flag)
            urk(:,:,2) = rk(1,1) * uc(:,:) + rk(1,2) * (urk(:,:,1) + dt * fp(:,:))

            !substep 3
            call c_to_p(up,urk(:,:,2))
            call euler_fluxes(ff,up)
            if (SI.eq.'WENO') then
                call numflux(fp, fm, urk(:,:,2), ff)
                call recon(fp,fm,f)
                call rhs_WENO(fp, f)
            else if (SI.eq.'FV') then
                call rhs_FV(fp,ff,urk(:,:,2))
            end if    
            if (enable_f_viscosity.eq.'True') call viscosity(urk(:,:,2),fp,flag)
            uc(:,:) = rk(2,1) * uc(:,:) + rk(2,2) * (urk(:,:,2) + dt * fp(:,:))
        end if

        ucsol(nt,:,:) = uc(:,:)

        call output(uc(1,:),x,'rho',nt)
        call output(uc(2,:),x,'mom',nt)
        call output(uc(3,:),x,'E',nt)
        if (maxval(abs(uc)).gt.100 ) stop 'Crash Forward'
    end do
    J0 = 0.5*dx*sum( uc(:,:)**2. )

    call output(uc(1,:),x,'rho',-1)
    call output(uc(2,:),x,'mom',-1)
    call output(uc(3,:),x,'E',-1)
    !call pas(uc(1,:))

    !adjoint solution
    if (enable_adjoint.eq.'False') stop

    pc = uc; flag='A'
    dJdPhi = dJdPhi+dt*pc
    if (enable_diag.eq.'True') call p_to_q(pc,ucsol(nsteps,:,:))
    do nt = nsteps-2,0,-1
        if (mod(nt,outfile).eq.0) print*, 'Adjoint: Timestep =', nt
        
        if (enable_diag .eq. 'False') then
            call getA(A,ucsol(nt,:,:))
        else if (enable_diag .eq. 'True') then
            call c_to_p(up,ucsol(nt+1,:,:))
            call c_to_p(up0,ucsol(nt,:,:))
            dcdt(:,:) = (ucsol(nt+1,:,:) - ucsol(nt,:,:))/dt
            !call getpt(dpdt,dcdt,up)
            call getpx(dpdx,-1.*dcdt,up)
            dpdt(:,:) = (up(:,:) - up0(:,:))/dt
            call getBfull(B,ucsol(nt,:,:),dpdt,dpdx)
            call getL(Lambda,ucsol(nt,:,:))       
        else if (enable_diag .eq. 'Split') then
            call getA(A,ucsol(nt+1,:,:))
            call c_to_p(up,ucsol(nt+1,:,:))
            dcdt = 0
            call getpx(dpdx,-1.*dcdt,up)
            call getdAT(Bx,ucsol(nt+1,:,:),dpdx)
        end if

        if (BT.eq.'Euler') then
            if (enable_diag.eq.'False') then
                if (SI.eq.'WENO') then
                    call adjointflux(fp, fm, pc)
                    call recon(fp,fm,f)
                    call rhs_WENO(fp, f)
                else if (SI.eq.'FV') then
                    call rhs_FVb(fp,pc)
                end if
                call mv(A,fp,matvec)
                if (enable_viscosity.eq.'True') call viscosity(pc,matvec,flag)
                pc = pc + dt * matvec
            else if (enable_diag.eq.'True') then
                call rhs_diag(fp, Lambda, pc)
                !if (enable_viscosity.eq.'True') call viscosity(pc,fp,flag)
                pc = pc + dt*fp
                call mv(B,pc,matvec)
                pc = pc + dt*matvec
                
                if (enable_viscosity.eq.'True') then
                    !apply viscosity to p, not q
                    fp = 0; pctemp = pc
                    call q_to_p(pctemp,ucsol(nt,:,:))
                    call viscosity(pctemp,fp,flag)
                    pctemp = pctemp + dt*fp
                    call p_to_q(pctemp,ucsol(nt,:,:))
                    pc = pctemp
                end if
            else if (enable_diag.eq.'Split') then
                call mv(A,pc,Bp)
                call mv(Bx,pc,Bxp)
                call rhs_split(fp,Bp,Bp)
                if (enable_viscosity.eq.'True') call viscosity(pc,fp,flag)
                pc = pc + dt*(fp-Bxp)
            end if
        else if (BT.eq.'RK3') then
            !substep 1
            if (enable_diag.eq.'False') then
                if (SI.eq.'WENO') then
                    call adjointflux(fp, fm, pc)
                    call recon(fp,fm,f)
                    call rhs_WENO(fp, f)
                else if (SI.eq.'FV') then
                    call rhs_FVb(fp,pc)
                end if
                call mv(A,fp,matvec)
                if (enable_viscosity.eq.'True') call viscosity(pc,matvec,flag)
                prk(:,:,1)  = pc(:,:) + dt * matvec(:,:)
            else if (enable_diag.eq.'True') then
                call rhs_diag(fp, Lambda, pc)
                if (enable_viscosity.eq.'True') call viscosity(pc,fp,flag)
                prk(:,:,1)  = pc(:,:) + dt * fp(:,:)
                call mv(B,prk(:,:,1),matvec)
                prk(:,:,1)  = pc(:,:) + dt * matvec(:,:)
            end if

            !substep 2
            if (enable_diag.eq.'False') then
                if (SI.eq.'WENO') then
                    call adjointflux(fp, fm, prk(:,:,1))
                    call recon(fp,fm,f)
                    call rhs_WENO(fp, f)
                else if (SI.eq.'FV') then
                    call rhs_FVb(fp,pc)
                end if
                call mv(A,fp,matvec)
                if (enable_viscosity.eq.'True') call viscosity(prk(:,:,1),matvec,flag)
                prk(:,:,2)  = rk(1,1) * pc(:,:) + rk(1,2) * (prk(:,:,1) + dt * matvec(:,:))
            else if (enable_diag.eq.'True') then
                call rhs_diag(fp, Lambda, prk(:,:,1))
                if (enable_viscosity.eq.'True') call viscosity(prk(:,:,1),fp,flag)
                prk(:,:,2)  = rk(1,1) * pc(:,:) + rk(1,2) * (prk(:,:,1) + dt * fp(:,:))
                call mv(B,prk(:,:,2),matvec)
                prk(:,:,2)  = rk(1,1) * pc(:,:) + rk(1,2) * (prk(:,:,1) + dt * matvec(:,:)) 
            end if

            !substep 3
            if (enable_diag.eq.'False') then
                if (SI.eq.'WENO') then
                    call adjointflux(fp, fm, prk(:,:,2))
                    call recon(fp,fm,f)
                    call rhs_WENO(fp, f)
                else if (SI.eq.'FV') then
                    call rhs_FVb(fp,pc)
                end if
                call mv(A,fp,matvec)
                if (enable_viscosity.eq.'True') call viscosity(prk(:,:,2),matvec,flag)
                pc(:,:)     = rk(2,1) * pc(:,:) + rk(2,2) * (prk(:,:,2) + dt * matvec(:,:)) 
            else if (enable_diag.eq.'True') then
                call rhs_diag(fp, Lambda, prk(:,:,2))
                if (enable_viscosity.eq.'True') call viscosity(prk(:,:,2),fp,flag)
                pc(:,:)     = rk(2,1) * pc(:,:) + rk(2,2) * (prk(:,:,2) + dt * fp(:,:))
                call mv(B,prk(:,:,2),matvec)
                pc(:,:)     = rk(2,1) * pc(:,:) + rk(2,2) * (prk(:,:,2) + dt * matvec(:,:))             
            end if
        end if
        
        dJdPhi = dJdPhi + dt*pc
        if (enable_diag.eq.'True') then 
            call output_diag(pc,ucsol(nt,:,:),x,'adj_',nt)
        else
            call output(pc(1,:),x,'adjrho',nt); call output(pc(2,:),x,'adjvel',nt); call output(pc(3,:),x,'adjE',nt)
        end if
        if (maxval(abs(pc)).gt.10 ) stop 'Crash Adj'
    end do
    if (enable_diag.eq.'True') call output_diag(pc,ucsol(0,:,:),x,'adj_',nt)
    if (enable_diag.eq.'True') call q_to_p(pc,ucsol(0,:,:))
    call output(pc(1,:),x,'adjrho',-1); call output(pc(2,:),x,'adjvel',-1); call output(pc(3,:),x,'adjE'  ,-1)
    !call out_targ(pc(1,:),x,'plateauL.',-0.75); call out_targ(pc(1,:),x,'plateauR.',1.0)

    print*, 'Check sensitivity'; flag='F'
    if (enable_acc.eq.'True') then

    call logspace(1.E-8,1.E-1,perts)
    do i = 1,npert
        x = 0.; uc = 0.; up = 0.; ucsol = 0.
        call getic(x,uc,up,ucsol)
        src = -1.*perts(i)*dJdPhi
        do nt = 1,nsteps,1
            if (FT.eq.'Euler') then
                call c_to_p(up,uc)
                call euler_fluxes(ff,up)
                if (SI.eq.'WENO') then
                    call numflux(fp, fm, uc, ff)
                    call recon(fp,fm,f)
                    call rhs_WENO(fp, f)
                else if (SI.eq.'FV') then
                    call rhs_FV(fp,ff,uc)
                end if    
                if (enable_f_viscosity.eq.'True') call viscosity(uc,fp,flag)
                uc = uc + dt * (fp+src)
            else if (FT.eq.'RK3') then
                stop 'rk3 not yet supported'
            end if
            !if (maxval(abs(uc)).gt.100 ) then
            !    print*, 'Crash', perts(i)
            !    exit
            !end if
        end do
        Js(i) = 0.5*dx*sum( uc(:,:)**2. )
        Fd(i) = (Js(i)-J0)/perts(i)
        errs(i) = Abs( Abs( Fd(i) ) - Abs( dx*sum(dJdPhi(:,:)**2.)) )  
    end do

    open(1,file='O/acc.dat')
    do i = 1,npert
        write(1,*) perts(i), errs(i)
    end do        

    open(1,file='O/dJdPhi1.dat')
    do i = 0,n
        write(1,*) x(i), djdphi(1,i)
    end do        
    close(1)

    end if
 
    stop
END PROGRAM main

!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine recon(fp,fm,f)
  use prms
  IMPLICIT NONE
  real, INTENT(IN) :: fp(3,0:n), fm(3,0:n)
  real :: frp(3,0:n), frm(3,0:n)
  real, INTENT(OUT) :: f(3,0:n)
  integer :: i

    do i = 1,3
        call RECONSTRUCT_FUNCTION (frp(i,:), fp(i,:), n,  1)
        call RECONSTRUCT_FUNCTION (frm(i,:), fm(i,:), n, -1)
        f(i,:) = frp(i,:) + frm(i,:)
    end do
end subroutine recon

SUBROUTINE moveback(u)
  USE prms
  IMPLICIT NONE
  REAL, INTENT(INOUT)    :: u(3,0:n)
  REAL                   :: udum(3,0:n)
  INTEGER :: i, loc

    udum = u
    !find middle shock point
    loc = minloc( abs(u(1,:)-1.5*amp), dim=1 )
    print*, 'loc = ', loc, u(:,loc)
    loc = loc-n/2

    do i = 1,3
        u(i,:) = udum(i,n)
    end do

    do i = 0,n-loc
        u(:,i) = udum(:,i+loc) !shift 'u' leftward by 'loc' gridpoints
    end do
END SUBROUTINE moveback

subroutine getpx(dpdx,ff,up)
    use prms
    implicit none
    real, intent(out) :: dpdx(3,0:n)
    real, intent(in) :: ff(3,0:n), up(3,0:n)
    real, dimension(0:n) :: rv2, num, den
    real :: kappa
    integer :: i

    rv2(:) = up(1,:)*up(2,:)*up(2,:)
    kappa = g/(g-1.)

    ! drho/dx

    !num(:) = 2*ff(1,:)*kappa*up(3,:) - 2.*ff(3,:)*up(1,:) + 2.*ff(2,:)*kappa*up(1,:)*up(2,:) + ff(1,:)*(3-4*kappa)*rv2(:)
    !den(:) = 2*kappa*up(3,:)*up(2,:) - 2.*(kappa-1)*up(2,:)*rv2(:)

    !do i = 0,n
    !    if (abs(den(i)) > 0.000001) then
    !        dpdx(1,i) = num(i)/den(i)
    !    else
    !        dpdx(1,i) = 0
    !    end if
    !end do

    do i = 1,n-1
        dpdx(:,i) = (1/(2.*dx))*(up(:,i+1)-up(:,i-1))
        !dpdx(:,i) = (1/(12.*dx))*(-1*up(:,i+2)+8*up(:,i+1)-8*up(:,i-1)+up(:,i-2))
    end do

    !dpdx(1,:) = ( 2*ff(1,:)*kappa*up(3,:) - 2.*ff(3,:)*up(1,:) + 2.*ff(2,:)*kappa*up(1,:)*up(2,:) + ff(1,:)*(3-4*kappa)*rv2(:) ) /&
    !        & ( 2*kappa*up(3,:)*up(2,:) - 2.*(kappa-1)*up(2,:)*rv2(:))
    ! dv/dx
    !dpdx(2,:) = ( 2*ff(3,:) + up(2,:)*(-2.*ff(2,:)*kappa + ff(1,:)*up(2,:)*(2*kappa-1.) )) /&
    !        & ( 2*kappa*up(3,:) - 2.*(kappa-1)*rv2(:) )
    ! dP/dx
    !dpdx(3,:) = ( 2*ff(2,:)*(kappa*up(3,:) + rv2(:)) - up(2,:)*(2*ff(1,:)*kappa*up(3,:) + 2*ff(3,:)*up(1,:) + ff(1,:)*rv2(:) ) ) /&
    !        & ( 2*kappa*up(3,:) - 2.*(kappa-1)*rv2(:) )


    !print*, minval(abs( 2*kappa*up(3,:)*up(2,:) - 2.*(kappa-1)*up(2,:)*rv2(:)) )

    !dpdx(:,0:4) = 0.
    !dpdx(:,n-4:n) = 0.
end subroutine getpx

!need to get time derivatives of primative variables (from conserved ones)
subroutine getpt(dpdt,ff,up)
    use prms
    implicit none
    real, intent(out) :: dpdt(3,0:n)
    real, intent(in) :: ff(3,0:n), up(3,0:n)
    real, dimension(0:n) :: rv2

    ! drho/dt
    dpdt(1,:) = ff(1,:)
    ! dv/dt
    dpdt(2,:) = ( ff(2,:) - ff(1,:)*up(2,:) )/up(1,:)
    ! dP/dt
    dpdt(3,:) = 0.5*(g-1)*(2*ff(3,:) + up(2,:)*(ff(1,:)*up(2,:) - 2*ff(2,:)) )
end subroutine getpt

subroutine linspace(xmin,xmax,xx)
    use prms
    implicit none
    real,intent(in) :: xmin,xmax
    real,intent(out) :: xx(npert)
    integer :: i
    do i=1,npert
        xx(i) = (xmax-xmin) * real(i-1) / real(npert-1) + xmin
    end do
end subroutine linspace

subroutine logspace(xmin,xmax,xx)
    use prms
    implicit none
    real,intent(in)  :: xmin,xmax
    real,intent(out) :: xx(npert)
    call linspace(log10(xmin),log10(xmax),xx)
    xx(:) = 10.**xx(:)
end subroutine logspace
