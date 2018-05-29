PROGRAM main

    USE doublePrecision
    USE prms
    USE conv
    USE rhs
    IMPLICIT NONE

    real(KIND = dp), DIMENSION(:),   ALLOCATABLE :: x 
    real(KIND = dp), DIMENSION(:,:),   ALLOCATABLE :: up, uc, frp, frm, fp, fm, f, ff, pc, matvec
    real(KIND = dp), DIMENSION(:,:,:),   ALLOCATABLE :: urk, prk, ucsol, A
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
            fm(3,0:n), ff(3,0:n), pc(3,0:n), matvec(3,0:n) )
    allocate( urk(3,0:n,1:2), prk(3,0:n,1:2), A(3,3,0:n), ucsol(0:nsteps,3,0:n) )

    call getic(x,uc,up,ucsol)

    flag='F'
    !bootstrap simulation
    if (.false.) then
    do nt = 1,int(nsteps/tend),1
    !do nt = 1,int(nsteps/tend),1
        if (mod(nt,outfile).eq.0) print*, 'Bootstrap: Timestep = ', nt

        if (FT.eq.'Euler') then
            call c_to_p(up,uc)
            call fluxes(ff,up)
            call numflux(fp, fm, uc, ff)
            call recon(fp,fm,f)
            call rhside(fp, f, n, hi)
            if (enable_viscosity.eq.'True') call viscosity(uc,fp,flag)
            uc = uc + dt * fp
        else if (FT.eq.'RK3') then
            !substep 1
            call c_to_p(up,uc)
            call fluxes(ff,up)
            call numflux(fp, fm, uc, ff)
            call recon(fp,fm,f)
            call rhside(fp, f, n, hi)
            if (enable_viscosity.eq.'True') call viscosity(uc,fp,flag)
            urk(:,:,1)= uc(:,:) + dt * fp(:,:)
            
            !substep 2
            call c_to_p(up,urk(:,:,1))
            call fluxes(ff,up)
            call numflux(fp, fm, urk(:,:,1), ff)
            call recon(fp,fm,f)
            call rhside(fp, f, n, hi)
            if (enable_viscosity.eq.'True') call viscosity(urk(:,:,1),fp,flag)
            urk(:,:,2) = rk(1,1) * uc(:,:) + rk(1,2) * (urk(:,:,1) + dt * fp(:,:))

            !substep 3
            call c_to_p(up,urk(:,:,2))
            call fluxes(ff,up)
            call numflux(fp, fm, urk(:,:,2), ff)
            call recon(fp,fm,f)
            call rhside(fp, f, n, hi)
            if (enable_viscosity.eq.'True') call viscosity(urk(:,:,2),fp,flag)
            uc(:,:) = rk(2,1) * uc(:,:) + rk(2,2) * (urk(:,:,2) + dt * fp(:,:))
        end if
    end do
    
    call moveback(uc)
    call output(uc(1,:),x,'rho',0)
    end if
    
    print*, 'Forward simulation'; flag='F'
    do nt = 1,nsteps,1
        if (mod(nt,outfile).eq.0) print*, 'Forward: Timestep = ', nt

        if (FT.eq.'Euler') then
            if (SI.eq.'WENO') then
                call c_to_p(up,uc)
                call fluxes(ff,up)
                call numflux(fp, fm, uc, ff)
                call recon(fp,fm,f)
                call rhside(fp, f, n, hi)
                if (enable_viscosity.eq.'True') call viscosity(uc,fp,flag)
            uc = uc + dt * fp
        else if (FT.eq.'RK3') then
            !substep 1
            call c_to_p(up,uc)
            call fluxes(ff,up)
            call numflux(fp, fm, uc, ff)
            call recon(fp,fm,f)
            call rhside(fp, f, n, hi)
            if (enable_viscosity.eq.'True') call viscosity(uc,fp,flag)
            urk(:,:,1)= uc(:,:) + dt * fp(:,:)
            
            !substep 2
            call c_to_p(up,urk(:,:,1))
            call fluxes(ff,up)
            call numflux(fp, fm, urk(:,:,1), ff)
            call recon(fp,fm,f)
            call rhside(fp, f, n, hi)
            if (enable_viscosity.eq.'True') call viscosity(urk(:,:,1),fp,flag)
            urk(:,:,2) = rk(1,1) * uc(:,:) + rk(1,2) * (urk(:,:,1) + dt * fp(:,:))

            !substep 3
            call c_to_p(up,urk(:,:,2))
            call fluxes(ff,up)
            call numflux(fp, fm, urk(:,:,2), ff)
            call recon(fp,fm,f)
            call rhside(fp, f, n, hi)
            if (enable_viscosity.eq.'True') call viscosity(urk(:,:,2),fp,flag)
            uc(:,:) = rk(2,1) * uc(:,:) + rk(2,2) * (urk(:,:,2) + dt * fp(:,:))
        end if
        
        ucsol(nt,:,:) = uc(:,:)
        call output(uc(1,:),x,'rho',nt)
        if (maxval(abs(uc)).gt.100 ) stop 'Crash Forward'
    end do
    call output(uc(1,:),x,'rho',-1)
    call pas(uc(1,:))

    !print*, 'L1 Error = ', sum(    abs( uc(1,:) - 3*(exp(-( x(:)-tend)**2.)+1) )) / real(n)
    !print*, 'L0 Error = ', maxval( abs( uc(1,:) - 3*(exp(-( x(:)-tend)**2.)+1) ))
    !open(1,file='D/l1err.out')
    !    write(1,*) sum(    abs( uc(1,:) - 3*(exp(-( x(:)-tend)**2.)+1) ))/real(n)
    !close(1)

    !adjoint solution
    if (enable_adjoint.eq.'False') stop

    pc = uc; flag='A'
    do nt = nsteps,0,-1
        if (mod(nt,outfile).eq.0) print*, 'Adjoint: Timestep =', nt
        call getA(A,ucsol(nt,:,:))
       
        if (BT.eq.'Euler') then
            call adjointflux(fp, fm, pc)
            call recon(fp,fm,f)
            call rhside(fp, f, n, hi)
            call mv(A,fp,matvec)
            if (enable_viscosity.eq.'True') call viscosity(pc,matvec,flag)
            pc = pc + dt * matvec
        else if (BT.eq.'RK3') then
            !substep 1
            call adjointflux(fp, fm, pc)
            call recon(fp,fm,f)
            call rhside(fp, f, n, hi)
            call mv(A,fp,matvec)
            if (enable_viscosity.eq.'True') call viscosity(pc,matvec,flag)
            prk(:,:,1) = pc(:,:) + dt * matvec(:,:)

            !substep 2
            call adjointflux(fp, fm, prk(:,:,1))
            call recon(fp,fm,f)
            call rhside(fp, f, n, hi)
            call mv(A,fp,matvec)
            if (enable_viscosity.eq.'True') call viscosity(prk(:,:,1),matvec,flag)
            prk(:,:,2) = rk(1,1) * pc(:,:) + rk(1,2) * (prk(:,:,1) + dt * matvec(:,:))
            
            !substep 3
            call adjointflux(fp, fm, prk(:,:,2))
            call recon(fp,fm,f)
            call rhside(fp, f, n, hi)
            call mv(A,fp,matvec)
            if (enable_viscosity.eq.'True') call viscosity(prk(:,:,2),matvec,flag)
            pc(:,:) = rk(2,1) * pc(:,:) + rk(2,2) * (prk(:,:,2) + dt * matvec(:,:)) 
        end if

        call output(pc(1,:),x,'adjrho',nt)
        call output(pc(2,:),x,'adjvel',nt)
        call output(pc(3,:),x,'adjE',nt)
        if (maxval(abs(pc)).gt.1000 ) stop 'Crash Adj'
    end do
    call output(pc(1,:),x,'adjrho',-1)
    call output(pc(2,:),x,'adjvel',-1)
    call output(pc(3,:),x,'adjE'  ,-1)
    
    call out_targ(pc(1,:),x,'plateauL.',-0.75)
    call out_targ(pc(1,:),x,'plateauR.',1.0)

    stop
END PROGRAM main

!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine recon(fp,fm,f)
  use prms
  use doubleprecision
  IMPLICIT NONE
  real(KIND = dp), INTENT(IN) :: fp(3,0:n), fm(3,0:n)
  real(KIND = dp) :: frp(3,0:n), frm(3,0:n)
  real(KIND = dp), INTENT(OUT) :: f(3,0:n)
  integer :: i

    do i = 1,3
        call RECONSTRUCT_FUNCTION (frp(i,:), fp(i,:), n,  1)
        call RECONSTRUCT_FUNCTION (frm(i,:), fm(i,:), n, -1)
        f(i,:) = frp(i,:) + frm(i,:)
    end do
end subroutine recon

SUBROUTINE moveback(u)
  USE doublePrecision
  USE prms
  IMPLICIT NONE
  REAL(KIND = dp), INTENT(INOUT)    :: u(3,0:n)
  REAL(KIND = dp)                   :: udum(3,0:n)
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

SUBROUTINE out_targ(u, x, str, targ)
  USE doublePrecision
  USE prms
  IMPLICIT NONE
  REAL(KIND = dp), INTENT(IN)   :: u(0:n)
  REAL(KIND = dp), INTENT(IN)   :: x(0:n)
  REAL, INTENT(IN)              :: targ
  CHARACTER(LEN = *), INTENT(IN) :: str
  character(len=30) :: fn
  INTEGER :: i, loc


    loc = minloc( abs(x-targ), dim=1 )
    write(fn,"('O/targ_',A,'out')") str
    open(2,file=fn,position="append")
    write(2,*) n,u(loc)
    close(2)
END SUBROUTINE out_targ

SUBROUTINE pas(u)
  USE doublePrecision
  USE prms
  IMPLICIT NONE
  REAL(KIND = dp), INTENT(IN) :: u(0:n)
  character(len=30) :: fn
  INTEGER :: i, counter

    counter = 0
    do i = 0,n
        if ( (2*amp*0.95 > u(i)) .and. (u(i) > amp*1.05) ) counter = counter + 1
    end do

    open(2,file='O/PAS.out',position="append")
    write(2,*) n,counter,counter*h
    close(2)
END SUBROUTINE pas
