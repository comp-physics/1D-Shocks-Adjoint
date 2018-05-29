MODULE prms

    REAL, PARAMETER :: pi      = 4.*ATAN(1.)
    REAL, PARAMETER :: rk(2,2) = RESHAPE( &
        [0.75, (1./3.), (0.25), (2./3.)], SHAPE(rk))

    real, parameter :: g=1.4    !ratio of specific heats
    INTEGER :: n                !number of spatial grid points
    REAL :: L                !length of computational domain
    REAL :: tend                !dimensionless end time of simulation
    REAL :: cfl                 !courant nymber
    REAL :: beta                !exponent of dissipation coefficient
    REAL :: amp                 !amplitude of initial condition
    INTEGER :: outfile          !how often to output files and print timesteps
    CHARACTER(LEN = 20) :: FT, BT           !time marching schemes for forward and backward simulations
    CHARACTER(LEN = 20) :: SI               !spatial integration
    CHARACTER(LEN = 20) :: FFlux, AFlux     !flux functions for forward and adjoint simulations
    CHARACTER(LEN = 20) :: enable_adjoint   !flag for adjoint simulation
    CHARACTER(LEN = 20) :: enable_viscosity !flag for added viscosity
    CHARACTER(LEN = 20) :: enable_f_viscosity !flag for added viscosity to forward equations
    CHARACTER(LEN = 20) :: enable_diag      !flag for adjoint diagonalization
    CHARACTER(LEN = 20) :: enable_acc       !flag for checking adjoint sensitivity
 
    integer :: nt, nsteps, nfor, which_ic, npert
    real    :: dx, dxi, t, dt
    real :: vinf, pres
    REAL :: maxwavespeed        !typically called alpha .. = abs(v) + c
    REAL :: eps                 !dissipation coefficeint 

contains

subroutine getprms
    integer :: i

    open(1,file='Input/weno.in')
    read(1,*) n;                print*, 'Nx = ',n
    read(1,*) L;                print*, 'L = ', L
    read(1,*) tend;             print*, 'Tend = ', Tend
    read(1,*) cfl;              print*, 'CFL = ', CFL
    read(1,*) beta;             print*, 'beta = ', beta
    read(1,*) amp;              print*, 'amp = ', amp
    read(1,*) outfile;          print*, 'Output frequency = ', outfile
    read(1,*) FT;               print*, 'FT: ', FT
    read(1,*) BT;               print*, 'BT: ', BT
    read(1,*) FFlux;            print*, 'Forward Flux: ', FFlux
    read(1,*) AFlux;            print*, 'Adjoint Flux: ', AFlux
    read(1,*) enable_adjoint;   print*, 'Adjoint enabled: ', enable_adjoint
    read(1,*) enable_viscosity; print*, 'Viscosity enabled: ', enable_viscosity
    read(1,*) enable_f_viscosity; print*, 'Viscosity enabled: ', enable_viscosity
    read(1,*) which_ic;         print*, 'Initial condition: ', which_ic
    read(1,*) SI;               print*, 'SI: ', SI
    read(1,*) enable_diag;      print*, 'Diagonalization enabled: ', enable_diag
    read(1,*) enable_acc;       print*, 'Sensitivity check enabled: ', enable_acc
    close(1)

    if (which_ic .eq. 1 .or. which_ic .eq. 2 .or. which_ic .eq.5) then
        vinf = 1; pres=1./g
        maxwavespeed = vinf + sqrt(g*pres/amp) !for entropy wave of speed 1
    else if (which_ic .eq. 3 .or. which_ic .eq. 33) then
        maxwavespeed = 1 + sqrt(g/0.125)
    else if (which_ic .eq. 4) then
        maxwavespeed = 0.2+sqrt(g*5.1/amp)
    else if (which_ic .eq. 6) then
        maxwavespeed = 2+sqrt(1./amp)
    else if (which_ic .eq. 7) then
        maxwavespeed = 1.1+sqrt(g*0.2/0.5)
    else if (which_ic .eq. 8) then
        maxwavespeed = 1.1+sqrt(g/0.125)
    end if

    npert = 10

    dx = L / real(n)
    dxi = 1.d0 / dx

    !eps = 0.0005
    !eps = 0.5*dx*maxwavespeed
    eps = dx**beta
    !eps = maxval( (/ dx**beta, 0.5*dx*maxwavespeed /) )
    !viscosity coefficient
    
    if ( (AFlux.eq.'AV') .or. (FFlux.eq.'AV') .or. (enable_viscosity.eq.'True' .and. eps>0.00001) &
                & .or. (enable_f_viscosity.eq.'True' .and. eps>0.000001) ) then
        !dt = minval( (/ cfl*dx/(maxwavespeed*eps), cfl*dx/maxwavespeed, cfl*dx/(maxwavespeed+eps) /) )
        dt = minval( (/ cfl*dx*dx/(maxwavespeed*eps), cfl*dx/maxwavespeed /) )
        !dt = cfl*dx/maxwavespeed
    else 
        dt = cfl * dx/maxwavespeed
    end if

    nsteps = INT(tend / dt)
    dt = tend / real(nsteps)
    print '(" Number of time steps = ", I7, " with dt = ", F16.12)', nsteps, dt
end subroutine getprms

END MODULE prms
