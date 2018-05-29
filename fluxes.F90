module fluxes

    contains

SUBROUTINE euler_fluxes(ff, up)
  USE prms
  IMPLICIT NONE
  REAL, INTENT(IN)  :: up(3,0:n)
  REAL, INTENT(OUT) :: ff(3,0:n)
  REAL              :: alpha
  integer :: i

  ff(1,:) = up(1,:)*up(2,:)
  ff(2,:) = up(1,:)*up(2,:)*up(2,:) + up(3,:)
  ff(3,:) = up(2,:)*( up(3,:)*g/(g-1.) + 0.5*up(1,:)*up(2,:)*up(2,:))

END SUBROUTINE euler_fluxes

SUBROUTINE adjointflux(fp, fm, pc)
  USE prms
  IMPLICIT NONE
  REAL, INTENT(IN)  :: pc(3,0:n)
  REAL, INTENT(OUT) :: fp(3,0:n), fm(3,0:n)
  REAL              :: alpha

  if ( aflux.eq.'GLF') then
    !alpha = maxwavespeed
    !alpha =1.
    alpha =0.
    fp    = 0.5 * (-pc + alpha*pc)
    fm    = 0.5 * (-pc - alpha*pc)
  else if (aflux.eq.'AV') then
    alpha = eps+maxwavespeed
    fp    = 0.5 * (-pc + alpha*pc)
    fm    = 0.5 * (-pc - alpha*pc)
  end if

END SUBROUTINE adjointflux

SUBROUTINE numflux(fp, fm, uc, f)
    USE prms
    IMPLICIT NONE
    REAL, INTENT(IN)  :: uc(3,0:n), f(3,0:n)
    REAL, INTENT(OUT) :: fp(3,0:n), fm(3,0:n)
    REAL              :: alpha
    integer :: i

    if (FFlux .eq. 'GLF') then
        alpha = maxwavespeed  !MAXVAL(ABS(u))
        fp    = 0.5 * (f + alpha*uc)
        fm    = 0.5 * (f - alpha*uc)
    else if( FFLUX .eq. 'LLF') then
        do i = 0,n
            alpha   = 1+sqrt(1./uc(1,i))
            fp(:,i) = 0.5 * (f(:,i) + alpha*uc(:,i))
            fm(:,i) = 0.5 * (f(:,i) - alpha*uc(:,i))
        end do
    else if (FFlux .eq. 'AV') then
        alpha = eps+maxwavespeed
        fp    = 0.5 * (f + alpha*uc)
        fm    = 0.5 * (f - alpha*uc)
    end if
END SUBROUTINE numflux

end module fluxes
