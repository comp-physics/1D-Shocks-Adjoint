module io

    contains

SUBROUTINE output (u, x, str, lt)
  USE prms
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: lt
  REAL, INTENT(IN) :: u(0:n)
  REAL, INTENT(IN) :: x(0:n)
  CHARACTER(LEN = *), INTENT(IN) :: str
  character(len=30) :: fn
  character(len=80) :: header
  INTEGER :: i, iter

  if ( n.gt.1000) then
      iter=10
    else
      iter=1
  end if

  if ( ( mod(lt,outfile).eq.0 ) .or. lt.eq.nsteps) then
    write(fn,"('D/',A,I10.10,'.csv')") str,lt
    !write(header,"('ZONE STRANDID=1, SOLUTIONTIME=',I4.4 )") lt
    write(header,"('xdat,   udat')") 
    OPEN  (2, FILE = fn)
    write(2,*) header
    do i = 0,n,iter
        !WRITE (2,*) x(i),',',u(i)
        WRITE (2,"(F10.5,A,F10.5)") x(i),',',u(i)
    end do
    CLOSE (2)
  end if

  if (lt .eq. -1) then
    write(fn,"('O/final_',A,'_',I5.5,'.csv')") str, n
    OPEN  (2, FILE = fn)
    write(header,"('ColumnX   ColumnY')") 
    write(2,*) header
    do i = 0,n,iter
        !WRITE (2,*) x(i),u(i)
        WRITE (2,"(2F10.5)") x(i),u(i)
    end do
    CLOSE (2)
 
  end if
END SUBROUTINE output

SUBROUTINE output_diag(pc, uc, x, str, lt)
  USE prms
  USE matrices
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: lt
  REAL, INTENT(INOUT) :: pc(3,0:n)
  REAL, INTENT(IN) :: uc(3,0:n)
  REAL, INTENT(IN) :: x(0:n)
  real :: psave(3,0:n)
  CHARACTER(LEN = *), INTENT(IN) :: str
  character(len=30) :: fn
  character(len=80) :: header
  INTEGER :: i, j, iter

  psave = pc
  call q_to_p(pc,uc)

  if ( n.gt.1000) then
      iter=10
    else
      iter=1
  end if

  if ( ( mod(lt,outfile).eq.0 ) .or. lt.eq.nsteps) then
    do j = 1,3
    write(fn,"('D/',A,I1.1,'_',I10.10,'.csv')") str,j,lt
    write(header,"('xdat,   udat')") 
    OPEN  (2, FILE = fn)
    write(2,*) header
    do i = 0,n,iter
        WRITE (2,"(F10.5,A,F10.5)") x(i),',',pc(j,i)
    end do
    CLOSE (2)
    end do
  end if
  pc = psave
END SUBROUTINE output_diag

SUBROUTINE out_targ(u, x, str, targ)
  USE prms
  IMPLICIT NONE
  REAL, INTENT(IN)   :: u(0:n)
  REAL, INTENT(IN)   :: x(0:n)
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
  USE prms
  IMPLICIT NONE
  REAL, INTENT(IN) :: u(0:n)
  character(len=30) :: fn
  INTEGER :: i, counter

    counter = 0
    do i = 0,n
        if ( (2*amp*0.95 > u(i)) .and. (u(i) > amp*1.05) ) counter = counter + 1
    end do

    open(2,file='O/PAS.out',position="append")
    write(2,*) n,counter,counter*dx
    close(2)
END SUBROUTINE pas


end module io
