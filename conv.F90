MODULE conv
    use prms

contains

subroutine c_to_p(up,uc)
    IMPLICIT NONE
    REAL, INTENT(IN)  :: uc(3,0:n)
    REAL, INTENT(OUT) :: up(3,0:n)

    up(1,:) = uc(1,:)
    up(2,:) = uc(2,:)/uc(1,:)
    up(3,:) = (g-1.)*(uc(3,:) - 0.5*(uc(2,:)**2.)/uc(1,:))

end subroutine c_to_p

subroutine p_to_c(up,uc)
    IMPLICIT NONE
    REAL, INTENT(IN)  :: up(3,0:n)
    REAL, INTENT(OUT) :: uc(3,0:n)

    uc(1,:) = up(1,:)
    uc(2,:) = up(1,:)*up(2,:)
    uc(3,:) = up(3,:)/(g-1.) + 0.5*up(1,:)*(up(2,:)**2.)

end subroutine p_to_c

END MODULE conv
