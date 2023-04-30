subroutine ghostcell(imax,jmax,centermat,ghostmat)
    implicit none
    
    integer, intent(in) :: imax, jmax
    integer :: i, j
    real(kind=8), intent(in) :: centermat(imax,jmax)
    real(kind=8), intent(out) :: ghostmat(-1:imax+2,-1:jmax+2)

    ghostmat = 0.0_8
    ghostmat(1:imax,1:jmax) = centermat

end subroutine ghostcell