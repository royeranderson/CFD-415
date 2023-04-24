subroutine cellarea(xmatg,ymatg,imax,jmax,Amat)
    implicit none

    integer, intent(in) :: imax, jmax
    integer :: i, j
    real(kind=8), intent(in) :: xmatg(-1:imax+2,-1:jmax+2), ymatg(-1:imax+2,-1:jmax+2)
    real(kind=8), intent(out) :: Amat(-1:imax+1,-1:jmax+1)

    Amat = 0.0

    do i=-1,(imax+1)
        do j=-1,(jmax+1)
            Amat(i,j)=abs(0.5*((xmatg(i+1,j+1)-xmatg(i,j))*(ymatg(i,j+1)-ymatg(i+1,j)) &
            -(ymatg(i+1,j+1)-ymatg(i,j))*(xmatg(i,j+1)-xmatg(i+1,j))))
            if ((i<1 .or. i>(imax-1)) .and. (j<1 .or. j>(jmax-1))) then
                Amat(i,j) = 0.0
            endif
        enddo
    enddo
    print*,'After Area'
    print*,xmatg

end subroutine cellarea