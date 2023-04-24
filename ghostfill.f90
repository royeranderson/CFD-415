subroutine ghostfill(xmatge,ymatge,xmat,ymat,imax,jmax,delx,dely,xmatg,ymatg)
    implicit none
    integer, intent(in) :: imax,jmax
    real(kind=8), intent(in) :: xmatge(-1:imax+2,-1:jmax+2), ymatge(-1:imax+2,-1:jmax+2),&
                                xmat(imax,jmax),ymat(imax,jmax), delx, dely
    real(kind=8), intent(out) :: xmatg(-1:imax+2,-1:jmax+2), ymatg(-1:imax+2,-1:jmax+2)
    integer :: ghostfiller, i, j

    ghostfiller = 2
    xmatg = xmatge
    ymatg = ymatge
    do i=-1,0
        do j = 1,jmax
            xmatg(i,j) = -1*ghostfiller*delx
            ymatg(i,j) = ymat(1,j)
        enddo
        ghostfiller = ghostfiller-1
    enddo
    ghostfiller = 1
    do i=imax+1,imax+2
        do j = 1,jmax
            xmatg(i,j) = (ghostfiller*delx) + xmat(imax,j)
            ymatg(i,j) = ymat(imax,j)
        enddo
        ghostfiller = ghostfiller+1
    enddo
    ghostfiller = 2
    do j=-1,0
        do i = 1,imax
            ymatg(i,j) = ymat(i,1)-1*ghostfiller*dely
            xmatg(i,j) = xmat(i,1)
        enddo
        ghostfiller = ghostfiller-1
    enddo
    ghostfiller = 1
    do j=jmax+1,jmax+2
        do i = 1,imax
            ymatg(i,j) = (ghostfiller*dely) + ymat(i,jmax)
            xmatg(i,j) = xmat(i,jmax)
        enddo
        ghostfiller = ghostfiller+1
    enddo

    xmatg(-1,-1) = -2*delx
    xmatg(-1,0) = -2*delx
    xmatg(0,-1) = -1*delx
    xmatg(0,0) = -1*delx
    xmatg(imax+2,-1) = xmat(imax,1)+2*delx
    xmatg(imax+2,0) = xmat(imax,1)+2*delx
    xmatg(imax+1,-1) = xmat(imax,1)+1*delx
    xmatg(imax+1,0) = xmat(imax,1)+1*delx
    xmatg(-1,jmax+2) = -2*delx
    xmatg(-1,jmax+1) = -2*delx
    xmatg(0,jmax+2) = -1*delx
    xmatg(0,jmax+1) = -1*delx
    xmatg(imax+2,jmax+2) = xmat(imax,1)+2*delx
    xmatg(imax+2,jmax+1) = xmat(imax,1)+2*delx
    xmatg(imax+1,jmax+2) = xmat(imax,1)+1*delx
    xmatg(imax+1,jmax+1) = xmat(imax,1)+1*delx

    ymatg(-1,-1) = -2*dely
    ymatg(-1,0) = -1*dely
    ymatg(0,-1) = -2*dely
    ymatg(0,0) = -1*dely
    ymatg(imax+2,-1) = -2*dely
    ymatg(imax+2,0) = -1*dely
    ymatg(imax+1,-1) = -2*dely
    ymatg(imax+1,0) = -1*dely
    ymatg(-1,jmax+2) = ymat(imax,jmax)+2*dely
    ymatg(-1,jmax+1) = ymat(imax,jmax)+1*dely
    ymatg(0,jmax+2) = ymat(imax,jmax)+2*dely
    ymatg(0,jmax+1) = ymat(imax,jmax)+1*dely
    ymatg(imax+2,jmax+2) = ymat(imax,jmax)+2*dely
    ymatg(imax+2,jmax+1) = ymat(imax,jmax)+1*dely
    ymatg(imax+1,jmax+2) = ymat(imax,jmax)+2*dely
    ymatg(imax+1,jmax+1) = ymat(imax,jmax)+1*dely

    print*,'In Ghostfill'
    print*,xmatg

end subroutine ghostfill