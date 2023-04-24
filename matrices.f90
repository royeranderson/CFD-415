subroutine matrices(xmat,ymat,imax,jmax,M_in,po_inf,rho_inf,T_inf)
    implicit none

    integer, intent(in) :: imax, jmax
    integer :: i, j
    real(kind=8), intent(in) :: xmat(imax,jmax),ymat(imax,jmax), M_in,po_inf,rho_inf,T_inf
    real(kind=8) :: xmatge(-1:imax+2,-1:jmax+2), ymatge(-1:imax+2,-1:jmax+2),&
                    xmatg(-1:imax+2,-1:jmax+2), ymatg(-1:imax+2,-1:jmax+2)
    real(kind=8) :: Amat(-1:imax+1,-1:jmax+1),qmat(-1:imax+1,-1:jmax+1,4),&
                    Rmat(-1:imax+1,-1:jmax+1),Dmat(-1:imax+1,-1:jmax+1),&
                    fmat(-1:imax+1,-1:jmax+1,4),gmat(-1:imax+1,-1:jmax+1,4),&
                    BCmat(-1:imax+1,-1:jmax+1,4)
    real(kind=8) :: delxi, deleta, delx, dely

    ! Find delxi and deleta
    delxi = 1.0/(imax-1)
    deleta = 1.0/(jmax-1)
    delx = delxi*xmat(imax,1)
    dely = deleta*ymat(1,jmax)

    ! Create Matrix With Ghost Cells
    call ghostcell(imax,jmax,xmat,xmatge)
    call ghostcell(imax,jmax,ymat,ymatge)

    ! Fill in Ghost Nodes
    call ghostfill(xmatge,ymatge,xmat,ymat,imax,jmax,delx,dely,xmatg,ymatg)

    ! Create a Matrix of Cell Areas
    call cellarea(xmatg,ymatg,imax,jmax,Amat)
    print*,'Matrices'
    print*,xmatg

    ! Create Q-matrix Initial Guess
    call initial(xmatg,ymatg,imax,jmax,M_in,qmat,gmat,fmat)

    ! Apply BCs to Q Matrix
    call eulerBCs(imax,jmax,qmat,gmat,fmat,po_inf,rho_inf,T_inf)

    ! Create a Matrix of Cell Residuals
    call residuals()

    ! Create a Matrix of Cell Dissapations
    call dissapation()

end subroutine matrices