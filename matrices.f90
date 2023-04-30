subroutine matrices(xmat,ymat,imax,jmax,M_in,po_inf,p_inf,p_ex,rho_inf,T_inf,a_inf)
    implicit none

    integer, intent(in) :: imax, jmax
    integer :: i, j, iterations
    real(kind=8), intent(in) :: xmat(imax,jmax),ymat(imax,jmax), M_in,po_inf,rho_inf,T_inf,a_inf,p_inf,p_ex
    real(kind=8) :: xmatge(-1:imax+2,-1:jmax+2), ymatge(-1:imax+2,-1:jmax+2),&
                    xmatg(-1:imax+2,-1:jmax+2), ymatg(-1:imax+2,-1:jmax+2)
    real(kind=8) :: Amat(-1:imax+1,-1:jmax+1),qmat(-1:imax+1,-1:jmax+1,4),&
                    Rmat(1:imax-1,1:jmax-1,4),Dmat(1:imax-1,1:jmax-1,4),&
                    fmat(-1:imax+1,-1:jmax+1,4),gmat(-1:imax+1,-1:jmax+1,4),&
                    BCmat(-1:imax+1,-1:jmax+1,4),wall_ang(-1:imax+1,2),alfmat(-1:imax+1,-1:jmax+1),&
                    facemat(-1:imax+1,-1:jmax+1,8),normmat(-1:imax+1,-1:jmax+1,4,2)

    real(kind=8) :: delxi, deleta, delx, dely

    ! Find delxi and deleta
    delxi = 1.0_8/(imax-1)
    deleta = 1.0_8/(jmax-1)
    delx = delxi*xmat(imax,1)
    dely = deleta*ymat(1,jmax)

    ! Create Matrix With Ghost Cells
    call ghostcell(imax,jmax,xmat,xmatge)
    call ghostcell(imax,jmax,ymat,ymatge)

    ! Fill in Ghost Nodes
    call ghostfill(xmatge,ymatge,xmat,ymat,imax,jmax,delx,dely,xmatg,ymatg)

    ! Create a Matrix of Cell Areas
    call cellarea(xmatg,ymatg,imax,jmax,Amat)

    ! Create Q-matrix Initial Guess
    call initial(xmatg,ymatg,imax,jmax,M_in,qmat,gmat,fmat,wall_ang,alfmat,facemat,normmat)

    ! call eulerBCs(imax,jmax,qmat,gmat,fmat,po_inf,p_inf,p_ex,rho_inf,T_inf,M_in,wall_ang,a_inf,alfmat)

    ! call dissapation(imax,jmax,qmat,facemat,normmat,po_inf,rho_inf,a_inf,Dmat)

    ! Iterate to Solution
    iterations = 1
    do while (iterations<100)
        print*,iterations
        call RK(imax,jmax,po_inf,rho_inf,a_inf,p_inf,p_ex,T_inf,qmat,Amat,Rmat,Dmat,normmat,facemat,alfmat,fmat,gmat,wall_ang,M_in)
        
        iterations = iterations+1
        ! print*,size(qmat,1)
        ! print*,size(qmat,2)
        print*,minval(Rmat(:,:,1))
        print*,minval(Rmat(:,:,2))
        print*,minval(Rmat(:,:,3))
        print*,minval(Rmat(:,:,4))
        ! print*,qmat(:,:,1)
    enddo
    !print*,iterations

    ! Plot
    call plot(xmatg,ymatg,imax,jmax,qmat)

end subroutine matrices