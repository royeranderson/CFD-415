subroutine matrices(xmat,ymat,imax,jmax,M_in,po_inf,p_inf,p_ex,rho_inf,T_inf,a_inf)
    implicit none

    integer, intent(in) :: imax, jmax
    integer :: i, j
    real(kind=8), intent(in) :: xmat(imax,jmax),ymat(imax,jmax), M_in,po_inf,rho_inf,T_inf,a_inf,p_inf,p_ex
    real(kind=8) :: xmatge(-1:imax+2,-1:jmax+2), ymatge(-1:imax+2,-1:jmax+2),&
                    xmatg(-1:imax+2,-1:jmax+2), ymatg(-1:imax+2,-1:jmax+2)
    real(kind=8) :: Amat(-1:imax+1,-1:jmax+1),qmat(-1:imax+1,-1:jmax+1,4),&
                    Rmat(-1:imax+1,-1:jmax+1),Dmat(-1:imax+1,-1:jmax+1),&
                    fmat(-1:imax+1,-1:jmax+1,4),gmat(-1:imax+1,-1:jmax+1,4),&
                    BCmat(-1:imax+1,-1:jmax+1,4),wall_ang(1:imax-1,2),alfmat(-1:imax+1,-1:jmax+1)
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

    ! Create Q-matrix Initial Guess
    call initial(xmatg,ymatg,imax,jmax,M_in,qmat,gmat,fmat,wall_ang,alfmat)

    ! Apply BCs to Q Matrix
    call eulerBCs(imax,jmax,qmat,gmat,fmat,po_inf,p_inf,p_ex,rho_inf,T_inf,M_in,wall_ang,a_inf,alfmat)

    open(unit=997, file='test2.plt', status='replace')
    write(997, *) ' VARIABLES = "X", "Y", "U", "V"'
    write(997, *) ' ZONE I=', imax+4, ', J=', jmax+4, ', DATAPACKING=BLOCK', &
        ', VARLOCATION=([3,4]=CELLCENTERED)'
    write(997, *) ((xmatg(i, j), i=-1,imax+2), j=-1, jmax+2)
    write(997, *) ((ymatg(i, j), i=-1,imax+2), j=-1, jmax+2)
    write(997, *) ((qmat(i,j,2), i=-1, imax+1), j=-1, jmax+1)
    write(997, *) ((qmat(i,j,3), i=-1, imax+1), j=-1, jmax+1)
    close(997)

    open(unit=10, file='grid.x', status='replace') ! Plot3D format data
    write(10,*) imax+4, jmax+4
    write(10,*) ((xmatg(i,j), i=-1,imax+2), j=-1,jmax+2), &
        ((ymatg(i,j), i=-1,imax+2), j=-1,jmax+2)

    close(10)

    open(unit=3, file='test.plt', status='replace')
    write(3, *) 'VARIABLES', ' = "X"', ', "Y"', ', "RHO"', ', "RHO_U"' , ', "RHO_V"' , &
            ', "RHO_E"'
    write(3, *) 'ZONE', ' I=', imax , ', J=', jmax, ', K=', 1, ', DATAPACKING=BLOCK', &
                ', VARLOCATION=', '([3, 4, 5, 6]=CELLCENTERED)'

    ! Node Locations
    write(3, *) ((xmatg(i, j), i = 1, imax), j = 1, jmax)
    write(3, *) ((ymatg(i, j), i = 1, imax), j = 1, jmax)
    ! q values
    write(3, *) ((qmat(i, j, 1), i = 1, imax-1), j = 1, jmax-1)
    write(3, *) ((qmat(i, j, 2), i = 1, imax-1), j = 1, jmax-1)
    write(3, *) ((qmat(i, j, 3), i = 1, imax-1), j = 1, jmax-1)
    write(3, *) ((qmat(i, j, 4), i = 1, imax-1), j = 1, jmax-1)

    close(3)

    ! Create a Matrix of Cell Residuals
    call residuals()

    ! Create a Matrix of Cell Dissapations
    call dissapation()

end subroutine matrices