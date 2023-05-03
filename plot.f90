subroutine plot(xmatg,ymatg,imax,jmax,qmat,resplot,itplot)
    implicit none
    real(kind=8),intent(in) :: xmatg(-1:imax+2,-1:jmax+2), ymatg(-1:imax+2,-1:jmax+2),&
                                qmat(-1:imax+1,-1:jmax+1,4)
    real(kind=8),intent(in) :: resplot(199),itplot(199)
    integer,intent(in) :: imax,jmax
    integer :: i,j


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
    write(3, *) 'VARIABLES', ' = "X"', ', "Y"', ', "RHO"', ', "U"' , ', "V"' , &
            ', "RHO_E"'
    write(3, *) 'ZONE', ' I=', imax , ', J=', jmax, ', K=', 1, ', DATAPACKING=BLOCK', &
                ', VARLOCATION=', '([3, 4, 5, 6]=CELLCENTERED)'

    ! Node Locations
    write(3, *) ((xmatg(i, j), i = 1, imax), j = 1, jmax)
    write(3, *) ((ymatg(i, j), i = 1, imax), j = 1, jmax)
    ! q values
    write(3, *) ((qmat(i, j, 1), i = 1, imax-1), j = 1, jmax-1)
    write(3, *) (((qmat(i, j, 2)/qmat(i, j, 1)), i = 1, imax-1), j = 1, jmax-1)
    write(3, *) (((qmat(i, j, 3)/qmat(i, j, 1)), i = 1, imax-1), j = 1, jmax-1)
    write(3, *) ((qmat(i, j, 4), i = 1, imax-1), j = 1, jmax-1)

    close(3)

    ! Open a pipe to Gnuplot and send commands
    open(unit=10, file='res.x', status='replace') ! Plot3D format data
    write(10,*) 199
    write(10,*) ((itplot(i)), i=1,199), &
        ((resplot(i)), i=1,199)
    close(10)

end subroutine plot