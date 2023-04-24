subroutine initial(xmatg,ymatg,imax,jmax,M_in,qmat,gmat,fmat)
    implicit none

    integer, intent(in) :: imax,jmax
    integer :: i,j
    real(kind=8), intent(in) :: xmatg(-1:imax+2,-1:jmax+2),ymatg(-1:imax+2,-1:jmax+2),&
                                M_in
    real(kind=8), intent(out) :: qmat(-1:imax+1,-1:jmax+1,4),fmat(-1:imax+1,-1:jmax+1,4),&
                                gmat(-1:imax+1,-1:jmax+1,4)
    real(kind=8) :: angle,const1,Mnew,angle1,angle2,i1,i2,i3,i4,height,const2

    const1 = 1.7857142857
    const2 = 0.7142857143
    print*,'Begin of Initial'
    print*,xmatg(-1,-1)

    ! Interior Initial
    do i = 1,(imax+1)
        do j = -1,(jmax+1)
            angle1 = atan2(ymatg(i+1,j)-ymatg(i,j),xmatg(i+1,j)-xmatg(i,j))
            print*,i
            print*,j
            print*,xmatg(-1,-1)
            angle2 = atan2(ymatg(i+1,j+1)-ymatg(i,j+1),xmatg(i+1,j+1)-xmatg(i,j+1))
            angle = 0.5*(angle1+angle2)
            height = 1/(0.5*((ymatg(i,jmax)-ymatg(i,1))+(ymatg(i+1,jmax)-ymatg(i+1,1))))
            Mnew = height*M_in

            qmat(i,j,1) = 1.0
            qmat(i,j,2) = Mnew*cos(angle)
            qmat(i,j,3) = Mnew*sin(angle)
            qmat(i,j,4) = const1 + 0.5*(Mnew**2)

            fmat(i,j,1) = Mnew*cos(angle)
            fmat(i,j,2) = (Mnew*cos(angle))**2+const2
            fmat(i,j,3) = (Mnew**2)*cos(angle)*sin(angle)
            fmat(i,j,4) = (const1 + 0.5*(Mnew**2) + const2)*Mnew*cos(angle)

            gmat(i,j,1) = Mnew*sin(angle)
            gmat(i,j,2) = (Mnew**2)*cos(angle)*sin(angle)
            gmat(i,j,3) = (Mnew*sin(angle))**2+const2
            gmat(i,j,4) = (const1 + 0.5*(Mnew**2) + const2)*Mnew*sin(angle)

            if (i > imax-1 .and. j > jmax-1) then
                qmat(i,j,1) = 0.0
                qmat(i,j,2) = 0.0
                qmat(i,j,3) = 0.0
                qmat(i,j,4) = 0.0

                fmat(i,j,1) = 0.0
                fmat(i,j,2) = 0.0
                fmat(i,j,3) = 0.0
                fmat(i,j,4) = 0.0

                gmat(i,j,1) = 0.0
                gmat(i,j,2) = 0.0
                gmat(i,j,3) = 0.0
                gmat(i,j,4) = 0.0
            endif

        if (MOD(jmax,2) .eq. 0) then 
            qmat(i,jmax/2,2) = Mnew
            qmat(i,jmax/2,3) = 0.0
        endif
        enddo
    enddo

    print*,'Before Ghost Init'
    print*,xmatg

    ! Ghost Initial
    !   Inlet
    do i = -1,0
        do j = -1,jmax+1

            qmat(i,j,1) = 1.0
            qmat(i,j,2) = M_in
            qmat(i,j,3) = 0.0
            qmat(i,j,4) = const1 + 0.5*(M_in**2)

            fmat(i,j,1) = M_in
            fmat(i,j,2) = (M_in)**2+const2
            fmat(i,j,3) = 0.0
            fmat(i,j,4) = (const1 + 0.5*(M_in**2) + const2)*M_in

            gmat(i,j,1) = 0.0
            gmat(i,j,2) = 0.0
            gmat(i,j,3) = const2
            gmat(i,j,4) = 0.0

            if (j<1 .or. j>jmax-1) then
                qmat(i,j,1) = 0.0
                qmat(i,j,2) = 0.0
                qmat(i,j,3) = 0.0
                qmat(i,j,4) = 0.0

                fmat(i,j,1) = 0.0
                fmat(i,j,2) = 0.0
                fmat(i,j,3) = 0.0
                fmat(i,j,4) = 0.0

                gmat(i,j,1) = 0.0
                gmat(i,j,2) = 0.0
                gmat(i,j,3) = 0.0
                gmat(i,j,4) = 0.0

            endif
        enddo
    enddo

    open(unit=3, file='test.plt', status='replace')
    write(3, *) 'VARIABLES', ' = "X"', ', "Y"', ', "RHO"', ', "RHO_U"' , ', "RHO_V"' , &
            ', "RHO_E"'
    write(3, *) 'ZONE', ' I=', imax+4 , ', J=', jmax+4, ', K=', 1, ', DATAPACKING=BLOCK', &
                ', VARLOCATION=', '([3, 4, 5, 6]=CELLCENTERED)'

    ! Node Locations
    write(3, *) ((xmatg(i, j), i = -1, imax+2), j = -1, jmax+2)
    write(3, *) ((ymatg(i, j), i = -1, imax+2), j = -1, jmax+2)
    ! q values
    write(3, *) ((qmat(i, j, 1), i = -1, imax+1), j = -1, jmax+1)
    write(3, *) ((qmat(i, j, 2), i = -1, imax+1), j = -1, jmax+1)
    write(3, *) ((qmat(i, j, 3), i = -1, imax+1), j = -1, jmax+1)
    write(3, *) ((qmat(i, j, 4), i = -1, imax+1), j = -1, jmax+1)

    close(3)

    open(unit=997, file='test2.plt', status='replace')
    write(997, *) ' VARIABLES = "X", "Y", "U", "V"'
    write(997, *) ' ZONE I=', imax+4, ', J=', jmax+4, ', DATAPACKING=BLOCK', &
        ', VARLOCATION=([3,4]=CELLCENTERED)'
    write(997, *) ((xmatg(i, j), i=-1,imax+2), j=-1, jmax+2)
    write(997, *) ((ymatg(i, j), i=-1,imax+2), j=-1, jmax+2)
    write(997, *) ((qmat(i,j,2), i=-1, imax+1), j=-1, jmax+1)
    write(997, *) ((qmat(i,j,3), i=-1, imax+1), j=-1, jmax+1)
    close(997)
    print*,'In Initial'
    print*,xmatg

    open(unit=10, file='grid.x', status='replace') ! Plot3D format data
    write(10,*) imax+4, jmax+4
    write(10,*) ((xmatg(i,j), i=-1,imax+2), j=-1,jmax+2), &
        ((ymatg(i,j), i=-1,imax+2), j=-1,jmax+2)

    close(10)


end subroutine initial