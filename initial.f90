subroutine initial(xmatg,ymatg,imax,jmax,M_in,qmat,gmat,fmat,wall_ang,alfmat)
    implicit none

    integer, intent(in) :: imax,jmax
    integer :: i,j
    real(kind=8), intent(in) :: xmatg(-1:imax+2,-1:jmax+2),ymatg(-1:imax+2,-1:jmax+2),&
                                M_in
    real(kind=8), intent(out) :: qmat(-1:imax+1,-1:jmax+1,4),fmat(-1:imax+1,-1:jmax+1,4),&
                                gmat(-1:imax+1,-1:jmax+1,4),wall_ang(1:imax-1,2),alfmat(-1:imax+1,-1:jmax+1)
    real(kind=8) :: angle,const1,Mnew,angle1,angle2,i1,i2,i3,i4,height,const2

    const1 = 1.7857142857
    const2 = 0.7142857143

    ! Interior Initial
    do i = -1,(imax+1)
        do j = -1,(jmax+1)
            angle1 = atan2(ymatg(i+1,j)-ymatg(i,j),xmatg(i+1,j)-xmatg(i,j))
            if ((i>=1 .and. i<imax) .and. j==1) then
                wall_ang(i,1) = angle1
            elseif ((i>=1 .and. i<imax) .and. j==(jmax-1)) then
                wall_ang(i,2) = angle1
            endif
            angle2 = atan2(ymatg(i+1,j+1)-ymatg(i,j+1),xmatg(i+1,j+1)-xmatg(i,j+1))
            angle = 0.5*(angle1+angle2)
            alfmat(i,j) = angle
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

            if ((i > imax-1 .or. i<1) .and. (j > jmax-1 .or. j<1)) then
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

end subroutine initial