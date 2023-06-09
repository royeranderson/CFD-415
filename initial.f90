subroutine initial(xmatg,ymatg,imax,jmax,M_in,qmat,gmat,fmat,wall_ang,alfmat,facemat,normmat,rho_inf,po_inf,a_inf)
    implicit none

    integer, intent(in) :: imax,jmax
    integer :: i,j
    real(kind=8), intent(in) :: xmatg(-1:imax+2,-1:jmax+2),ymatg(-1:imax+2,-1:jmax+2),&
                                M_in,rho_inf,po_inf,a_inf
    real(kind=8), intent(out) :: qmat(-1:imax+1,-1:jmax+1,4),fmat(-1:imax+1,-1:jmax+1,4),&
                                gmat(-1:imax+1,-1:jmax+1,4),wall_ang(-1:imax+1,2),alfmat(-1:imax+1,-1:jmax+1)
    real(kind=8),intent(out) :: facemat(-1:imax+1,-1:jmax+1,8),normmat(-1:imax+1,-1:jmax+1,4,2)
    real(kind=8) :: angle,const1,Mnew,angle1,angle2,i1,i2,i3,i4,height,const2,&
                    recipn,recipw,recips,recipe,mag

    const1 = 1.7857142857_8
    const2 = 0.7142857143_8
    

    ! Interior Initial
    do i = -1,(imax+1)
        do j = -1,(jmax+1)
            facemat(i,j,1) = sqrt((xmatg(i,j+1) - xmatg(i+1,j+1))**2 + (ymatg(i,j+1) - ymatg(i+1,j+1))**2)
            facemat(i,j,2) = sqrt((xmatg(i,j) - xmatg(i,j+1))**2 + (ymatg(i,j) - ymatg(i,j+1))**2)
            facemat(i,j,3) = sqrt((xmatg(i+1,j) - xmatg(i,j))**2 + (ymatg(i+1,j) - ymatg(i,j))**2)
            facemat(i,j,4) = sqrt((xmatg(i+1,j+1) - xmatg(i+1,j))**2 + (ymatg(i+1,j+1) - ymatg(i+1,j))**2)

            normmat(i,j,1,1) = (ymatg(i,j+1) - ymatg(i+1,j+1))/facemat(i,j,1)
            normmat(i,j,1,2) = -1*(xmatg(i,j+1) - xmatg(i+1,j+1))/facemat(i,j,1)
            normmat(i,j,2,1) = (ymatg(i,j) - ymatg(i,j+1))/facemat(i,j,2)
            normmat(i,j,2,2) = -1*(xmatg(i,j) - xmatg(i,j+1))/facemat(i,j,2)
            normmat(i,j,3,1) = (ymatg(i+1,j) - ymatg(i,j))/facemat(i,j,3)
            normmat(i,j,3,2) = -1*(xmatg(i+1,j) - xmatg(i,j))/facemat(i,j,3)
            normmat(i,j,4,1) = (ymatg(i+1,j+1) - ymatg(i+1,j))/facemat(i,j,4)
            normmat(i,j,4,2) = -1*(xmatg(i+1,j+1) - xmatg(i+1,j))/facemat(i,j,4)
            
            facemat(i,j,5) = (ymatg(i,j+1) - ymatg(i+1,j+1))/(xmatg(i,j+1) - xmatg(i+1,j+1))
            recipn = -1.0/facemat(i,j,5)
            facemat(i,j,6) = (ymatg(i,j) - ymatg(i,j+1))/(xmatg(i,j) - xmatg(i,j+1))
            recipw = -1.0/facemat(i,j,6)
            facemat(i,j,7) = (ymatg(i+1,j) - ymatg(i,j))/(xmatg(i+1,j) - xmatg(i,j))
            recips = -1.0/facemat(i,j,7)
            facemat(i,j,8) = (ymatg(i+1,j+1) - ymatg(i+1,j))/(xmatg(i+1,j+1) - xmatg(i+1,j))
            recipe = -1.0/facemat(i,j,8)

            angle1 = atan2(ymatg(i+1,j)-ymatg(i,j),xmatg(i+1,j)-xmatg(i,j))
            if (j==1) then
                wall_ang(i,1) = atan2(ymatg(i+1,j)-ymatg(i,j),xmatg(i+1,j)-xmatg(i,j))
            elseif (j==(jmax-1)) then
                wall_ang(i,2) = atan2(ymatg(i+1,j+1)-ymatg(i,j+1),xmatg(i+1,j+1)-xmatg(i,j+1))
            endif
            angle2 = atan2(ymatg(i+1,j+1)-ymatg(i,j+1),xmatg(i+1,j+1)-xmatg(i,j+1))
            angle = 0.5_8*(angle1+angle2)
            alfmat(i,j) = angle
            height = 1/(0.5_8*((ymatg(i,jmax)-ymatg(i,1))+(ymatg(i+1,jmax)-ymatg(i+1,1))))
            Mnew = height*M_in

            ! qmat(i,j,1) = 1.0
            ! qmat(i,j,2) = Mnew*cos(angle)
            ! qmat(i,j,3) = Mnew*sin(angle)
            ! qmat(i,j,4) = const1 + 0.5*(Mnew**2)

            ! fmat(i,j,1) = Mnew*cos(angle)
            ! fmat(i,j,2) = (Mnew*cos(angle))**2+const2
            ! fmat(i,j,3) = (Mnew**2)*cos(angle)*sin(angle)
            ! fmat(i,j,4) = (const1 + 0.5*(Mnew**2) + const2)*Mnew*cos(angle)

            ! gmat(i,j,1) = Mnew*sin(angle)
            ! gmat(i,j,2) = (Mnew**2)*cos(angle)*sin(angle)
            ! gmat(i,j,3) = (Mnew*sin(angle))**2+const2
            ! gmat(i,j,4) = (const1 + 0.5*(Mnew**2) + const2)*Mnew*sin(angle)

            qmat(i,j,1) = 1.0_8
            qmat(i,j,2) = M_in
            qmat(i,j,3) = 0.0_8
            qmat(i,j,4) = const1 + 0.5_8*(M_in**2)

            call fg(imax,jmax,qmat,fmat,gmat,po_inf,rho_inf,a_inf)

            ! fmat(i,j,1) = M_in
            ! fmat(i,j,2) = (M_in)**2+const2
            ! fmat(i,j,3) = 0.0_8
            ! fmat(i,j,4) = (const1 + 0.5_8*(M_in**2) + const2)*M_in

            ! gmat(i,j,1) = 0.0_8
            ! gmat(i,j,2) = 0.0_8
            ! gmat(i,j,3) = const2
            ! gmat(i,j,4) = 0.0_8

            ! if ((i > imax-1 .or. i<1) .and. (j > jmax-1 .or. j<1)) then
            !     qmat(i,j,1) = 0.0
            !     qmat(i,j,2) = 0.0
            !     qmat(i,j,3) = 0.0
            !     qmat(i,j,4) = 0.0

            !     fmat(i,j,1) = 0.0
            !     fmat(i,j,2) = 0.0
            !     fmat(i,j,3) = 0.0
            !     fmat(i,j,4) = 0.0

            !     gmat(i,j,1) = 0.0
            !     gmat(i,j,2) = 0.0
            !     gmat(i,j,3) = 0.0
            !     gmat(i,j,4) = 0.0
            ! endif

        ! if (MOD(jmax,2) .eq. 0) then 
        !     qmat(i,jmax/2,2) = Mnew
        !     qmat(i,jmax/2,3) = 0.0
        ! endif
        enddo
    enddo
    

end subroutine initial