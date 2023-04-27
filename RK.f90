subroutine RK(imax,jmax,po_inf,rho_inf,a_inf,p_inf,p_ex,T_inf,qmat,Amat,Rmat,Dmat,normmat,facemat,alfmat,fmat,gmat,wall_ang,M_in)
    implicit none

    integer,intent(in) :: imax,jmax
    real(kind=8),intent(inout) :: qmat(-1:imax+1,-1:jmax+1,4),Rmat(1:imax-1,1:jmax-1,4),&
                                    fmat(-1:imax+1,-1:jmax+1,4),gmat(-1:imax+1,-1:jmax+1,4)
    real(kind=8),intent(in) :: Amat(-1:imax+1,-1:jmax+1),&
                                Dmat(1:imax-1,1:jmax-1,4),normmat(1:imax-1,1:jmax-1,4,2),&
                                facemat(1:imax-1,1:jmax-1,8),po_inf,rho_inf,a_inf,&
                                p_inf,p_ex,T_inf,M_in,&
                                wall_ang(1:imax-1,2),alfmat(-1:imax+1,-1:jmax+1)
    real(kind=8) :: qmat1(-1:imax+1,-1:jmax+1,4),qmat2(-1:imax+1,-1:jmax+1,4),&
                    qmat3(-1:imax+1,-1:jmax+1,4),qmat4(-1:imax+1,-1:jmax+1,4),&
                    fmat1(-1:imax+1,-1:jmax+1,4),&
                    gmat1(-1:imax+1,-1:jmax+1,4),&
                    tmat(1:imax-1,1:jmax-1),eigen(4),denom,delt,CFL,al1,al2,al3,al4,&
                    qmatold(-1:imax+1,-1:jmax+1,4),Rmatnew(1:imax-1,1:jmax-1,4)
    integer :: i,j

    CFL = 1.0

    ! Calculate the Delta t Matrix
    do i=1,imax-1
        do j=1,jmax-1
            call eig(i,j,imax,jmax,po_inf,rho_inf,a_inf,qmat,normmat,eigen)
            denom = (eigen(1)*facemat(i,j,1))+(eigen(2)*facemat(i,j,2))+(eigen(3)*facemat(i,j,3))+(eigen(4)*facemat(i,j,4))
            tmat(i,j) = (2.0*Amat(i,j))/denom
        enddo
    enddo
    delt = minval(tmat) * CFL

    ! Begin Runge-Kutta Steps

    ! Step 1
    fmat1 = fmat
    gmat1 = gmat
    call eulerBCs(imax,jmax,qmat,gmat1,fmat1,po_inf,p_inf,p_ex,rho_inf,T_inf,M_in,wall_ang,a_inf,alfmat)
    call fg(imax,jmax,qmat,fmat1,gmat1,po_inf,rho_inf,a_inf)
    call residuals(imax,jmax,facemat,normmat,fmat1,gmat1,Rmat)

    do i=1,imax-1
        do j=1,jmax-1
            qmat1(i,j,:) = qmat(i,j,:) - ((al1*delt)/Amat(i,j))*(Rmat(i,j,:))
        enddo
    enddo

    ! Step 2
    call fg(imax,jmax,qmat1,fmat1,gmat1,po_inf,rho_inf,a_inf)
    call eulerBCs(imax,jmax,qmat1,gmat1,fmat1,po_inf,p_inf,p_ex,rho_inf,T_inf,M_in,wall_ang,a_inf,alfmat)
    call residuals(imax,jmax,facemat,normmat,fmat1,gmat1,Rmatnew)

    do i=1,imax-1
        do j=1,jmax-1
            qmat2(i,j,:) = qmat(i,j,:) - ((al2*delt)/Amat(i,j))*(Rmatnew(i,j,:))
        enddo
    enddo

    ! Step 3
    call fg(imax,jmax,qmat2,fmat1,gmat1,po_inf,rho_inf,a_inf)
    call eulerBCs(imax,jmax,qmat2,gmat1,fmat1,po_inf,p_inf,p_ex,rho_inf,T_inf,M_in,wall_ang,a_inf,alfmat)
    call residuals(imax,jmax,facemat,normmat,fmat1,gmat1,Rmatnew)

    do i=1,imax-1
        do j=1,jmax-1
            qmat3(i,j,:) = qmat(i,j,:) - ((al3*delt)/Amat(i,j))*(Rmatnew(i,j,:))
        enddo
    enddo

    ! Step 4
    call fg(imax,jmax,qmat3,fmat1,gmat1,po_inf,rho_inf,a_inf)
    call eulerBCs(imax,jmax,qmat3,gmat1,fmat1,po_inf,p_inf,p_ex,rho_inf,T_inf,M_in,wall_ang,a_inf,alfmat)
    call residuals(imax,jmax,facemat,normmat,fmat1,gmat1,Rmatnew)

    do i=1,imax-1
        do j=1,jmax-1
            qmat4(i,j,:) = qmat(i,j,:) - ((al4*delt)/Amat(i,j))*(Rmatnew(i,j,:))
        enddo
    enddo
    !print*,maxval(qmat1-qmat)
    qmatold = qmat
    qmat = qmat4
    call fg(imax,jmax,qmat4,fmat1,gmat1,po_inf,rho_inf,a_inf)
    fmat = fmat1
    gmat = gmat1
    call eulerBCs(imax,jmax,qmat,gmat,fmat,po_inf,p_inf,p_ex,rho_inf,T_inf,M_in,wall_ang,a_inf,alfmat)
    Rmat = Rmatnew

end subroutine RK