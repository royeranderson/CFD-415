subroutine RK(imax,jmax,po_inf,rho_inf,a_inf,p_inf,p_ex,T_inf,qmat,Amat,Rmat,Dmat,normmat,facemat,alfmat,fmat,gmat,wall_ang,M_in)
    implicit none

    integer,intent(in) :: imax,jmax
    real(kind=8),intent(inout) :: qmat(-1:imax+1,-1:jmax+1,4),Rmat(1:imax-1,1:jmax-1,4),&
                                    fmat(-1:imax+1,-1:jmax+1,4),gmat(-1:imax+1,-1:jmax+1,4),&
                                    Dmat(1:imax-1,1:jmax-1,4)
    real(kind=8),intent(in) :: Amat(-1:imax+1,-1:jmax+1),&
                                normmat(-1:imax+1,-1:jmax+1,4,2),&
                                facemat(-1:imax+1,-1:jmax+1,8),po_inf,rho_inf,a_inf,&
                                p_inf,p_ex,T_inf,M_in,&
                                wall_ang(-1:imax+1,2),alfmat(-1:imax+1,-1:jmax+1)
    real(kind=8) :: tmat(1:imax-1,1:jmax-1),eigen(4),denom,delt,CFL,al1,al2,al3,al4,&
                    qmatold(-1:imax+1,-1:jmax+1,4),Rmatnew(1:imax-1,1:jmax-1,4),&
                    fmatold(-1:imax+1,-1:jmax+1,4),gmatold(-1:imax+1,-1:jmax+1,4)
    integer :: i,j

    CFL = 1.0_8

    ! Calculate the Delta t Matrix

    call eulerBCs(imax,jmax,qmat,gmat,fmat,po_inf,p_inf,p_ex,rho_inf,T_inf,M_in,wall_ang,a_inf,alfmat,normmat)

    do i=1,imax-1
        do j=1,jmax-1
            call eig(i,j,imax,jmax,po_inf,rho_inf,a_inf,qmat,normmat,eigen)
            denom = (eigen(1)*facemat(i,j,1))+(eigen(2)*facemat(i,j,2))+(eigen(3)*facemat(i,j,3))+(eigen(4)*facemat(i,j,4))
            tmat(i,j) = (2.0_8*Amat(i,j))/denom
        enddo
    enddo
    delt = minval(tmat) * CFL
    ! print*,delt

    ! Begin Runge-Kutta Steps
    al1 = 0.25_8
    al2 = (1.0_8)/(3.0_8)
    al3 = 0.5_8
    al4 = 1.0_8

    ! Step 1
    call eulerBCs(imax,jmax,qmat,gmat,fmat,po_inf,p_inf,p_ex,rho_inf,T_inf,M_in,wall_ang,a_inf,alfmat,normmat)
    ! print*,'s1'
    call fg(imax,jmax,qmat,fmat,gmat,po_inf,rho_inf,a_inf)

    qmatold = qmat
    fmatold = fmat
    gmatold = gmat

    call residuals(imax,jmax,facemat,normmat,fmat,gmat,Rmat)
    call dissapation(imax,jmax,qmat,facemat,normmat,po_inf,rho_inf,a_inf,Dmat)

    ! print*,'oc'
    do i=1,imax-1
        do j=1,jmax-1
            qmat(i,j,:) = qmatold(i,j,:) - ((al1*delt)/Amat(i,j))*(Rmat(i,j,:)-Dmat(i,j,:))
            ! if (j==9) then
            !     print*,i,j
            !     print*,'R'
            !     print*,Rmat(i,j,3)
            !     print*,'term'
            !     print*,((al1*delt)/Amat(i,j))*(Rmat(i,j,3)-Dmat(i,j,3))
            ! endif
        enddo
    enddo

    ! Step 2
    call eulerBCs(imax,jmax,qmat,gmat,fmat,po_inf,p_inf,p_ex,rho_inf,T_inf,M_in,wall_ang,a_inf,alfmat,normmat)
    ! print*,'s2'
    call fg(imax,jmax,qmat,fmat,gmat,po_inf,rho_inf,a_inf)
    call residuals(imax,jmax,facemat,normmat,fmat,gmat,Rmatnew)

    do i=1,imax-1
        do j=1,jmax-1
            qmat(i,j,:) = qmatold(i,j,:) - ((al2*delt)/Amat(i,j))*(Rmatnew(i,j,:)-Dmat(i,j,:))
        enddo
    enddo

    ! Step 3
    call eulerBCs(imax,jmax,qmat,gmat,fmat,po_inf,p_inf,p_ex,rho_inf,T_inf,M_in,wall_ang,a_inf,alfmat,normmat)
    ! print*,'s3'
    call fg(imax,jmax,qmat,fmat,gmat,po_inf,rho_inf,a_inf)
    call residuals(imax,jmax,facemat,normmat,fmat,gmat,Rmatnew)

    do i=1,imax-1
        do j=1,jmax-1
            qmat(i,j,:) = qmatold(i,j,:) - ((al3*delt)/Amat(i,j))*(Rmatnew(i,j,:)-Dmat(i,j,:))
        enddo
    enddo

    ! Step 4
    call eulerBCs(imax,jmax,qmat,gmat,fmat,po_inf,p_inf,p_ex,rho_inf,T_inf,M_in,wall_ang,a_inf,alfmat,normmat)
    ! print*,'s4'
    call fg(imax,jmax,qmat,fmat,gmat,po_inf,rho_inf,a_inf)
    call residuals(imax,jmax,facemat,normmat,fmat,gmat,Rmatnew)

    do i=1,imax-1
        do j=1,jmax-1
            qmat(i,j,:) = qmatold(i,j,:) - ((al4*delt)/Amat(i,j))*(Rmatnew(i,j,:)-Dmat(i,j,:))
        enddo
    enddo
    !print*,maxval(qmat1-qmat)

    qmatold = qmat
    fmatold = fmat
    gmatold = gmat

    call eulerBCs(imax,jmax,qmat,gmat,fmat,po_inf,p_inf,p_ex,rho_inf,T_inf,M_in,wall_ang,a_inf,alfmat,normmat)
    ! print*,'post'
    call fg(imax,jmax,qmat,fmat,gmat,po_inf,rho_inf,a_inf)
    Rmat = Rmatnew

end subroutine RK