subroutine fg(imax,jmax,qmat,fmat,gmat,po_inf,rho_inf,a_inf)
    implicit none
    
    integer,intent(in) :: imax,jmax
    real(kind=8),intent(in) :: qmat(-1:imax+1,-1:jmax+1,4),po_inf,rho_inf,a_inf
    real(kind=8),intent(inout) :: fmat(-1:imax+1,-1:jmax+1,4),gmat(-1:imax+1,-1:jmax+1,4)
    integer :: i,j
    real(kind=8) :: p,M,u,v

    do i = 1,imax-1
        do j = 1,jmax-1
            M = sqrt((qmat(i,j,2)**2+qmat(i,j,3)**2)/(qmat(i,j,1)))
            p = ((po_inf*101325.0)*(1.0 + 0.2 * (M**2)))/(rho_inf*a_inf**2)
            u = qmat(i,j,2)/qmat(i,j,1)
            v = qmat(i,j,3)/qmat(i,j,1)

            fmat(0,j,1) = qmat(i,j,2)
            fmat(0,j,2) = qmat(i,j,2)*u + p
            fmat(0,j,3) = qmat(i,j,2)*v
            fmat(0,j,4) = (((2.5*p) + 0.5*qmat(i,j,1)*(u**2 + v**2))+p)*u

            gmat(0,j,1) = qmat(i,j,3)
            gmat(0,j,2) = qmat(i,j,2)*u
            gmat(0,j,3) = qmat(i,j,2)*v + p
            gmat(0,j,4) = (((2.5*p) + 0.5*qmat(i,j,1)*(u**2 + v**2))+p)*v
        enddo
    enddo

end subroutine fg