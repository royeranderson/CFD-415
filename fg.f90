subroutine fg(imax,jmax,qmat,fmat,gmat,po_inf,rho_inf,a_inf)
    implicit none
    
    integer,intent(in) :: imax,jmax
    real(kind=8),intent(in) :: qmat(-1:imax+1,-1:jmax+1,4),po_inf,rho_inf,a_inf
    real(kind=8),intent(inout) :: fmat(-1:imax+1,-1:jmax+1,4),gmat(-1:imax+1,-1:jmax+1,4)
    integer :: i,j
    real(kind=8) :: p,M,u,v

    do i = -1,imax+1
        do j = -1,jmax+1
            u = qmat(i,j,2)/qmat(i,j,1)
            v = qmat(i,j,3)/qmat(i,j,1)
            p = .4_8*(qmat(i,j,4)-.5_8*qmat(i,j,1)*(u**2+v**2))

            if (j==8) then
                !print*,v
            endif


            fmat(i,j,1) = qmat(i,j,2)
            fmat(i,j,2) = qmat(i,j,2)*u + p
            fmat(i,j,3) = qmat(i,j,2)*v
            fmat(i,j,4) = (((2.5_8*p) + 0.5_8*qmat(i,j,1)*(u**2 + v**2))+p)*u

            gmat(i,j,1) = qmat(i,j,3)
            gmat(i,j,2) = qmat(i,j,2)*u
            gmat(i,j,3) = qmat(i,j,2)*v + p
            gmat(i,j,4) = (((2.5_8*p) + 0.5_8*qmat(i,j,1)*(u**2 + v**2))+p)*v
        enddo
    enddo

end subroutine fg