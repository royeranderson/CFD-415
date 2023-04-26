subroutine eig(ie,je,imax,jmax,po_inf,rho_inf,a_inf,qmat,normmat,eigen)
    implicit none
    integer,intent(in) :: ie,je,imax,jmax
    real(kind=8),intent(in) :: qmat(-1:imax+1,-1:jmax+1,4),normmat(1:imax-1,1:jmax-1,4,2),&
                                po_inf,rho_inf,a_inf
    real(kind=8),intent(out) :: eigen(4)
    real(kind=8) :: c,ue,un,q,p,u,v

    q = sqrt((qmat(ie,je,2)**2+qmat(ie,je,3)**2)/(qmat(ie,je,1)))
    u = qmat(ie,je,2)/qmat(ie,je,1)
    v = qmat(ie,je,3)/qmat(ie,je,1)
    p = ((po_inf*101325.0)*(1.0 + 0.2 * (q**2)))/(rho_inf*a_inf**2)
    c = sqrt((1.4*p)/qmat(ie,je,1))

    eigen(1) = (abs(u*normmat(ie,je,1,1)) + abs(v*normmat(ie,je,1,2)))+c
    eigen(2) = (abs(u*normmat(ie,je,2,2)) + abs(v*normmat(ie,je,2,1)))+c
    eigen(3) = (abs(u*normmat(ie,je,3,1)) + abs(v*normmat(ie,je,3,2)))+c
    eigen(4) = (abs(u*normmat(ie,je,4,2)) + abs(v*normmat(ie,je,4,1)))+c

end subroutine