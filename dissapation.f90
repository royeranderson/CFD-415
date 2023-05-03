subroutine dissapation(imax,jmax,qmat,facemat,normmat,po_inf,rho_inf,a_inf,Dmat)
    implicit none

    integer,intent(in) :: imax,jmax
    real(kind=8),intent(in) :: facemat(-1:imax+1,-1:jmax+1,8),normmat(-1:imax+1,-1:jmax+1,4,2),&
                                po_inf,rho_inf,a_inf,qmat(-1:imax+1,-1:jmax+1,4)
    real(kind=8),intent(out) ::Dmat(1:imax-1,1:jmax-1,4)
    real(kind=8) :: q,qup,qdown,qleft,qright,nu2,nu4,s21,s22,s23,s24,p1,p2,p3,p4,&
                    p5,p6,s41,s42,s43,s44,e1a,e2a,e3a,e4a,e1b,e2b,e3b,e4b,e1,e2,e3,e4,&
                    l1,l2,l3,l4,D2(4),D4(4),D2l(4),D2r(4),D4l(4),D4r(4),del31(4),del32(4),del33(4),del34(4)
    real(kind=8) :: eigen(4)
    integer :: i,j

    nu2 = 0.0_8
    nu2 = 0.5_8
    nu4 = 0.01_8
    !print*,nu4

    do i = 1,imax-1
        do j = 1,jmax-1
            call eig(i,j,imax,jmax,po_inf,rho_inf,a_inf,qmat,normmat,eigen)
            e1a = eigen(4)
            call eig(i+1,j,imax,jmax,po_inf,rho_inf,a_inf,qmat,normmat,eigen)
            e1b = eigen(2)
            call eig(i,j,imax,jmax,po_inf,rho_inf,a_inf,qmat,normmat,eigen)
            e2a = eigen(2)
            call eig(i-1,j,imax,jmax,po_inf,rho_inf,a_inf,qmat,normmat,eigen)
            e2b = eigen(4)
            call eig(i,j,imax,jmax,po_inf,rho_inf,a_inf,qmat,normmat,eigen)
            e3a = eigen(1)
            call eig(i,j+1,imax,jmax,po_inf,rho_inf,a_inf,qmat,normmat,eigen)
            e3b = eigen(3)
            call eig(i,j,imax,jmax,po_inf,rho_inf,a_inf,qmat,normmat,eigen)
            e4a = eigen(3)
            call eig(i,j-1,imax,jmax,po_inf,rho_inf,a_inf,qmat,normmat,eigen)
            e4b = eigen(1)

            e1 = .5_8*(e1a+e1b)
            e2 = .5_8*(e2a+e2b)
            e3 = .5_8*(e3a+e3b)
            e4 = .5_8*(e4a+e4b)

            l1 = facemat(i,j,4)
            l2 = facemat(i,j,2)
            l3 = facemat(i,j,1)
            l4 = facemat(i,j,3)


            call pressure(i+1,j,qmat,rho_inf,a_inf,p1)
            call pressure(i,j,qmat,rho_inf,a_inf,p2)
            call pressure(i-1,j,qmat,rho_inf,a_inf,p3)
            call pressure(i+2,j,qmat,rho_inf,a_inf,p4)
            call pressure(i+1,j,qmat,rho_inf,a_inf,p5)
            call pressure(i,j,qmat,rho_inf,a_inf,p6)
            s21 = nu2*0.5_8*((abs(p1-2.0_8*p2+p3)/(p1+2.0_8*p2+p3))+(abs(p4-2.0_8*p5+p6)/(p4+2.0_8*p5+p6)))
            
            
            call pressure(i+1,j,qmat,rho_inf,a_inf,p1)
            call pressure(i,j,qmat,rho_inf,a_inf,p2)
            call pressure(i-1,j,qmat,rho_inf,a_inf,p3)
            call pressure(i,j,qmat,rho_inf,a_inf,p4)
            call pressure(i-1,j,qmat,rho_inf,a_inf,p5)
            call pressure(i-2,j,qmat,rho_inf,a_inf,p6)
            s22 = nu2*0.5_8*((abs(p1-2.0_8*p2+p3)/(p1+2.0_8*p2+p3))+(abs(p4-2.0_8*p5+p6)/(p4+2.0_8*p5+p6)))
            

            call pressure(i,j+1,qmat,rho_inf,a_inf,p1)
            call pressure(i,j,qmat,rho_inf,a_inf,p2)
            call pressure(i,j-1,qmat,rho_inf,a_inf,p3)
            call pressure(i,j+2,qmat,rho_inf,a_inf,p4)
            call pressure(i,j+1,qmat,rho_inf,a_inf,p5)
            call pressure(i,j,qmat,rho_inf,a_inf,p6)
            s23 = nu2*0.5_8*((abs(p1-2.0_8*p2+p3)/(p1+2.0_8*p2+p3))+(abs(p4-2.0_8*p5+p6)/(p4+2.0_8*p5+p6)))
            

            call pressure(i,j+1,qmat,rho_inf,a_inf,p1)
            call pressure(i,j,qmat,rho_inf,a_inf,p2)
            call pressure(i,j-1,qmat,rho_inf,a_inf,p3)
            call pressure(i,j,qmat,rho_inf,a_inf,p4)
            call pressure(i,j-1,qmat,rho_inf,a_inf,p5)
            call pressure(i,j-2,qmat,rho_inf,a_inf,p6)
            s24 = nu2*0.5_8*((abs(p1-2.0_8*p2+p3)/(p1+2.0_8*p2+p3))+(abs(p4-2.0_8*p5+p6)/(p4+2.0_8*p5+p6)))
            

            s41 = max(0.0_8,nu4-s21)
            s42 = max(0.0_8,nu4-s22)
            s43 = max(0.0_8,nu4-s23)
            s44 = max(0.0_8,nu4-s24)
            print*,s41,s42,s43,s44

            D2l(:) = (s21*l1*e1*(qmat(i+1,j,:) - qmat(i,j,:))-s22*l2*e2*(qmat(i,j,:)-qmat(i-1,j,:)))
            !print*,e1b
            D2r(:) = (s23*l3*e3*(qmat(i,j+1,:) - qmat(i,j,:))-s24*l4*e4*(qmat(i,j,:)-qmat(i,j-1,:)))
            D2(:) = D2l + D2r


            del31 = (qmat(i+2,j,:)-2.0_8*qmat(i+1,j,:)+qmat(i,j,:))-(qmat(i+1,j,:)-2.0_8*qmat(i,j,:)+qmat(i-1,j,:))
            del32 = (qmat(i+1,j,:)-2.0_8*qmat(i,j,:)+qmat(i+1,j,:))-(qmat(i,j,:)-2.0_8*qmat(i-1,j,:)+qmat(i-2,j,:))
            del33 = (qmat(i,j+2,:)-2.0_8*qmat(i,j+1,:)+qmat(i,j,:))-(qmat(i,j+1,:)-2.0_8*qmat(i,j,:)+qmat(i,j-1,:))
            del34 = (qmat(i,j+1,:)-2.0_8*qmat(i,j,:)+qmat(i,j-1,:))-(qmat(i,j,:)-2.0_8*qmat(i,j-1,:)+qmat(i,j-2,:))
            D4l(:) = (s41*l1*e1*del31-s42*l2*e2*del32)
            D4r(:) = (s43*l3*e3*del33-s44*l4*e4*del34)
            D4(:) = D4l + D4r
            ! print*,'4'
            ! print*,D4

            Dmat(i,j,:) = D2(:)-D4(:)
            ! print*,'tot'
            ! print*,Dmat(i,j,:)


        enddo
    enddo







contains
    subroutine pressure(ip,jp,qmat,rho_inf,a_inf,p)
        integer, intent(in) :: ip,jp
        real(kind=8),intent(in) :: qmat(-1:imax+1,-1:jmax+1,4),rho_inf,a_inf
        real(kind=8), intent(out) :: p
        real(kind=8) :: q,u,v
        q = sqrt((qmat(ip,jp,2)**2+qmat(ip,jp,3)**2)/(qmat(ip,jp,1)**2))
        u = qmat(ip,jp,2)/qmat(ip,jp,1)
        v = qmat(ip,jp,3)/qmat(ip,jp,1)
        p = .4*(qmat(ip,jp,4)-.5*qmat(ip,jp,1)*(u**2+v**2))
    end subroutine pressure

end subroutine dissapation