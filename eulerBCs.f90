subroutine eulerBCs(imax,jmax,qmat,gmat,fmat,po_inf,p_inf,p_ex,rho_inf,T_inf,M_in,wall_ang,a_inf,alfmat)
    implicit none

    integer,intent(in) :: imax,jmax
    real(kind=8),intent(inout) :: qmat(-1:imax+1,-1:jmax+1,4),fmat(-1:imax+1,-1:jmax+1,4),&
                                    gmat(-1:imax+1,-1:jmax+1,4)
    real(kind=8),intent(in) :: po_inf,rho_inf,T_inf,M_in,wall_ang(1:imax-1,2),&
                                alfmat(-1:imax+1,-1:jmax+1),a_inf,p_inf,p_ex
    integer :: i,j,rowmove,rowmove2
    real(kind=8) :: const1,const2,Riem1,Riem2,theta,r,phi,p2,c2,q2,q,p,c,rho,M,u,v,s

    const1 = 1.7857142857
    const2 = 0.7142857143
    Riem1 = M_in +5.0 ! for gam = 1.4

    ! Upper Ghost Cells
    rowmove = 0
    do j = jmax,(jmax+1)
        rowmove2 = 1+(rowmove*2)
        do i = 1,(imax-1)

            theta = atan2(qmat(i,j-rowmove2,3), qmat(i,j-rowmove2,2)) - wall_ang(i,2)
            phi = 2.0*atan2(tan(theta/2.0), real(-1.0,kind=8))
            r = sqrt(qmat(i,j-rowmove2,2)**2 + qmat(i,j-rowmove2,3)**2)

            qmat(i,j,1) = qmat(i,j-rowmove2,1)
            qmat(i,j,2) = r*cos(phi)
            qmat(i,j,3) = r*sin(phi)
            qmat(i,j,4) = qmat(i,j-rowmove2,4)

            fmat(i,j,1) = r*cos(phi)
            fmat(i,j,2) = (r*cos(phi))**2+const2
            fmat(i,j,3) = (r**2)*cos(phi)*sin(phi)
            fmat(i,j,4) = (const1 + 0.5*(r**2) + const2)*r*cos(phi)

            gmat(i,j,1) = r*sin(phi)
            gmat(i,j,2) = (r**2)*cos(phi)*sin(phi)
            gmat(i,j,3) = (r*sin(phi))**2+const2
            gmat(i,j,4) = (const1 + 0.5*(r**2) + const2)*r*sin(phi)
            
        enddo
        rowmove = rowmove+1
    enddo

    
    ! Lower Ghost Cells
    rowmove = 0
    do j = -1,0
        rowmove2 = 1+(rowmove*2)
        do i = 1,(imax-1)

            theta = atan2(qmat(i,j+rowmove2,3), qmat(i,j+rowmove2,2)) - wall_ang(i,1)
            phi = 2.0*atan2(tan(theta/2.0), real(-1.0,kind=8))
            r = sqrt(qmat(i,j+rowmove2,2)**2 + qmat(i,j+rowmove2,3)**2)

            qmat(i,j,1) = qmat(i,j+rowmove2,1)
            qmat(i,j,2) = r*cos(phi)
            qmat(i,j,3) = r*sin(phi)
            qmat(i,j,4) = qmat(i,j+rowmove2,4)

            fmat(i,j,1) = r*cos(phi)
            fmat(i,j,2) = (r*cos(phi))**2+const2
            fmat(i,j,3) = (r**2)*cos(phi)*sin(phi)
            fmat(i,j,4) = (const1 + 0.5*(r**2) + const2)*r*cos(phi)

            gmat(i,j,1) = r*sin(phi)
            gmat(i,j,2) = (r**2)*cos(phi)*sin(phi)
            gmat(i,j,3) = (r*sin(phi))**2+const2
            gmat(i,j,4) = (const1 + 0.5*(r**2) + const2)*r*sin(phi)
            
        enddo
        rowmove = rowmove+1
    enddo

    ! Inlet BC's 
    do j = 1,jmax-1
        q2 = sqrt((qmat(1,j,2)**2+qmat(1,j,3)**2)/(qmat(1,j,1)))
        p2 = ((po_inf*101325.0)*(1.0 + 0.2 * (q2**2)))/(rho_inf*a_inf**2)
        c2 = sqrt((1.4*p2)/qmat(1,j,1))
        Riem2 = q2 - (5.0*c2)
        q = (Riem1+Riem2)*0.5
        c = (Riem1-Riem2)*0.1
        M = q/c
        p = ((po_inf*101325.0)*(1.0 + 0.2 * (M**2)))/(rho_inf*a_inf**2)
        rho = (1.4*p)/(c**2)
        u = M*cos(alfmat(0,j))
        v = M*sin(alfmat(0,j))
        
        qmat(0,j,1) = rho
        qmat(0,j,2) = rho*u
        qmat(0,j,3) = rho*v
        qmat(0,j,4) = (2.5*p) + 0.5*rho*(u**2 + v**2)

        fmat(0,j,1) = rho*u
        fmat(0,j,2) = rho*u**2 + p
        fmat(0,j,3) = rho*u*v
        fmat(0,j,4) = (((2.5*p) + 0.5*rho*(u**2 + v**2))+p)*u

        gmat(0,j,1) = rho*v
        gmat(0,j,2) = rho*u*v
        gmat(0,j,3) = rho*v**2 + p
        gmat(0,j,4) = (((2.5*p) + 0.5*rho*(u**2 + v**2))+p)*v

    enddo


    ! Outlet BC's
    do j = jmax,jmax+1
        q2 = sqrt((qmat(imax-1,j,2)**2+qmat(imax-1,j,3)**2)/(qmat(imax-1,j,1)))
        p2 = ((po_inf*101325.0)*(1.0 + 0.2 * (q2**2)))/(rho_inf*a_inf**2)
        s = p2/(qmat(imax-1,j,1)**1.4)
        c2 = sqrt((1.4*p2)/qmat(imax-1,j,1))
        Riem1 = q2 + (5.0*c2)

        rho = (p_ex/s)**1.4
        c = sqrt((1.4*p_ex)/rho)
        q = Riem1 - (5.0*c)
        v = qmat(imax-1,j,3)/qmat(imax-1,j,1)
        u = sqrt(q**2 - v**2)
        
        qmat(1,j,1) = rho
        qmat(1,j,2) = rho*u
        qmat(1,j,3) = rho*v
        qmat(1,j,4) = (2.5*p_ex) + 0.5*rho*(u**2 + v**2)

        fmat(1,j,1) = rho*u
        fmat(1,j,2) = rho*u**2 + p_ex
        fmat(1,j,3) = rho*u*v
        fmat(1,j,4) = (((2.5*p_ex) + 0.5*rho*(u**2 + v**2))+p_ex)*u

        gmat(1,j,1) = rho*v
        gmat(1,j,2) = rho*u*v
        gmat(1,j,3) = rho*v**2 + p_ex
        gmat(1,j,4) = (((2.5*p_ex) + 0.5*rho*(u**2 + v**2))+p_ex)*v

    enddo
    

end subroutine eulerBCs