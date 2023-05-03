subroutine eulerBCs(imax,jmax,qmat,gmat,fmat,po_inf,p_inf,p_ex,rho_inf,T_inf,M_in,wall_ang,a_inf,alfmat,normmat)
    implicit none

    integer,intent(in) :: imax,jmax
    real(kind=8),intent(inout) :: qmat(-1:imax+1,-1:jmax+1,4),fmat(-1:imax+1,-1:jmax+1,4),&
                                    gmat(-1:imax+1,-1:jmax+1,4)
    real(kind=8),intent(in) :: po_inf,rho_inf,T_inf,M_in,wall_ang(-1:imax+1,2),&
                                alfmat(-1:imax+1,-1:jmax+1),a_inf,p_inf,p_ex,&
                                normmat(-1:imax+1,-1:jmax+1,4,2)
    integer :: i,j,rowmove,rowmove2
    real(kind=8) :: const1,const2,Riem1,Riem2,theta,r,phi,p2,c2,q2,q,p,c,rho,M,u,v,s,dot

    const1 = 1.7857142857_8
    const2 = 0.7142857143_8
    Riem1 = M_in +5.0_8 ! for gam = 1.4

    ! Upper Ghost Cells
    rowmove = 0
    do j = jmax,(jmax+1)
        rowmove2 = 1+(rowmove*2)
        do i = 1,(imax-1)

            ! theta = atan2(qmat(i,j-rowmove2,3), qmat(i,j-rowmove2,2)) - wall_ang(i,2)
            ! phi = 2.0_8*atan2(tan(theta/2.0_8), -1.0_8) + wall_ang(i,2)
            ! phi = wall_ang(i,2) - theta
            ! r = sqrt((qmat(i,j-rowmove2,2)**2 + qmat(i,j-rowmove2,3)**2)/qmat(i,j-rowmove2,1)**2)

            ! qmat(i,j,1) = qmat(i,j-rowmove2,1)
            ! qmat(i,j,2) = r*cos(phi)
            ! qmat(i,j,3) = r*sin(phi)
            ! u = qmat(i,j,2)/qmat(i,j,1)
            ! v = qmat(i,j,3)/qmat(i,j,1)
            ! qmat(i,j,4) = qmat(i,j-rowmove2,4)
            ! p = .4_8*(qmat(i,j,4)-.5_8*qmat(i,j,1)*(u**2+v**2))

            ! fmat(i,j,1) = r*cos(phi)*qmat(i,j-rowmove2,1)
            ! fmat(i,j,2) = qmat(i,j-rowmove2,1)*(r*cos(phi))**2+p
            ! fmat(i,j,3) = (r**2)*cos(phi)*sin(phi)*qmat(i,j-rowmove2,1)
            ! fmat(i,j,4) = (2.5_8*p + .5_8*qmat(i,j-rowmove2,1)*((r*cos(phi))**2+(r*sin(phi))**2))*r*cos(phi)

            ! gmat(i,j,1) = r*sin(phi)*qmat(i,j-rowmove2,1)
            ! gmat(i,j,2) = (r**2)*cos(phi)*sin(phi)*qmat(i,j-rowmove2,1)
            ! gmat(i,j,3) = qmat(i,j-rowmove2,1)*(r*sin(phi))**2+p
            ! gmat(i,j,4) = (2.5_8*p + .5_8*qmat(i,j-rowmove2,1)*((r*cos(phi))**2+(r*sin(phi))**2))*r*sin(phi)

            u = qmat(i,jmax-1,2)/qmat(i,jmax-1,1)
            v = qmat(i,jmax-1,3)/qmat(i,jmax-1,1)
            dot = u*normmat(i,jmax,3,2) + v*normmat(i,jmax,3,1)

            qmat(i,jmax,1) = qmat(i,jmax-1,1)
            qmat(i,jmax,2) = u
            qmat(i,jmax,3) = v-2*dot*normmat(i,jmax,3,1)
            qmat(i,jmax,4) = qmat(i,jmax-1,4)
            

            call fg(imax,jmax,qmat,fmat,gmat,po_inf,rho_inf,a_inf)
            
            u = qmat(i,jmax-2,2)/qmat(i,jmax-2,1)
            v = qmat(i,jmax-2,3)/qmat(i,jmax-2,1)
            dot = u*normmat(i,jmax,3,2) + v*normmat(i,jmax,3,1)
            ! print*,i,jmax-2
            ! print*,v
            ! print*,normmat(i,jmax,3,1)
            ! print*,normmat(i,jmax,3,2)
            qmat(i,jmax+1,1) = qmat(i,jmax-2,1)
            qmat(i,jmax+1,2) = u
            qmat(i,jmax+1,3) = v-2*dot*normmat(i,jmax,3,2)
            qmat(i,jmax+1,4) = qmat(i,jmax-2,4)
            print*,'rho v',i,jmax
            print*,dot
            print*,normmat(i,jmax,3,2)
            print*,normmat(i,jmax,3,1)
            print*,'v'
            print*,qmat(i,jmax,3)
            print*,'u'
            print*,qmat(i,jmax,2)
            print*,sqrt(qmat(i,jmax,3)**2+qmat(i,jmax,2)**2)

            call fg(imax,jmax,qmat,fmat,gmat,po_inf,rho_inf,a_inf)
            
        enddo
        rowmove = rowmove+1
    enddo

    
    ! Lower Ghost Cells
    rowmove = 1
    do j = -1,0
        rowmove2 = 1+(rowmove*2)
        do i = 1,(imax-1)

            ! theta = atan2(qmat(i,j+rowmove2,3), qmat(i,j+rowmove2,2)) - wall_ang(i,1)
            ! phi = 2.0_8*atan2(tan(theta/2.0_8), real(-1.0,kind=8)) + wall_ang(i,1)
            ! phi = wall_ang(i,1) - theta
            ! r = sqrt((qmat(i,j+rowmove2,2)**2 + qmat(i,j+rowmove2,3)**2)/qmat(i,j+rowmove2,1)**2)

            ! qmat(i,j,1) = qmat(i,j+rowmove2,1)
            ! qmat(i,j,2) = r*cos(phi)
            ! qmat(i,j,3) = r*sin(phi)
            ! u = qmat(i,j,2)/qmat(i,j,1)
            ! v = qmat(i,j,3)/qmat(i,j,1)
            ! qmat(i,j,4) = qmat(i,j+rowmove2,4)
            ! p = .4_8*(qmat(i,j,4)-.5_8*qmat(i,j,1)*(u**2+v**2))

            ! fmat(i,j,1) = r*cos(phi)*qmat(i,j+rowmove2,1)
            ! fmat(i,j,2) = qmat(i,j+rowmove2,1)*(r*cos(phi))**2+p
            ! fmat(i,j,3) = (r**2)*cos(phi)*sin(phi)*qmat(i,j+rowmove2,1)
            ! fmat(i,j,4) = (2.5_8*p + .5_8*qmat(i,j+rowmove2,1)*((r*cos(phi))**2+(r*sin(phi))**2))*r*cos(phi)

            ! gmat(i,j,1) = r*sin(phi)*qmat(i,j+rowmove2,1)
            ! gmat(i,j,2) = (r**2)*cos(phi)*sin(phi)*qmat(i,j+rowmove2,1)
            ! gmat(i,j,3) = qmat(i,j+rowmove2,1)*(r*sin(phi))**2+p
            ! gmat(i,j,4) = (2.5_8*p + .5_8*qmat(i,j+rowmove2,1)*((r*cos(phi))**2+(r*sin(phi))**2))*r*sin(phi)

            u = qmat(i,1,2)/qmat(i,1,1)
            v = qmat(i,1,3)/qmat(i,1,1)
            dot = u*normmat(i,0,1,2) + v*normmat(i,0,1,1)

            qmat(i,0,1) = qmat(i,1,1)
            qmat(i,0,2) = u
            qmat(i,0,3) = v-2*dot*normmat(i,0,1,1)
            qmat(i,0,4) = qmat(i,1,4)
            ! print*,'rho v',i,0
            ! print*,dot
            ! print*,normmat(i,0,1,2)
            ! print*,normmat(i,0,1,1)
            ! print*,qmat(i,0,3)
            ! print*,qmat(i,0,2)
            ! print*,sqrt(qmat(i,0,3)**2+qmat(i,0,2)**2)

            call fg(imax,jmax,qmat,fmat,gmat,po_inf,rho_inf,a_inf)
            
            u = qmat(i,2,2)/qmat(i,2,1)
            v = qmat(i,2,3)/qmat(i,2,1)
            dot = u*normmat(i,0,1,2) + v*normmat(i,0,1,1)

            qmat(i,-1,1) = qmat(i,2,1)
            qmat(i,-1,2) = u
            qmat(i,-1,3) = v-2*dot*normmat(i,0,1,1)
            qmat(i,-1,4) = qmat(i,2,4)

            call fg(imax,jmax,qmat,fmat,gmat,po_inf,rho_inf,a_inf)
            
            
        enddo
        rowmove = rowmove-1
    enddo

    ! Inlet BC's 
    do j = -1,jmax+1
        qmat(-1,j,1) = 1.0_8
        qmat(-1,j,2) = M_in
        qmat(-1,j,3) = 0.0_8
        qmat(-1,j,4) = const1 + 0.5*(M_in**2)

        fmat(-1,j,1) = M_in
        fmat(-1,j,2) = (M_in)**2+const2
        fmat(-1,j,3) = 0.0_8
        fmat(-1,j,4) = (const1 + 0.5*(M_in**2) + const2)*M_in

        gmat(-1,j,1) = 0.0_8
        gmat(-1,j,2) = 0.0_8
        gmat(-1,j,3) = const2
        gmat(-1,j,4) = 0.0_8

        ! qmat(0,j,:) = qmat(-1,j,:)
        ! fmat(0,j,:) = fmat(-1,j,:)
        ! gmat(0,j,:) = gmat(-1,j,:)

        q2 = sqrt((qmat(2,j,2)**2+qmat(2,j,3)**2)/(qmat(2,j,1)**2))
        u = qmat(2,j,2)/qmat(2,j,1)
        v = qmat(2,j,3)/qmat(2,j,1)
        p2 = .4_8*(qmat(2,j,4)-.5_8*qmat(2,j,1)*(u**2+v**2))
        c2 = sqrt((1.4_8*p2)/qmat(2,j,1))
        Riem2 = q2 - (5.0_8*c2)
        q = (Riem1+Riem2)*0.5_8
        c = (Riem1-Riem2)*0.1_8
        M = q/c
        p = ((po_inf)/((1.0_8 + 0.2_8 * (M**2)))**3.5_8)/(rho_inf*a_inf**2)
        rho = (1.4_8*p)/(c**2)
        u = M*cos(alfmat(0,j))
        v = M*sin(alfmat(0,j))
        
        qmat(0,j,1) = rho
        qmat(0,j,2) = rho*u
        qmat(0,j,3) = rho*v
        qmat(0,j,4) = (2.5_8*p) + 0.5_8*rho*(u**2 + v**2)

        fmat(0,j,1) = rho*u
        fmat(0,j,2) = rho*u**2 + p
        fmat(0,j,3) = rho*u*v
        fmat(0,j,4) = (((2.5_8*p) + 0.5_8*rho*(u**2 + v**2))+p)*u

        gmat(0,j,1) = rho*v
        gmat(0,j,2) = rho*u*v
        gmat(0,j,3) = rho*v**2 + p
        gmat(0,j,4) = (((2.5_8*p) + 0.5_8*rho*(u**2 + v**2))+p)*v

    enddo

    ! Outlet BC's
    do j = -1,jmax+1
        q2 = sqrt((qmat(imax-2,j,2)**2+qmat(imax-2,j,3)**2)/(qmat(imax-2,j,1)**2))
        u = qmat(imax-2,j,2)/qmat(imax-2,j,1)
        v = qmat(imax-2,j,3)/qmat(imax-2,j,1)
        p2 = .4_8*(qmat(imax-2,j,4)-.5_8*qmat(imax-2,j,1)*(u**2+v**2))
        s = p2/(qmat(imax-2,j,1)**1.4_8)
        c2 = sqrt((1.4_8*p2)/qmat(imax-2,j,1))
        Riem1 = q2 + (5.0_8*c2)
        ! print*,'out',imax-2,j
        ! print*,v
        ! print*,qmat(imax-2,j,3)

        if((q2/c2)<1.0) then
            p = p_ex
        else
            p = p2
        endif
        !print*,p

        rho = (p/s)**(1.0_8/1.4_8)
        !print*,rho
        c = sqrt((1.4_8*p)/rho)
        q = Riem1 - (5.0_8*c)
        v = qmat(imax-2,j,3)/qmat(imax-2,j,1)
        u = sqrt(q**2 - v**2)
        
        qmat(imax-1,j,1) = rho
        qmat(imax-1,j,2) = rho*u
        
        ! print*,'out',imax-1,j
        ! print*,rho*u
        ! print*,rho
        ! print*,u
        qmat(imax-1,j,3) = rho*v
        qmat(imax-1,j,4) = (2.5_8*p) + 0.5_8*rho*(u**2 + v**2)

        fmat(imax-1,j,1) = rho*u
        fmat(imax-1,j,2) = rho*u**2 + p
        fmat(imax-1,j,3) = rho*u*v
        fmat(imax-1,j,4) = (((2.5_8*p) + 0.5_8*rho*(u**2 + v**2))+p)*u

        gmat(imax-1,j,1) = rho*v
        gmat(imax-1,j,2) = rho*u*v
        gmat(imax-1,j,3) = rho*v**2 + p
        gmat(imax-1,j,4) = (((2.5_8*p) + 0.5_8*rho*(u**2 + v**2))+p)*v

        qmat(imax,j,1:3) =2*qmat(imax-1,j,1:3)-qmat(imax-2,j,1:3)
        u = qmat(imax,j,2)/qmat(imax,j,1)
        v = qmat(imax,j,3)/qmat(imax,j,1)
        !p = .4_8*(qmat(imax,j,4)-.5_8*qmat(imax,j,1)*(u**2+v**2))
        qmat(imax,j,4) = (2.5_8*p) + 0.5_8*qmat(imax,j,1)*(u**2 + v**2)

        fmat(imax,j,1) = qmat(imax,j,2)
        fmat(imax,j,2) = qmat(imax,j,2)*u + p
        fmat(imax,j,3) = qmat(imax,j,2)*v
        fmat(imax,j,4) = (((2.5_8*p) + 0.5_8*qmat(imax,j,1)*(u**2 + v**2))+p)*u

        gmat(imax,j,1) = qmat(imax,j,3)
        gmat(imax,j,2) = qmat(imax,j,2)*u
        gmat(imax,j,3) = qmat(imax,j,2)*v + p
        gmat(imax,j,4) = (((2.5_8*p) + 0.5_8*qmat(imax,j,1)*(u**2 + v**2))+p)*v

        qmat(imax+1,j,:) = qmat(imax,j,:)
        fmat(imax+1,j,:) = fmat(imax,j,:)
        gmat(imax+1,j,:) = gmat(imax,j,:) 

    enddo  

end subroutine eulerBCs