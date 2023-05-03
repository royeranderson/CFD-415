subroutine residuals(imax,jmax,facemat,normmat,fmat,gmat,Rmat)
    implicit none

    integer,intent(in) :: imax,jmax
    real(kind=8),intent(in) :: facemat(-1:imax+1,-1:jmax+1,8),normmat(-1:imax+1,-1:jmax+1,4,2),&
                                fmat(-1:imax+1,-1:jmax+1,4),gmat(-1:imax+1,-1:jmax+1,4)
    real(kind=8),intent(inout) :: Rmat(1:imax-1,1:jmax-1,4)
    real(kind=8) :: fn,gn,fw,gw,fs,gs,fe,ge,dyn,dxn,dyw,dxw,dys,dxs,dye,dxe
    integer :: i,j,k
    
    do i = 1,imax-1
        do j = 1,jmax-1
            do k = 1,4
                fn = 0.5_8*(fmat(i,j,k)+fmat(i,j+1,k))
                gn = 0.5_8*(gmat(i,j,k)+gmat(i,j+1,k))

                fw = 0.5_8*(fmat(i,j,k)+fmat(i-1,j,k))
                gw = 0.5_8*(gmat(i,j,k)+gmat(i-1,j,k))

                fs = 0.5_8*(fmat(i,j,k)+fmat(i,j-1,k))
                gs = 0.5_8*(gmat(i,j,k)+gmat(i,j-1,k))

                fe = 0.5_8*(fmat(i,j,k)+fmat(i+1,j,k))
                ge = 0.5_8*(gmat(i,j,k)+gmat(i+1,j,k))


                dxn = -1.0_8*(facemat(i,j,1)*normmat(i,j,1,2))
                dyn = (facemat(i,j,1)*normmat(i,j,1,1))
                

                dxw = -1.0_8*(facemat(i,j,2)*normmat(i,j,2,2))
                dyw = (facemat(i,j,2)*normmat(i,j,2,1))
                

                dxs = -1.0_8*(facemat(i,j,3)*normmat(i,j,3,2))
                dys = (facemat(i,j,3)*normmat(i,j,3,1))
                

                dxe = -1.0_8*(facemat(i,j,4)*normmat(i,j,4,2))
                dye = (facemat(i,j,4)*normmat(i,j,4,1))

                Rmat(i,j,k) = fn*dyn - gn*dxn+fw*dyw - gw*dxw+fs*dys - gs*dxs+fe*dye - ge*dxe
                if (j==9 .and. k==3 .and. i==19) then
                    ! print*,'inres'
                    ! print*,i,j
                    ! print*,fn*dyn
                    ! print*,dxw
                    ! print*,fs*dys
                    ! print*,gn*dxn
                    ! print*,gs*dxs
                    ! print*,fw*dyw
                    ! print*,fe*dye
                    ! print*,gw*dxw
                    ! print*,ge*dxe
                    ! print*,Rmat(i,j,k)
                endif
            enddo
        enddo
    enddo

end subroutine residuals