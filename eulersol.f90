program eulersol
    implicit none

    character(len=80) :: fname
    integer :: i, j, imax, jmax
    real(kind=8), allocatable :: xmat(:,:), ymat(:,:)
    real(kind=8) :: M_in,po_inf,rho_inf,T_inf,a_inf,p_ex,p_inf

    ! Read in the Data for the PDE Grid into the Euler Solver
    fname = 'pde_grid2.x'
    open(unit=10, file=fname, status='old')
    read(10,*) imax, jmax
    allocate(xmat(imax,jmax))
    allocate(ymat(imax,jmax))
    read(10,*) ((xmat(i,j), i=1,imax), j=1,jmax), &
        ((ymat(i,j), i=1,imax), j=1,jmax)
    close(10) 

    ! Set the desired inlet Mach #
    M_in = 0.7_8



    ! Freestream Conditions
    p_inf = 1.0_8*101325.0_8 !atm
    p_ex =  0.93946969_8 * p_inf ! must be updated with M_in
    p_ex = 1.0_8*101325.0_8
    po_inf = p_inf*(1.0_8 + 0.2_8 * (M_in**2))**(3.5_8)
    rho_inf = 1.225_8 !kg/m^3
    T_inf = 288.203_8 ! K
    a_inf = sqrt(1.4_8*T_inf*287.0_8)
    p_inf = p_inf/(rho_inf*a_inf**2)
    !print*,p_inf
    p_ex = p_ex/(rho_inf*a_inf**2)


    ! Create matrices of Area, q, flux, and dissapation
    call matrices(xmat,ymat,imax,jmax,M_in,po_inf,p_inf,p_ex,rho_inf,T_inf,a_inf)
    !print*,imax

end program eulersol