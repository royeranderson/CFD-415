program eulersol
    implicit none

    character(len=80) :: fname
    integer :: i, j, imax, jmax
    real(kind=8), allocatable :: xmat(:,:), ymat(:,:)
    real(kind=8) :: M_in,po_inf,rho_inf,T_inf

    ! Read in the Data for the PDE Grid into the Euler Solver
    fname = 'pde_grid.x'
    open(unit=10, file=fname, status='old')
    read(10,*) imax, jmax
    allocate(xmat(imax,jmax))
    allocate(ymat(imax,jmax))
    read(10,*) ((xmat(i,j), i=1,imax), j=1,jmax), &
        ((ymat(i,j), i=1,imax), j=1,jmax)
    close(10) 

    ! Set the desired inlet Mach #
    M_in = 0.3

    ! Freestream Conditions
    po_inf = 1.0 !atm
    rho_inf = 1.225 !kg/m^3
    T_inf = 298 ! K


    ! Create matrices of Area, q, flux, and dissapation
    call matrices(xmat,ymat,imax,jmax,M_in,po_inf,rho_inf,T_inf)

end program eulersol