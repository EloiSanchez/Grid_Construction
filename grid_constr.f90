program grid_construction
implicit none

    integer :: nx, ny, nz, n4
    real(kind = 8) :: hx, hy, hz
    character(len=30) :: fileout

    complex(kind=8), allocatable :: psi(:,:,:)

    real(kind = 8) :: x, y, z, xmin, xmax, ymin, ymax, zmin, zmax
    real(kind = 8) :: a, b
    integer :: i, j, k

    namelist /input/ nx, ny, nz, hx, hy, hz, n4, fileout

    open(unit=10, file="input.txt", status="old")
    print *, "opened"
    read(unit=10,nml=input)
    print *, "read"
    close(10)
    print *, "closed"

    xmin = - hx * dble(nx - 1) * 0.5d0
    ymin = - hy * dble(ny - 1) * 0.5d0
    zmin = - hz * dble(nz - 1) * 0.5d0

    xmax = - xmin
    ymax = - ymin
    zmax = - zmin

    print *, ""
    print *, "PARAMETERS=================================="
    print *, "nx  ", nx, "ny  ", ny, "nz  ", nz
    print *, "hx  ", hx, "hy  ", hy, "hz  ", hz
    print *, "xmin", xmin, "ymin", ymin, "zmin", zmin
    print *, "xmax", xmax, "ymax", ymax, "zmax", zmax
    print *, "============================================"
    print *, ""

    allocate(psi(nx,ny,nz))
    
    a = 5.d-3
    b = 1.d-3

    do i = 1, nx
        x = xmin + dble(i - 1) * hx 
        do j = 1, ny
            y = ymin + dble(j - 1) * hy 
            do k = 1, nz
                z = zmin + dble(k - 1) * hz 
                psi(i,j,k) = (exp(-a*(x**2 + y**2 + z**2)) - exp(b*(x**2 + y**2 + z**2))) * cmplx(1,1)
                ! if (j == ny/2 .and. k == nz/2) print *, x, y, z, psi(i,j,k)
            end do
        end do
    end do
    
    psi = psi * sqrt(dble(n4)/(sum(psi*conjg(psi))*hx*hy*hz))

    print *, "Final number of He atoms is", real(sum(psi*conjg(psi)))*hx*hy*hz

    open(10,file=fileout,form='UNFORMATTED', status="new")
    write(10) 0.d0,xmax,ymax,zmax,hx,hy,hz,nx,ny,nz,.true.
    write(10) psi
    close(10)
end program