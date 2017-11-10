program sbcalib

  use globals
  use gdal
  use mapio
  use slezanbear
  use geo

  implicit none

  type dataentry
    real(dp) :: lat, lng
    real(dp) :: sky
  end type

  real(dp) :: llatsrc, llngsrc, ulatsrc, ulngsrc
  type(source), dimension(:,:), allocatable :: src

  type(dataentry), dimension(1024), target :: dat_all
  type(dataentry), dimension(:), pointer :: dat
  integer :: ndata, i, j

  character(128), dimension(9) :: fndem

  integer, parameter :: npar = 5
  integer, parameter :: nprobes = npar + 1
  real(dp), dimension(npar), parameter :: parlo &
  &   = [ 3e-6,  0.0, 2e3, 0.0, 3.0e-7]
  real(dp), dimension(npar), parameter :: parhi &
  &   = [ 9e-6, 12.0, 9e3, 0.4, 4.0e-7]
  real(dp), dimension(npar, nprobes) :: par
  real(dp), dimension(nprobes) :: modelerr
  real(dp), dimension(:,:), allocatable :: modelsky

  real(sp), dimension(:,:), allocatable :: mapi, maph
  type(geotransform) :: gti, gth
  integer :: errno

  character(*), parameter :: fni = 'SVDNB_npp_20150101-20151231_75N060W_vcm-orm-ntl_v10_c201701311200.avg_rade9.tif'
  character(*), parameter :: datafmt = '(F9.4,1X,F9.4,3X,F5.2,1X,F4.2,3X,F5.2)'

  terrain_attenuation = .false.

  ndata = 0
  open(33, file = 'data.txt', action = 'read')
  do i = 1,size(dat_all)
    read(33, *, iostat = errno) dat_all(i) % lat, dat_all(i) % lng, dat_all(i) % sky
    if (errno .ne. 0) exit
    ndata = ndata + 1
  end do
  close(33)

  dat => dat_all(1:ndata)

  llat = minval(dat % lat)
  ulat = maxval(dat % lat)
  llng = minval(dat % lng)
  ulng = maxval(dat % lng)

  fndem(1) = "eudem_dem_5deg_n45e010.tif"
  fndem(2) = "eudem_dem_5deg_n45e015.tif"
  fndem(3) = "eudem_dem_5deg_n45e020.tif"
  fndem(4) = "eudem_dem_5deg_n50e010.tif"
  fndem(5) = "eudem_dem_5deg_n50e015.tif"
  fndem(6) = "eudem_dem_5deg_n50e020.tif"
  fndem(7) = "eudem_dem_5deg_n55e010.tif"
  fndem(8) = "eudem_dem_5deg_n55e015.tif"
  fndem(9) = "eudem_dem_5deg_n55e020.tif"

  print '("dmax = ",F0.2," km")', dmax / 1e3
  call estimate_data_limits(llat, llng, ulat, ulng, dmax * 1.05, &
  &       llatsrc, llngsrc, ulatsrc, ulngsrc)
  print '("data limits = ",F6.2,1X,F7.2)', llatsrc, llngsrc, ulatsrc, ulngsrc

  call GDALAllRegister

  call map_read_section(fni, llatsrc, llngsrc, ulatsrc, ulngsrc, mapi, gti)
  write (0, '("satellite raster size = ", 2I6)') size(mapi,1), size(mapi,2)

  read_dem: block
    integer nox, noy

    call map_gen_size(llatsrc, llngsrc, ulatsrc, ulngsrc, elev_grid_meters, &
    &              nox, noy, gth)

    allocate( maph(nox,noy) )
    write (*,'("dem raster size = ",2I6)') size(maph,1), size(maph,2)

    call map_read_tiles(fndem, maph, gth)
  end block read_dem

  call random_seed()
  allocate(modelsky(ndata, nprobes))

  allocate( src(size(mapi,1), size(mapi,2)) )
  call mksources(mapi, gti, maph, gth, src)

  call random_number(par)

  do i = 1, npar
    par(i,:) = parlo(i) + par(i,:) * (parhi(i) - parlo(i))
  end do

  call simplex(par, parlo, parhi, modelerr, getdeviation)

  write (0,*) '--------------------'
  do i = 1, nprobes
    write (*, '(F8.3,F7.3,F8.1,F6.3,F7.3,Es11.3)') par(1,i) * 1e6, &
    & par(2,i), par(3,i), par(4,i), par(5,i) * 1e7, modelerr(i)
  end do

  deallocate(src)

  deallocate(mapi, maph)

contains

  pure real(dp) function getdeviation(par) result(dev)
    real(dp), dimension(:), intent(in) :: par
    type(modelpar) :: mpar
    real(dp), dimension(ndata) :: modelsky
    real(dp) :: sky(5), hobs
    integer :: i

    mpar % alpha  = par(1)
    mpar % relabs = par(2)
    mpar % hscale = par(3)
    mpar % beta   = par(4)
    mpar % skybg  = par(5)

    iterate_datpoints: do i = 1, ndata
      call onepoint(mpar, dat(i) % lat, dat(i) % lng, &
      &       src, maph, gth, sky, hobs)
      modelsky(i) = sky(4)
    end do iterate_datpoints

    ! dev = sqrt(sum((ity2mag(modelsky) - dat % sky)**2) / ndata)
    dev = sum((modelsky - mag2ity(dat % sky))**2 / modelsky)
  end function

  subroutine simplex(x,xlo,xhi,y,f)
    interface
      pure real(dp) function simplexfun(x)
        import sp,dp
        real(dp), dimension(:), intent(in) :: x
      end function
    end interface
    real(dp), dimension(:,:), intent(inout) :: x
    real(dp), dimension(:), intent(in) :: xlo, xhi
    real(dp), dimension(:), intent(out) :: y
    procedure(simplexfun) :: f

    real(dp), dimension(*), parameter :: a = [ 1.0, 2.0, 0.5, 0.5 ]
    integer,  dimension(size(x,1) + 1) :: s
    real(dp), dimension(size(x,1)) :: x_0, x_r, x_e, x_c
    real(dp) :: y_r
    integer :: iter, i, n

    character(*), parameter :: actionfmt = '("simplex ", I2, ": ", &
          & A12, "(", I2, ")")'
    character(*), parameter :: statfmt = &
          & '(8X, "miny = ", Es10.3, " (", I2, ")", &
          &   3X, "maxy = ", Es10.3, " (", I2, ")")'

    n = size(x,1)

    if (size(x,2) /= n + 1) error stop "simplex: x must have shape (n,n+1)"
    if (size(y)   /= n + 1) error stop "simplex: size(y) /= size(x,1) + 1"
    if (size(xlo) /= n) error stop "simplex: size(xlo) /= size(x,1)"
    if (size(xhi) /= n) error stop "simplex: size(xhi) /= size(x,1)"

    !$omp parallel do
    do i = 1,n+1
      y(i) = f(x(:,i))
    end do
    !$omp end parallel do

    simplex_loop: do iter = 1, 72

      call isort(y,s)
      x_0 = sum(x(:,s(1:n)),2) / n

      write (0, statfmt) y(s(1)), s(1), y(s(n+1)), s(n+1)

      ! reflection
      x_r = x_0 + a(1) * (x_0 - x(:,s(n+1)))
      y_r = f(x_r)
      if (y_r < y(s(n))) then
        if (y_r > y(s(1))) then
          write (0,actionfmt) iter, "REFLECT", s(n+1)
          x(:,s(n+1)) = x_r
        else
          ! expansion
          x_e = x_0 + a(2) * (x_r - x_0)
          if (f(x_e) < y_r) then
            write (0,actionfmt) iter, "EXPAND", s(n+1)
            x(:,s(n+1)) = x_e
          else
            write (0,actionfmt) iter, "REFLECT", s(n+1)
            x(:,s(n+1)) = x_r
          end if
        end if
        y(s(n+1)) = f(x(:,s(n+1)))
      else
        ! contraction
        x_c = x_0 + a(3) * (x(:,s(n+1)) - x_0)
        if (f(x_c) < y(s(n+1))) then
          write (0,actionfmt) iter, "CONTRACT", s(n+1)
          x(:,s(n+1)) = x_c
          y(s(n+1)) = f(x(:,s(n+1)))
        else
          ! shrink
          write (0,actionfmt) iter, "SHRINK", s(1)
          !$omp parallel do
          do i = 2,n+1
            x(:,s(i)) = x(:,s(1)) + a(4) * (x(:,s(i)) - x(:,s(1)))
            y(s(i)) = f(x(:,s(i)))
          end do
          !$omp end parallel do
          ! do concurrent (i = 2:n+1)
          !   x(:,s(i)) = x(:,s(1)) + a(4) * (x(:,s(i)) - x(:,s(1)))
          !   y(s(i)) = f(x(:,s(i)))
          ! end do
        end if
      end if
    end do simplex_loop

    call isort(y,s)
    write (0, statfmt) y(s(1)), s(1), y(s(n+1)), s(n+1)

    x = x(:,s)
    y = y(s)

  end subroutine

  pure subroutine isort(x,s)
    real(dp), intent(in)  :: x(:)
    integer,  intent(out) :: s(:)
    integer :: i,j,sref

    if (size(x) /= size(s)) error stop "isort: size(x) /= size(s)"

    forall (i = 1:size(x)) s(i) = i

    do i = 2, size(x)
      sref = s(i)
      j = i - 1
      do while ( j > 0 .and. x(s(j)) > x(sref) )
        s(j + 1) = s(j)
        j = j - 1
      end do
      s(j + 1) = sref
    end do
  end subroutine

end program sbcalib
