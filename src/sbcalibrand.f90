program sbcalibrand

  use globals
  use gdal
  use mapio
  use slezanbear
  use geo

  implicit none

  type dataentry
    real(dp) :: lat, lng
    real(dp) :: sky, err
    real(dp) :: model
  end type

  real(dp) :: llatsrc, llngsrc, ulatsrc, ulngsrc
  type(modelpar) :: par
  type(dataentry), dimension(100) :: dat
  integer :: ndata, i, j

  character(128), dimension(9) :: fndem

  real(sp), dimension(:,:), allocatable :: mapi, maph
  type(geotransform) :: gti, gth
  integer :: errno
  real(dp) :: sky(4), hobs

  character(*), parameter :: fni = 'SVDNB_npp_20150101-20151231_75N060W_vcm-orm-ntl_v10_c201701311200.avg_rade9.tif'
  character(*), parameter :: datafmt = '(F9.4,1X,F9.4,3X,F5.2,1X,F4.2,3X,F5.2)'

  ndata = 0

  open(33, file = 'data.txt', action = 'read')
  do i = 1,size(dat)
    read(33, *, iostat = errno) dat(i) % lat, dat(i) % lng, &
      & dat(i) % sky, dat(i) % err
    if (errno .ne. 0) exit
    ndata = ndata + 1
  end do
  close(33)

  llat = minval(dat(1:ndata) % lat)
  ulat = maxval(dat(1:ndata) % lat)
  llng = minval(dat(1:ndata) % lng)
  ulng = maxval(dat(1:ndata) % lng)

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

  printcalib: block
    type(source), dimension(:,:), allocatable :: src
    real :: a(5)
    real(dp) :: err(100), erra, errb

    allocate( src(size(mapi,1), size(mapi,2)) )
    call mksources(mapi, gti, maph, gth, src)

    open(34, file = 'calib.txt', action = 'write')
    write (34,'(A8,A7,A8,A6,A7,3A9)') 'SCTR', 'ABS', 'HSCAL', 'SIDE', 'BG', &
    & 'DEV', 'D(21.5)', 'D(18.5)'

    iterate_parameters: do j = 1, 16000
      call random_number(a)

      par % alpha =  M(a(1), 1.0, 10.0) * 1e-6
      par % relabs = M(a(2), 0.00, 24.00)
      par % hscale = M(a(3), 0.50, 18.0) * 1e3
      par % beta =   M(a(4), 0.0, 0.8)
      par % skybg =  M(a(5), 2.8, 4.0) * 1e-7

      par % alpha = (par % alpha) / (par % hscale / 8500)**0.333

      !$omp parallel do private(sky,hobs)
      iterate_datpoints: do i = 1, ndata
        call onepoint(par, dat(i) % lat, dat(i) % lng, &
        &       src, maph, gth, sky, hobs)
        dat(i) % model = ity2mag(sky(4))
      end do iterate_datpoints
      !$omp end parallel do

      err(1:ndata) = dat(1:ndata) % model - dat(1:ndata) % sky
      call regr(dat(1:ndata) % sky, err(1:ndata), erra, errb)

      write (34,'(F8.3,F7.3,F8.1,F6.3,F7.3,3F9.4)') &
        & par % alpha * 1e6, par % relabs, &
        & par % hscale, par % beta, (par % skybg) * 1e7, &
        & sqrt(sum(err(1:ndata)**2) / ndata), &
        & errb + erra * 21.5, errb + erra * 18.5
      flush(34)
    end do iterate_parameters

    close(34)

    deallocate(src)
  end block printcalib

  deallocate(mapi, maph)

contains

  elemental real function m(a, l, u) result(b)
    real, intent(in) :: a
    real, intent(in) :: l, u
    b = l + a * (u - l)
  end function

end program sbcalibrand
