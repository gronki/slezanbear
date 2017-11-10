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
  type(dataentry), dimension(1024), target :: dat_all
  type(dataentry), dimension(:), pointer :: dat
  integer :: ndata, i, j

  character(128), dimension(9) :: fndem

  real(sp), dimension(:,:), allocatable :: mapi, maph
  type(geotransform) :: gti, gth
  integer :: errno
  real(dp) :: sky(4), hobs

  character(*), parameter :: fni = 'SVDNB_npp_20150101-20151231_75N060W_vcm-orm-ntl_v10_c201701311200.avg_rade9.tif'
  character(*), parameter :: datafmt = '(F9.4,1X,F9.4,3X,F5.2,1X,F4.2,3X,F5.2)'

  terrain_attenuation = .false.

  ndata = 0
  open(33, file = 'data.txt', action = 'read')
  do i = 1,size(dat_all)
    read(33, *, iostat = errno) dat_all(i) % lat, dat_all(i) % lng, &
      & dat_all(i) % sky, dat_all(i) % err
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

  write (0,'("dmax = ",F0.2," km")') dmax / 1e3
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
    write (0,'("dem raster size = ",2I6)') size(maph,1), size(maph,2)

    call map_read_tiles(fndem, maph, gth)
  end block read_dem

  call random_seed()

  printcalib: block
    type(source), dimension(:,:), allocatable :: src
    real :: a(5)
    real(dp) :: err(ndata)

    allocate( src(size(mapi,1), size(mapi,2)) )
    call mksources(mapi, gti, maph, gth, src)

    open(34, file = 'calib.txt', action = 'write')
    write (34,'(A8,A7,A8,A6,A7,A10,A10,A7)') 'SCTR', 'ABS', 'HSCAL', 'SIDE', &
    & 'BG', 'CHISQ', 'DEVITY', 'DEVMAG'

    iterate_parameters: do j = 1,16807
      call random_number(a)

      par % alpha =  M(a(1), 3.0, 24.0) * 1e-6
      par % relabs = M(a(2), 0.00, 12.00)
      par % hscale = M(a(3), 2.0, 12.0) * 1e3
      par % beta =   M(a(4), 0.0, 0.4)
      par % skybg =  mag2ity(M(a(5), 22.1, 21.6))

      !$omp parallel do private(sky,hobs)
      iterate_datpoints: do i = 1, ndata
        call onepoint(par, dat(i) % lat, dat(i) % lng, &
        &       src, maph, gth, sky, hobs)
        dat(i) % model = sky(4)
      end do iterate_datpoints
      !$omp end parallel do

      write (34,'(F8.3,F7.3,F8.1,F6.3,F7.3,ES10.3,ES10.3,F7.4)') &
        & par % alpha * 1e6, par % relabs, &
        & par % hscale, par % beta, (par % skybg) * 1e7, &
        &      sum((mag2ity(dat % sky) - dat % model)**2  / dat % model), &
        & sqrt(sum((mag2ity(dat % sky) - dat % model)**2) / ndata), &
        & sqrt(sum((dat % sky - ity2mag(dat % model))**2) / ndata)
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
