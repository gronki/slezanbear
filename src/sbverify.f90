program sbverify

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
  type(modelpar) :: mpar


  character(128), dimension(9) :: fndem

  integer, parameter :: npar = 5
  real(dp), dimension(npar), parameter :: parlo &
  &   = [ 1e-6,  0.0, 0.5e3, 0.0, 2.8e-7]
  real(dp), dimension(npar), parameter :: parhi &
  &   = [10e-6, 24.0,  18e3, 0.8, 4.0e-7]
  real(dp), dimension(npar, npar + 1) :: par
  real(dp), dimension(npar + 1) :: modelerr
  real(dp), dimension(:,:), allocatable :: modelsky

  real(sp), dimension(:,:), allocatable :: mapi, maph
  type(geotransform) :: gti, gth
  integer :: errno

  character(*), parameter :: fni = 'SVDNB_npp_20150101-20151231_75N060W_vcm-orm-ntl_v10_c201701311200.avg_rade9.tif'
  character(*), parameter :: datafmt = '(F9.4,1X,F9.4,3X,F5.2,1X,F4.2,3X,F5.2)'

  read (*,*) mpar % alpha, mpar % relabs, mpar % hscale, mpar % beta, &
  &       mpar % skybg

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
  allocate(modelsky(ndata, npar + 1))

  printcalib: block
    real(dp) :: sky(4), hobs

    allocate( src(size(mapi,1), size(mapi,2)) )
    call mksources(mapi, gti, maph, gth, src)

    iterate_datpoints: do i = 1, ndata
      call onepoint(mpar, dat(i) % lat, dat(i) % lng, &
      &       src, maph, gth, sky, hobs)
      write (*,'(2F10.4,3F10.3)') &
      &       dat(i) % lat, dat(i) % lng, dat(i) % sky, ity2mag(sky(4)), &
      &       ity2mag(sky(4)) - dat(i) % sky
    end do iterate_datpoints

    deallocate(src)
  end block printcalib

  deallocate(mapi, maph)

end program sbverify
