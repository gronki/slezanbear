program sbcalib

  use globals
  use gdal
  use mapio
  use slezanbear
  use geo

  implicit none

  type dataentry
    real(fp) :: lat, lng
    real(fp) :: sky, err
    real(fp) :: model
  end type

  real(fp) :: llatsrc, llngsrc, ulatsrc, ulngsrc
  type(modelpar), dimension(:), allocatable :: par
  type(dataentry), dimension(100) :: dat
  integer :: ndata, i, j

  character(128), dimension(9) :: fndem

  real(sp), dimension(:,:), allocatable :: mapi, maph
  type(geotransform) :: gti, gth
  integer :: errno
  real(fp) :: sky(4), hobs

  character(*), parameter :: fni = 'SVDNB_npp_20150101-20151231_75N060W_vcm-orm-ntl_v10_c201701311200.avg_rade9.tif'
  character(*), parameter :: datafmt = '(F9.4,1X,F9.4,3X,F5.2,1X,F4.2,3X,F5.2)'

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

  call initpar(par)

  printcalib: block
    type(source), dimension(:,:), allocatable :: src

    allocate( src(size(mapi,1), size(mapi,2)) )
    call mksources(mapi, gti, maph, gth, src)

    open(34, file = 'calib.txt', action = 'write')
    write (34,'(2A7,A8,A6,A6,A10)') 'SCTR', 'ABS', 'HSCAL', 'SIDE', 'BG', 'DEV'

    iterate_parameters: do j = 1, size(par)
      !$omp parallel do private(sky,hobs)
      iterate_datpoints: do i = 1, ndata
        call onepoint(par(j), dat(i) % lat, dat(i) % lng, &
        &       src, maph, gth, sky, hobs)
        dat(i) % model = ity2mag(sky(4))
      end do iterate_datpoints
      !$omp end parallel do

      write (34,'(2F7.3,F8.1,F6.3,F6.2,F10.2)') &
        & par(j) % alpha, par(j) % relabs, &
        & par(j) % hscale, par(j) % beta, ity2mag(par(j) % skybg), &
        & sum((dat(1:ndata) % sky - dat(1:ndata) % model)**2)
    end do iterate_parameters

    close(34)

    deallocate(src)
  end block printcalib
  !
  ! do i = 1, ndata
  !   write (*,datafmt) dat(i) % lat, dat(i) % lng, &
  !     & dat(i) % sky, dat(i) % err, ity2mag(dat(i) % model)
  ! end do

  deallocate(mapi, maph, par)

contains

  subroutine initpar(par)
    type(modelpar), dimension(:), intent(inout), allocatable :: par

    integer :: i1,i2,i3,i4,i5,i
    real(fp), dimension(*), parameter :: par_alpha = [0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040]
    real(fp), dimension(*), parameter :: par_beta = [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
    real(fp), dimension(*), parameter :: par_skybg = [2.5e-7, 3e-7, 2e-7, 3.6e-7]
    real(fp), dimension(*), parameter :: par_hscale = [8000, 6000, 10000, 4000, 2000]
    real(fp), dimension(*), parameter :: par_relabs = [0.0, 0.4, 0.7, 1.0, 2.0, 4.0, 8.0, 12.0]

    allocate(par(size(par_alpha) * size(par_beta) * size(par_skybg) &
      & * size(par_hscale) * size(par_relabs)))

    i = 1
    do i3 = 1,size(par_skybg)
      do i2 = 1,size(par_beta)
        do i4 = 1,size(par_hscale)
          do i5 = 1,size(par_relabs)
            do i1 = 1,size(par_alpha)
              par(i) % alpha = par_alpha(i1)
              par(i) % beta = par_beta(i2)
              par(i) % skybg = par_skybg(i3)
              par(i) % hscale = par_hscale(i4)
              par(i) % relabs = par_relabs(i5)
              i = i + 1
            end do
          end do
        end do
      end do
    end do
  end subroutine initpar

end program sbcalib
