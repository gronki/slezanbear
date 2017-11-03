program sbcalib

  use globals
  use gdal
  use mapio
  use slezanbear
  use geo

  implicit none

  type dataentry
    real(fp) :: lat, lng
    real(fp) :: sky
  end type

  real(fp) :: llatsrc, llngsrc, ulatsrc, ulngsrc
  type(source), dimension(:,:), allocatable :: src

  type(dataentry), dimension(1024), target :: dat_all
  type(dataentry), dimension(:), pointer :: dat
  integer :: ndata, i, j

  character(128), dimension(9) :: fndem

  integer, parameter :: npar = 5
  integer, parameter :: nprobes = npar + 1
  real(fp), dimension(npar), parameter :: parlo &
  &   = [ 3e-6,  0.0, 4e3, 0.0, 3.0e-7]
  real(fp), dimension(npar), parameter :: parhi &
  &   = [ 9e-6, 24.0, 9e3, 0.4, 4.0e-7]
  real(fp), dimension(npar, nprobes) :: par
  real(fp), dimension(nprobes) :: modelerr
  real(fp), dimension(:,:), allocatable :: modelsky

  real(sp), dimension(:,:), allocatable :: mapi, maph
  type(geotransform) :: gti, gth
  integer :: errno

  character(*), parameter :: fni = 'SVDNB_npp_20150101-20151231_75N060W_vcm-orm-ntl_v10_c201701311200.avg_rade9.tif'
  character(*), parameter :: datafmt = '(F9.4,1X,F9.4,3X,F5.2,1X,F4.2,3X,F5.2)'

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

  printcalib: block
    type(modelpar) :: mpar
    real(fp) :: sky(4), hobs, a(npar,nprobes)
    integer :: ibest, iworst, k

    allocate( src(size(mapi,1), size(mapi,2)) )
    call mksources(mapi, gti, maph, gth, src)

    call random_number(par)

    do i = 1, npar
      par(i,:) = parlo(i) + par(i,:) * (parhi(i) - parlo(i))
    end do

    open(34, file = 'calib.log.txt', action = 'write')

    simplex: do k = 1,16

      iterate_modelpoints: do j = 1, nprobes
        mpar % alpha  = par(1,j)
        mpar % relabs = par(2,j)
        mpar % hscale = par(3,j)
        mpar % beta   = par(4,j)
        mpar % skybg  = par(5,j)

        !$omp parallel do private(i, sky, hobs)
        iterate_datpoints: do i = 1, ndata
          call onepoint(mpar, dat(i) % lat, dat(i) % lng, &
          &       src, maph, gth, sky, hobs)
          modelsky(i,j) = ity2mag(sky(4))
        end do iterate_datpoints
        !$omp end parallel do

        modelerr(j) = sqrt(sum((modelsky(:,j) - dat % sky)**2) / ndata)
      end do iterate_modelpoints

      print '(*(F12.4))', modelerr
      ! if (any(modelerr < 0.22)) exit simplex

      iworst = maxloc(modelerr,1)
      ibest = minloc(modelerr,1)

      write (34, '(F12.3, F12.2, F12.1, F12.4, F12.3, F12.4)') &
      &       par(1,ibest) * 1e6, par(2,ibest), par(3,ibest), par(4,ibest), &
      &       par(5,ibest) * 1e7, modelerr(ibest)

      call random_number(a)
      do concurrent (i = 1:npar, j = 1:nprobes, j /= ibest)
        par(i,j) = par(i,ibest) * (1 + (2 * a(i,j) - 1) / sqrt(4.0*k))
      end do

    end do simplex

    close(34)

    print '(*(F12.3))', par(1,:) * 1e6
    print '(*(F12.2))', par(2,:)
    print '(*(F12.1))', par(3,:)
    print '(*(F12.4))', par(4,:)
    print '(*(F12.3))', par(5,:) * 1e7

    ibest = minloc(modelerr,1)
    verify_datpoints: do i = 1, ndata
      write (*,'(2F10.4,3F10.3)') &
      &       dat(i) % lat, dat(i) % lng, dat(i) % sky, modelsky(i,ibest), &
      &       modelsky(i,ibest) - dat(i) % sky
    end do verify_datpoints

    deallocate(src)
  end block printcalib

  deallocate(mapi, maph)

end program sbcalib
