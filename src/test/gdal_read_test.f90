program gdal_read_test

  use gdal
  use geo
  use iso_c_binding
  use gdal_test_common

  implicit none

  type(GDALDriverH) :: driver
  type(GDALDatasetH) :: gd
  type(GDALRasterBandH) :: band

  integer :: nx, ny, nz, errno, i
  real(8) :: g(6)
  real(8), dimension(:,:,:), allocatable :: dat
  type(geotransform) :: gt

  call GDALAllRegister

  driver = GDALGetDriverByName('GTiff' // char(0))
  gd = GdalOpen('test.tif' // char(0), GA_READONLY)

  if ( .not.GdalAssociated(gd) ) then
    error stop "error opening dataset"
  end if

  errno = GDALGetGeotransform(gd, g)
  print '("errno = ",I0)', errno
  print *,g
  call gt % import(g)

  nx = GDALGetRasterXSize(gd)
  ny = GDALGetRasterYSize(gd)
  nz = GDALGetRasterCount(gd)

  print '("nx = ",I0,"  ny = ",I0,"  nz = ",I0)', nx, ny, nz

  allocate(dat(nx,ny,nz))

  do i = 1,3
    band = GDALGetRasterBand(gd,i)
    errno = GDALRasterIO_f(band, GF_READ, 0, 0, dat(:,:,i))
    print '("band ",I0," errno = ",I0)', i,errno
  end do

  call checkdata(gt, dat, 1,   1)
  call checkdata(gt, dat, 120, 1)
  call checkdata(gt, dat, 120, 80)
  call checkdata(gt, dat, 1,   80)
  call checkdata(gt, dat, 60,  40)
  call checkdata(gt, dat, 61,  40)
  call checkdata(gt, dat, 61,  41)
  call checkdata(gt, dat, 60,  41)

  call GDALClose(gd)
  deallocate(dat)

end program
