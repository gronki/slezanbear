program gdal_write_test

  use gdal
  use geo
  use iso_c_binding
  use gdal_test_common

  implicit none

  type(GDALDriverH) :: driver
  type(GDALDatasetH) :: gd
  type(GDALRasterBandH) :: band
  integer, parameter :: nx = 120, ny = 80, nz = 3
  real(8), dimension(nx,ny,nz) :: dat
  real(8) :: g(6)
  integer :: errno, i

  g = [16.9, (17.1 - 16.9) / (nx - 1), 0.0, &
  &    50.9, 0.0, (51.1 - 50.9) / (ny - 1)]

  call GDALAllRegister

  driver = GDALGetDriverByName('GTiff' // char(0))
  gd = GDALCreate(driver, 'test.tif' // char(0), nx, ny, nz, GDT_FLOAT64, c_null_ptr)

  print '("nx = ",I0,"  ny = ",I0,"  nz = ",I0)', nx, ny, nz

  gen_image : block

    real(8) :: lat, lng, dist
    type(geotransform) :: gt
    integer :: i,j

    call gt % import(g)

    do concurrent (i = 1:nx, j = 1:ny)
      call arr2geo(gt, real(i,8), real(j,8), lat, lng)
      dist = radearth * angdist(lat, lng, 51.0d0, 17.0d0)
      dat(i,j,1) = lat
      dat(i,j,2) = lng
      dat(i,j,3) = dist
    end do

    call checkdata(gt, dat, 1,   1)
    call checkdata(gt, dat, 120, 1)
    call checkdata(gt, dat, 120, 80)
    call checkdata(gt, dat, 1,   80)
    call checkdata(gt, dat, 60,  40)
    call checkdata(gt, dat, 61,  40)
    call checkdata(gt, dat, 61,  41)
    call checkdata(gt, dat, 60,  41)

  end block gen_image


  do i = 1,3
    band = GDALGetRasterBand(gd, i)
    errno = GDALRasterIO_f(band, GF_WRITE, 0, 0, dat(:,:,i))
    print '("band ",I0," errno = ",I0)', i,errno
  end do
  ! errno = GDALSetProjection(gd, "WGS 84" // char(0))
  ! print '("errno = ",I0)', errno
  errno = GDALSetGeoTransform(gd, g)
  print '("errno = ",I0)', errno
  call GDALClose(gd)

end program gdal_write_test
