program sbmap_p

  use globals
  use gdal
  use mapio
  use slezanbear
  use geo

  implicit none

  real(dp) :: llatsrc, llngsrc, ulatsrc, ulngsrc
  type(modelpar) :: par

  character(128), dimension(9) :: fndem

  real(sp), dimension(:,:), allocatable :: mapi, maph
  type(geotransform) :: gti, gth

  character(*), parameter :: fni = 'SVDNB_npp_20150101-20151231_75N060W_vcm-orm-ntl_v10_c201701311200.avg_rade9.tif'

  fndem(1) = "eudem_dem_5deg_n45e010.tif"
  fndem(2) = "eudem_dem_5deg_n45e015.tif"
  fndem(3) = "eudem_dem_5deg_n45e020.tif"
  fndem(4) = "eudem_dem_5deg_n50e010.tif"
  fndem(5) = "eudem_dem_5deg_n50e015.tif"
  fndem(6) = "eudem_dem_5deg_n50e020.tif"
  fndem(7) = "eudem_dem_5deg_n55e010.tif"
  fndem(8) = "eudem_dem_5deg_n55e015.tif"
  fndem(9) = "eudem_dem_5deg_n55e020.tif"

  !! dmax = 1.05 * radearth * acos(1 / (1 + 5 * par % hscale / radearth))
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


  compute_model: block
    type(GDALDatasetH) :: gd
    type(GDALRasterBandH) :: band
    type(GDALDriverH) :: driver

    integer :: errno
    type(geotransform) :: gtm
    real(dp), dimension(:,:), allocatable :: hobs
    real(dp), dimension(:,:,:), allocatable :: sky
    integer nox,noy,i
    real(dp) :: t0,t1

    call map_gen_size(llat, llng, ulat, ulng, model_grid_meters, nox, noy, gtm)
    allocate( sky(nox,noy,4), hobs(nox,noy) )
    write (*,'("model raster size = ",2I6)') nox, noy

    par % alpha = 5.678e-6
    par % relabs = 9.072
    par % hscale = 7026
    par % beta = 0.112
    par % skybg = 3.56e-7

    call cpu_time(t0)
    call sbmap(par, mapi, gti, maph, gth, sky, hobs, gtm)
    call cpu_time(t1)

    write (*, '("CPUTIME = ",F0.2," mins")') (t1 - t0) / 60

    driver = GDALGetDriverByName('GTiff' // char(0))
    gd = GDALCreate(driver, trim(outfn) // ".tif" // char(0), &
    &               nox, noy, 5, GDT_FLOAT64, c_null_ptr)

    band = GDALGetRasterBand(gd,1)
    errno = GDALRasterIO_f(band, GF_WRITE, 0, 0, hobs)
    print '("GDALRasterIO_f, errno = ",I0)', errno

    do i = 1,4
      band = GDALGetRasterBand(gd,i + 1)
      errno = GDALRasterIO_f(band, GF_WRITE, 0, 0, sky(:,:,i))
      print '("GDALRasterIO_f, errno = ",I0)', errno
    end do

    errno = GDALSetGeoTransform(gd, gtexport(gtm))
    print '("GDALSetGeoTransform, errno = ",I0)', errno
    call GDALClose(gd)

    deallocate( sky, hobs )
  end block compute_model

  deallocate(mapi, maph)


end program sbmap_p
