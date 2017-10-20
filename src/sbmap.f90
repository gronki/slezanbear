program sbmap_p

  use globals
  use gdal, only: GDALDriverH, GDALAllRegister, GDALGetDriverByName
  use maputils
  use slezanbear
  use geo

  implicit none

  real(fp) :: llatsrc, llngsrc, ulatsrc, ulngsrc
  type(modelpar) :: par

  integer :: ntiles
  type(tileinfo), dimension(64) :: tiles
  type(GDALDriverH) :: driver

  real(sp), dimension(:,:), allocatable :: mapi, maph
  type(geotransform) :: gti, gth

  character(*), parameter :: fni = 'SVDNB_npp_20150101-20151231_75N060W_vcm-orm-ntl_v10_c201701311200.avg_rade9.tif'

  !! dmax = 1.05 * radearth * acos(1 / (1 + 5 * par % hscale / radearth))
  print '("dmax = ",F0.2," km")', dmax / 1e3
  call estimate_data_limits(llat, llng, ulat, ulng, dmax * 1.05, &
  &       llatsrc, llngsrc, ulatsrc, ulngsrc)

  print '("data limits = ",F6.2,1X,F7.2)', llatsrc, llngsrc, ulatsrc, ulngsrc

  readtiles: block
    integer :: errno
    real(fp) o1,o2
    type(tileinfo) :: til
    open(30, file = 'demlist.txt', action = 'read')
    ntiles = 0
    do while (.true.)
      call til % read_(30, errno)
      if ( errno .ne. 0 ) exit
      o1 = overlap(til % lx, til % ux, llngsrc, ulngsrc)
      o2 = overlap(til % ly, til % uy, llatsrc, ulatsrc)
      if ( o1 > 0 .and. o2 > 0) then
        ntiles = ntiles + 1
        write (0,'(I4,1X,A32," overlap = ",F0.3)') ntiles, til % fn, o1 * o2
        tiles(ntiles) = til
      end if
    end do
    close(30)
  end block readtiles

  call GDALAllRegister
  driver = GDALGetDriverByName('GTiff' // char(0))

  call gdal_read_section(fni, llatsrc, llngsrc, ulatsrc, ulngsrc, mapi, gti)
  !call binmap(mapi, gti, 3)

  write (0, '("satellite raster size = ", 2I6)') size(mapi,1), size(mapi,2)

  allocate( maph( &
  &   nint(abs(ulngsrc - llngsrc) / m2d(elev_grid_meters) &
  &                               * cos(llat/2 + ulat/2)), &
  &   nint(abs(ulatsrc - llatsrc) / m2d(elev_grid_meters)) &
  &  ))

  write (*,'("dem raster size = ",2I6)') size(maph,1), size(maph,2)

  gth % x0 = llngsrc
  gth % xi = (ulngsrc - llngsrc) / (size(maph,1) - 1)
  gth % xj = 0
  gth % y0 = ulatsrc
  gth % yi = 0
  gth % yj = -(ulatsrc - llatsrc) / (size(maph,2) - 1)

  loadtiles: block
    integer i
    real(sp), dimension(:,:), allocatable :: maptile
    type(geotransform) :: gttile
    do i = 1,ntiles
      call gdal_read_section(tiles(i) % fn, llatsrc, llngsrc, &
            & ulatsrc, ulngsrc, maptile, gttile)
      write (0,'("geotransform trimmed = ",3F15.4)') gttile % x0, gttile % xi, &
       gttile % xj, gttile % y0, gttile % yi, gttile % yj
      call projectmap(maptile, gttile, maph, gth)
      if (allocated(maptile)) deallocate(maptile)
    end do
  end block loadtiles

  compute_model: block
    use gdal
    type(GDALDatasetH) :: gd
    type(GDALRasterBandH) :: band
    integer :: errno
    type(geotransform) :: gtm
    real(fp), dimension(:,:), allocatable :: I1, I2, hobs
    integer nox,noy
    real(fp) :: t0,t1

    nox = nint(abs(ulng - llng) / m2d(model_grid_meters) * cos(llat/2 + ulat/2))
    noy = nint(abs(ulat - llat) / m2d(model_grid_meters))
    allocate( I1(nox,noy), I2(nox,noy), hobs(nox,noy) )

    write (*,'("model raster size = ",2I6)') nox, noy

    gtm % x0 = llng
    gtm % xi = (ulng - llng) / (nox - 1)
    gtm % xj = 0
    gtm % y0 = ulat
    gtm % yi = 0
    gtm % yj = (llat - ulat) / (noy - 1)

    call cpu_time(t0)
    call sbmap(par, mapi, gti, maph, gth, I1, I2, hobs, gtm)
    call cpu_time(t1)

    write (*, '("CPUTIME = ",F0.2)') t1 - t0

    gd = GDALCreate(driver, trim(outfn) // char(0), &
    &               nox, noy, 3, GDT_FLOAT64, c_null_ptr)

    band = GDALGetRasterBand(gd,1)
    errno = GDALRasterIO_f(band, GF_WRITE, 0, 0, I1)
    print '("GDALRasterIO_f, errno = ",I0)', errno

    band = GDALGetRasterBand(gd,2)
    errno = GDALRasterIO_f(band, GF_WRITE, 0, 0, I2)
    print '("GDALRasterIO_f, errno = ",I0)', errno

    band = GDALGetRasterBand(gd,3)
    errno = GDALRasterIO_f(band, GF_WRITE, 0, 0, hobs)
    print '("GDALRasterIO_f, errno = ",I0)', errno

    errno = GDALSetGeoTransform(gd, [gtm % x0, gtm % xi, gtm % xj, &
          & gtm % y0, gtm % yi, gtm % yj])
    print '("GDALSetGeoTransform, errno = ",I0)', errno
    call GDALClose(gd)

    deallocate( I1, I2, hobs )
  end block compute_model

  deallocate(mapi, maph)

contains

  elemental real(fp) function overlap(x1l,x1u,x2l,x2u)
    real(fp), intent(in) :: x1l, x1u
    real(fp), intent(in) :: x2l, x2u
    real(fp) :: c1, r1, c2, r2
    c1 =    (x1l + x1u) / 2
    r1 = abs(x1u - x1l) / 2
    c2 =    (x2l + x2u) / 2
    r2 = abs(x2u - x2l) / 2
    overlap = ((r1 + r2) - abs(c1 - c2)) / (r1 + r2)
  end function

  elemental real(fp) function m2d(m) result(d)
    real(fp), intent(in) :: m
    d = (m / radearth) * (180 / pi)
  end function

end program sbmap_p
