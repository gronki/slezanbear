module mapio

  use globals
  implicit none

contains

  !-----------------------------------------------------------------------------

  subroutine map_read_section(fn, llat, llng, ulat, ulng, map, gt)
    use iso_fortran_env, only: dp => real64
    use gdal
    use geo, only: geotransform, geo2arr, gtinit, arr2geo

    character(*), intent(in) :: fn
    real(dp), intent(in) :: llat, llng, ulat, ulng
    type(geotransform), intent(out) :: gt
    real(sp), intent(out), allocatable, dimension(:,:) :: map
    type(GDALDatasetH) :: gd
    type(GDALRasterBandH) :: band
    real(dp) :: g(6)
    integer :: nx, ny, nz, errno
    integer :: offx, offy, cnx, cny

    gd = GdalOpen(trim(fn) // char(0), GA_READONLY)
    if (.not. GDALAssociated(gd)) then
      error stop "could not open satellite map"
    end if

    write (0, '("opened ",A)') trim(fn)

    errno = GDALGetGeotransform(gd, g)
    print '("GDALGetGeotransform, errno = ",I0)', errno

    call gtinit(gt,g)

    write (0, '("geotransform = ",3F12.5)') g

    nx = GDALGetRasterXSize(gd)
    ny = GDALGetRasterYSize(gd)
    nz = GDALGetRasterCount(gd)

    write (0, '("nx = ",I0," ny = ",I0)') nx, ny

    if (nz .ne. 1) error stop "expected only one layer"

    determine_geometry: block
      real(dp) :: lli, lui, uui, uli
      real(dp) :: llj, luj, uuj, ulj
      real(dp) :: li , lj , ui , uj

      call geo2arr(gt, llat, llng, lli, llj)
      call geo2arr(gt, llat, ulng, lui, luj)
      call geo2arr(gt, ulat, ulng, uui, uuj)
      call geo2arr(gt, ulat, llng, uli, ulj)

      li = min(lli, lui, uui, uli)
      ui = max(lli, lui, uui, uli)
      lj = min(llj, luj, uuj, ulj)
      uj = max(llj, luj, uuj, ulj)

      if (     ui .le. 1   &
      &   .or. uj .le. 1   &
      &   .or. li .ge. nx  &
      &   .or. lj .ge. ny  ) then
        write (0,'(A)') "warning: this map does not contain any&
            & pixel in given range"
        return
      end if

      offx = max(0, floor(li) - 1)
      offy = max(0, floor(lj) - 1)

      cnx = min(nx, ceiling(ui)) - offx
      cny = min(ny, ceiling(uj)) - offy
    end block determine_geometry

    write (0, '("offset = ",2I6)') offx, offy
    write (0, '("window = ",2I6)') cnx, cny

    allocate( map(cnx,cny) )

    if (offx .ne. 0 .or. offy .ne. 0) then
      correct_geotransform: block
        real(dp) :: x0, y0
        call arr2geo(gt, real(offx + 1,dp), real(offy + 1,dp), y0, x0)
        gt % x0 = x0
        gt % y0 = y0
      end block correct_geotransform
    end if

    band = GDALGetRasterBand(gd,1)
    errno = GDALRasterIO_f(band, GF_READ, offx, offy, map)
    print '("GDALRasterIO_f, errno = ",I0)', errno
    if (errno .ne. 0) error stop "error reading GeoTIFF"

    call GDALClose(gd)
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine map_quicksave_r32(fn, map, gt, errno)
    use gdal
    use geo, only: geotransform, gtexport
    use iso_c_binding, only: c_null_ptr
    use iso_fortran_env, only: real32, real64

    character(*), intent(in) :: fn
    real(sp), dimension(:,:), contiguous, intent(inout) :: map
    type(geotransform), intent(in) :: gt
    integer, intent(inout) :: errno

    type(GDALDatasetH) :: gd
    type(GDALRasterBandH) :: band
    type(GDALDriverH) :: driver

    driver = GDALGetDriverByName('GTiff' // char(0))
    gd = GDALCreate(driver, trim(fn) // char(0), &
    &               size(map,1), size(map,2), 1, GDT_FLOAT32, c_null_ptr)
    band = GDALGetRasterBand(gd,1)

    try: block
      errno = GDALRasterIO_f(band, GF_WRITE, 0, 0, map)
      if (errno .ne. 0) exit try
      errno = GDALSetGeoTransform(gd, gtexport(gt))
      if (errno .ne. 0) exit try
    end block try

    call GDALClose(gd)
  end subroutine

  !-----------------------------------------------------------------------------

  ! subroutine read_map_data(mapi, gti, maph, gth)
  !   use gdal
  !   use geo
  !   real(sp), dimension(:,:), allocatable :: mapi, maph
  !   type(geotransform) :: gti, gth
  ! end subroutine read_map_data

  !-----------------------------------------------------------------------------

  subroutine map_read_tiles(fnlist, map, gt)
    use gdal
    use geo
    character(*), dimension(:) :: fnlist
    real(sp), dimension(:,:) :: map
    type(GDALDatasetH) :: gd
    type(geotransform) :: gt

    integer :: errno, i
    real(dp) o1,o2
    real(dp) :: llatsrc, llngsrc, ulatsrc, ulngsrc
    real(dp) :: llattil, llngtil, ulattil, ulngtil
    real(sp), dimension(:,:), allocatable :: maptile
    type(geotransform) :: gttile
    real(dp) :: g(6)

    call map_limits(size(map,1), size(map,2), gt, &
    &               llatsrc, llngsrc, ulatsrc, ulngsrc)

    do i = 1, size(fnlist)

      gd = GDALOpen(trim(fnlist(i)) // char(0), GA_READONLY)
      errno = GDALGetGeotransform(gd, g)
      call gttile % import(g)
      call map_limits(GDALGetRasterXSize(gd), GDALGetRasterYSize(gd), gttile, &
      &               llattil, llngtil, ulattil, ulngtil)
      call GDALClose(gd)

      o1 = overlap(llngtil, ulngtil, llngsrc, ulngsrc)
      o2 = overlap(llattil, ulattil, llatsrc, ulatsrc)

      if ( o1 > 0 .and. o2 > 0) then
        write (0,'(A32," overlap = ",F0.3)') fnlist(i), o1 * o2

        call map_read_section(fnlist(i), llatsrc, llngsrc, &
              & ulatsrc, ulngsrc, maptile, gttile)

        write (0,'("geotransform trimmed = ",3F15.4)')  gtexport(gttile)

        call projectmap(maptile, gttile, map, gt)
        if (allocated(maptile)) deallocate(maptile)
      else
        write (0,'(A32," overlap = ",A)') fnlist(i), 'none'
      end if
    end do

  contains

    elemental real(dp) function overlap(x1l,x1u,x2l,x2u)
      real(dp), intent(in) :: x1l, x1u
      real(dp), intent(in) :: x2l, x2u
      real(dp) :: c1, r1, c2, r2
      c1 =    (x1l + x1u) / 2
      r1 = abs(x1u - x1l) / 2
      c2 =    (x2l + x2u) / 2
      r2 = abs(x2u - x2l) / 2
      overlap = ((r1 + r2) - abs(c1 - c2)) / (r1 + r2)
    end function

  end subroutine map_read_tiles

end module
