module maputils

  use globals
  use geo
  use gdal
  implicit none

  character(*), parameter :: tileinfo_fmt = '(A32,1X,I6,1X,I6,4(1X,F9.4))'

  type tileinfo
    character(32) :: fn
    integer :: nx = 0, ny = 0
    real(fp) :: lx, ux
    real(fp) :: ly, uy
  contains
    procedure :: read_ => tileinfo_read
  end type tileinfo

contains

  subroutine tileinfo_read(ti,u,errno)
    class(tileinfo) :: ti
    integer :: u,errno
    read (u, tileinfo_fmt, iostat = errno) ti % fn, &
          & ti % nx, ti % ny, &
          & ti % lx, ti % ux, &
          & ti % ly, ti % uy
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine estimate_data_limits(lat1, lng1, lat2, lng2, dmax, llat, llng, ulat, ulng)

    use geo, only: radearth, angdist

    real(fp), intent(in) :: lat1, lng1
    real(fp), intent(in) :: lat2, lng2
    real(fp), intent(in) :: dmax
    real(fp), intent(out) :: llat, llng, ulat, ulng
    integer, parameter :: p = 5
    real(fp), dimension(:,:), allocatable :: latm, lngm
    logical,  dimension(:,:), allocatable :: mask
    real(fp), dimension(4,2) :: llims, ulims
    integer i,j

    allocate( latm(p*180,p*360), lngm(p*180,p*360), mask(p*180,p*360) )

    forall (i = 1:size(latm,1), j = 1:size(latm,2))
      latm(i,j) = real(i - 1, fp) / p - 90
      lngm(i,j) = (lng1 + lng2) / 2 + real(j - 1, fp) / p - 180
    end forall

    !$omp parallel sections private(mask)

    !$omp section
    mask(:,:) = angdist(latm, lngm, lat1, lng1) .le. (dmax / radearth)
    llims(1,1) = minval(latm, mask)
    ulims(1,1) = maxval(latm, mask)
    llims(1,2) = minval(lngm, mask)
    ulims(1,2) = maxval(lngm, mask)

    !$omp section
    mask(:,:) = angdist(latm, lngm, lat1, lng2) .le. (dmax / radearth)
    llims(2,1) = minval(latm, mask)
    ulims(2,1) = maxval(latm, mask)
    llims(2,2) = minval(lngm, mask)
    ulims(2,2) = maxval(lngm, mask)

    !$omp section
    mask(:,:) = angdist(latm, lngm, lat2, lng2) .le. (dmax / radearth)
    llims(3,1) = minval(latm, mask)
    ulims(3,1) = maxval(latm, mask)
    llims(3,2) = minval(lngm, mask)
    ulims(3,2) = maxval(lngm, mask)

    !$omp section
    mask(:,:) = angdist(latm, lngm, lat2, lng1) .le. (dmax / radearth)
    llims(4,1) = minval(latm, mask)
    ulims(4,1) = maxval(latm, mask)
    llims(4,2) = minval(lngm, mask)
    ulims(4,2) = maxval(lngm, mask)

    !$omp end parallel sections

    llat = minval(llims(:,1))
    llng = minval(llims(:,2))
    ulat = maxval(ulims(:,1))
    ulng = maxval(ulims(:,2))

    deallocate( latm, lngm, mask )

  end subroutine estimate_data_limits

  !-----------------------------------------------------------------------------

  subroutine gdal_read_section(fn, llat, llng, ulat, ulng, map, gt)
    use iso_fortran_env, only: dp => real64
    use gdal
    character(*), intent(in) :: fn
    real(fp), intent(in) :: llat, llng, ulat, ulng
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
        real(fp) :: x0, y0
        call arr2geo(gt, real(offx + 1,fp), real(offy + 1,fp), y0, x0)
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

end module maputils
