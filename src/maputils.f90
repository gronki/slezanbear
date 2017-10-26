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
    call growmask(mask)
    llims(1,1) = minval(latm, mask)
    ulims(1,1) = maxval(latm, mask)
    llims(1,2) = minval(lngm, mask)
    ulims(1,2) = maxval(lngm, mask)

    !$omp section
    mask(:,:) = angdist(latm, lngm, lat1, lng2) .le. (dmax / radearth)
    call growmask(mask)
    llims(2,1) = minval(latm, mask)
    ulims(2,1) = maxval(latm, mask)
    llims(2,2) = minval(lngm, mask)
    ulims(2,2) = maxval(lngm, mask)

    !$omp section
    mask(:,:) = angdist(latm, lngm, lat2, lng2) .le. (dmax / radearth)
    call growmask(mask)
    llims(3,1) = minval(latm, mask)
    ulims(3,1) = maxval(latm, mask)
    llims(3,2) = minval(lngm, mask)
    ulims(3,2) = maxval(lngm, mask)

    !$omp section
    mask(:,:) = angdist(latm, lngm, lat2, lng1) .le. (dmax / radearth)
    call growmask(mask)
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

  contains

    subroutine growmask(mask)
      logical, dimension(:,:), intent(inout) :: mask
      logical, dimension(:,:), allocatable :: mask2
      integer i, j, li, ui, lj, uj
      mask2 = mask
      do j = 1, size(mask,2)
        do i = 1, size(mask,1)
          li = max(1,            i - 1)
          ui = min(size(mask,1), i + 1)
          lj = max(1,            j - 1)
          uj = min(size(mask,2), j + 1)
          mask(i,j) = any(mask2(li:ui, lj:uj))
        end do
      end do
    end subroutine

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

  !-----------------------------------------------------------------------------

  subroutine gdal_quicksave_r32(fn, map, gt, errno)
    use gdal
    use iso_c_binding, only: c_null_ptr
    use iso_fortran_env, only: real32, real64

    character(*), intent(in) :: fn
    real(sp), dimension(:,:), intent(inout) :: map
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
      errno = GDALSetGeoTransform(gd, [gt % x0, gt % xi, gt % xj, &
                                    &  gt % y0, gt % yi, gt % yj])
      if (errno .ne. 0) exit try
    end block try

    call GDALClose(gd)
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine binmap(map, gt, b)
    real(sp), dimension(:,:), allocatable, intent(inout) :: map
    type(geotransform), intent(inout) :: gt
    integer, intent(in) :: b
    real(sp), dimension(:,:), allocatable :: buf
    real(fp) :: x0, y0
    integer :: i,j

    allocate(buf(floor(size(map,1) / real(b)), &
    &            floor(size(map,2) / real(b))))

    forall (i = 1:size(buf,1), j = 1:size(buf,2))
      buf(i,j) = sum(map(1 + (i-1) * b : i * b, 1 + (j-1) * b : j * b)) / b**2
    end forall

    call arr2geo(gt, (b + 1d0) / 2, (b + 1d0) / 2, y0, x0)

    gt % x0 = x0
    gt % y0 = y0
    gt % xi = b * (gt % xi)
    gt % xj = b * (gt % xj)
    gt % yi = b * (gt % yi)
    gt % yj = b * (gt % yj)

    deallocate(map)
    call move_alloc(buf, map)

  end subroutine

end module maputils
