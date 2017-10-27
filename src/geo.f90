module geo

  use globals
  use iso_fortran_env, only: r32 => real32, r64 => real64
  implicit none

  type geotransform
    !             g(1)    g(2)    g(3)
    real(r64) :: x0 = 0, xi = 1, xj = 0
    !             g(4)    g(5)    g(6)
    real(r64) :: y0 = 0, yi = 0, yj = 1
  contains
    procedure :: det => gtdet
    procedure :: export => gtexport
    procedure, private :: gtimport_32
    procedure, private :: gtimport_64
    generic :: import => gtimport_32, gtimport_64
  end type

  interface gtinit
    module procedure :: gtimport_32
    module procedure :: gtimport_64
  end interface

  type map_rect
    real(r64) :: llat, ulat
    real(r64) :: llng, ulng
  end type

contains

  elemental pure real(fp) function gtdet(gt) result(det)
    class(geotransform), intent(in) :: gt
    det = (gt % xi) * (gt % yj) - (gt % xj) * (gt % yi)
  end function

  pure subroutine gtimport_32(gt,g)
    class(geotransform), intent(inout) :: gt
    real(r32), dimension(6), intent(in) :: g
    gt % x0 = g(1)
    gt % xi = g(2)
    gt % xj = g(3)
    gt % y0 = g(4)
    gt % yi = g(5)
    gt % yj = g(6)
  end subroutine

  function gtexport(gt) result(g)
    class(geotransform), intent(in) :: gt
    real(fp), dimension(6) :: g
    g(1) = gt % x0
    g(2) = gt % xi
    g(3) = gt % xj
    g(4) = gt % y0
    g(5) = gt % yi
    g(6) = gt % yj
  end function

  pure subroutine gtimport_64(gt,g)
    class(geotransform), intent(inout) :: gt
    real(r64), dimension(6), intent(in) :: g
    gt % x0 = g(1)
    gt % xi = g(2)
    gt % xj = g(3)
    gt % y0 = g(4)
    gt % yi = g(5)
    gt % yj = g(6)
  end subroutine

  !-----------------------------------------------------------------------------

  elemental real(fp) function m2d(m) result(d)
    real(fp), intent(in) :: m
    d = (m / radearth) * (180 / pi)
  end function

  !-----------------------------------------------------------------------------

  elemental subroutine arr2geo(gt, i, j, lat, lng)
    class(geotransform), intent(in) :: gt
    real(fp), intent(in) :: i, j
    real(fp), intent(out) :: lat, lng
    lng = (gt % x0) + (gt % xi) * (i - 1) + (gt % xj) * (j - 1)
    lat = (gt % y0) + (gt % yi) * (i - 1) + (gt % yj) * (j - 1)
  end subroutine

  elemental subroutine geo2arr(gt, lat, lng, i, j)
    class(geotransform), intent(in) :: gt
    real(fp), intent(in) :: lat, lng
    real(fp), intent(out) :: i, j
    real(fp) :: det
    det = (gt % xi) * (gt % yj) - (gt % xj) * (gt % yi)
    i = (   (gt % yj) * (lng - (gt % x0)) &
    &     - (gt % xj) * (lat - (gt % y0)) ) / det + 1
    j = ( - (gt % yi) * (lng - (gt % x0)) &
    &     + (gt % xi) * (lat - (gt % y0)) ) / det + 1
  end subroutine

  !-----------------------------------------------------------------------------

  elemental subroutine geo2xyz(lat, lng, x, y, z)
    real(fp), intent(in) :: lat, lng
    real(fp), intent(out) :: x, y, z
    x = cos(lat / deg_in_rad) * cos(lng / deg_in_rad)
    y = cos(lat / deg_in_rad) * sin(lng / deg_in_rad)
    z = sin(lat / deg_in_rad)
  end subroutine

  elemental subroutine xyz2geo(x, y, z, lat, lng)
    real(fp), intent(in) :: x, y, z
    real(fp), intent(out) :: lat, lng
    lat = atan2(z, sqrt(x**2 + y**2)) * deg_in_rad
    lng = atan2(y, x) * deg_in_rad
  end subroutine

  !-----------------------------------------------------------------------------

  elemental real(fp) function angdist(lat1, lng1, lat2, lng2) result(rads)
    real(fp), intent(in) :: lat1, lng1
    real(fp), intent(in) :: lat2, lng2
    rads = acos(cos(lat1 / deg_in_rad) * cos(lat2 / deg_in_rad) &
          & * cos((lng1 - lng2) / deg_in_rad) &
          & + sin(lat1 / deg_in_rad) * sin(lat2 / deg_in_rad))
  end function

  elemental real function raydist(lat1, lng1, h1, lat2, lng2, h2) result(l)
    real(fp), intent(in) :: lat1, lng1, h1
    real(fp), intent(in) :: lat2, lng2, h2
    real(fp), dimension(3) :: x1, x2
    call geo2xyz(lat1, lng1, x1(1), x1(2), x1(3))
    call geo2xyz(lat2, lng2, x2(1), x2(2), x2(3))
    l = sqrt(sum((x1 * (1 + h1) - x2 * (1 + h2))**2))
  end function

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

  !-----------------------------------------------------------------------------

  pure subroutine interpolmap(map, gt, lat, lng, mout)
    real, intent(in), dimension(:,:) :: map
    real(fp), intent(in) :: lat, lng
    type(geotransform), intent(in) :: gt
    real(fp), intent(out) :: mout
    integer :: ki, kj
    real(fp) :: ri,rj
    real(fp) :: xi,xj

    call geo2arr(gt, lat, lng, xi, xj)

    ki = max(1, min(floor(xi), size(map,1) - 1))
    kj = max(1, min(floor(xj), size(map,2) - 1))

    ri = xi - ki
    rj = xj - kj

    if ( abs(2 * ri - 1) .le. 2 .and. abs(2 * rj - 1) .le. 2 ) then
      mout  = map(ki,   kj  ) * (1 - ri)  * (1 - rj)  &
      &     + map(ki+1, kj  ) * ri        * (1 - rj)  &
      &     + map(ki,   kj+1) * (1 - ri)  * rj        &
      &     + map(ki+1, kj+1) * ri        * rj
    end if

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine projectmap(mapsrc, gtsrc, mapdest, gtdest)
    real(sp), intent(in), dimension(:,:) :: mapsrc
    real(sp), intent(out), dimension(:,:) :: mapdest
    type(geotransform), intent(in) :: gtsrc, gtdest
    integer i,j
    integer ki, kj
    real(fp) lat, lng, ri, rj, xi, xj

    do concurrent (i = 1:size(mapdest,1), j = 1:size(mapdest,2))
      call arr2geo(gtdest, real(i,fp), real(j,fp), lat, lng)
      call geo2arr(gtsrc, lat, lng, xi, xj)

      ki = max(1, min(floor(xi), size(mapsrc,1) - 1))
      kj = max(1, min(floor(xj), size(mapsrc,2) - 1))

      ri = xi - ki
      rj = xj - kj

      if ( abs(2 * ri - 1) .le. 2 .and. abs(2 * rj - 1) .le. 2 ) then
        mapdest(i,j) = mapsrc(ki,   kj  ) * (1 - ri)  * (1 - rj)  &
        &            + mapsrc(ki+1, kj  ) * ri        * (1 - rj)  &
        &            + mapsrc(ki,   kj+1) * (1 - ri)  * rj        &
        &            + mapsrc(ki+1, kj+1) * ri        * rj
      end if
    end do

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine estimate_data_limits(lat1, lng1, lat2, lng2, dmax, llat, llng, ulat, ulng)

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

  subroutine map_limits(nx, ny, gt, llat, llng, ulat, ulng)
    integer, intent(in) :: nx, ny
    type(geotransform) :: gt
    real(fp), intent(out) :: llat, llng, ulat, ulng
    real(fp), dimension(4) :: lats, lngs

    call arr2geo(gt, real( 1,fp), real( 1,fp), lats(1), lngs(1))
    call arr2geo(gt, real(nx,fp), real( 1,fp), lats(2), lngs(2))
    call arr2geo(gt, real(nx,fp), real(ny,fp), lats(3), lngs(3))
    call arr2geo(gt, real( 1,fp), real(ny,fp), lats(4), lngs(4))

    llat = minval(lats)
    ulat = maxval(lats)
    llng = minval(lngs)
    ulng = maxval(lngs)

  end subroutine map_limits

  !-----------------------------------------------------------------------------

  subroutine map_gen_size(llat, llng, ulat, ulng, grid, nx, ny, gt)
    real(fp), intent(in) :: llat, llng, ulat, ulng
    real(fp), intent(in) :: grid
    integer, intent(out) :: nx, ny
    type(geotransform), intent(out) :: gt

    nx = nint(abs(ulng - llng) / m2d(grid) &
    &      * cos((llat + ulat) / (2 * deg_in_rad)))
    ny = nint(abs(ulat - llat) / m2d(grid))

    gt % x0 = llng
    gt % xi = (ulng - llng) / (nx - 1)
    gt % xj = 0
    gt % y0 = ulat
    gt % yi = 0
    gt % yj = (llat - ulat) / (ny - 1)
  end subroutine map_gen_size

end module
