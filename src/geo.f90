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
    procedure, private :: gtimport_32
    procedure, private :: gtimport_64
    generic :: import => gtimport_32, gtimport_64
  end type

  interface gtinit
    module procedure :: gtimport_32
    module procedure :: gtimport_64
  end interface

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


  pure subroutine georay(lat1, lng1, h1, lat2, lng2, h2, t, latx, lngx, hx)
    real(fp), intent(in) :: lat1, lng1, lat2, lng2, h1, h2
    real(fp), dimension(:), intent(in) :: t
    real(fp), dimension(:), intent(out) :: latx, lngx, hx
    real(fp), dimension(3) :: x1, x2, xi
    integer :: i

    call geo2xyz(lat1, lng1, x1(1), x1(2), x1(3))
    call geo2xyz(lat2, lng2, x2(1), x2(2), x2(3))

    do concurrent (i = 1:size(t))
      xi = (1 - t(i)) * x1 * (1 + h1) + t(i) * x2 * (1 + h2)
      hx(i) = sqrt(xi(1)**2 + xi(2)**2 + xi(3)**2) - 1
      call xyz2geo(xi(1), xi(2), xi(3), latx(i), lngx(i))
    end do

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

    print '("PROJECTMAP",2I10)',  shape(mapsrc), shape(mapdest)

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

end module
