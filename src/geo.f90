module geo

  use globals
  implicit none

  type geotransform
    !       g(1)  g(2)  g(3)
    real :: x0,   xi,    xj
    !       g(4)  g(5)  g(6)
    real :: y0,   yi,    yj
  contains
    procedure :: det => geotransform_det
    procedure :: import => geotransform_import
  end type

contains

  elemental pure real function geotransform_det(gt) result(det)
    class(geotransform), intent(in) :: gt
    det = (gt % xi) * (gt % yj) - (gt % xj) * (gt % yi)
  end function

  pure subroutine geotransform_import(gt,g)
    class(geotransform), intent(inout) :: gt
    real, dimension(6), intent(in) :: g
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
    real, intent(in) :: i, j
    real, intent(out) :: lat, lng
    lng = (gt % x0) + (gt % xi) * (i - 1) + (gt % xj) * (j - 1)
    lat = (gt % y0) + (gt % yi) * (i - 1) + (gt % yj) * (j - 1)
  end subroutine

  elemental subroutine geo2arr(gt, lat, lng, i, j)
    class(geotransform), intent(in) :: gt
    real, intent(in) :: lat, lng
    real, intent(out) :: i, j
    real :: det
    det = (gt % xi) * (gt % yj) - (gt % xj) * (gt % yi)
    i = (   (gt % yj) * (lng - (gt % x0)) &
    &     - (gt % xj) * (lat - (gt % y0)) ) / det + 1
    j = ( - (gt % yi) * (lng - (gt % x0)) &
    &     + (gt % xi) * (lat - (gt % y0)) ) / det + 1
  end subroutine

  !-----------------------------------------------------------------------------

  elemental subroutine geo2xyz(lat, lng, x, y, z)
    real, intent(in) :: lat, lng
    real, intent(out) :: x, y, z
    x = cos(lat / deg_in_rad) * cos(lng / deg_in_rad)
    y = cos(lat / deg_in_rad) * sin(lng / deg_in_rad)
    z = sin(lat / deg_in_rad)
  end subroutine

  elemental subroutine xyz2geo(x, y, z, lat, lng)
    real, intent(in) :: x, y, z
    real, intent(out) :: lat, lng
    lat = atan2(z, sqrt(x**2 + y**2)) * deg_in_rad
    lng = atan2(y, x) * deg_in_rad
  end subroutine

  !-----------------------------------------------------------------------------

  elemental real function angdist(lat1, lng1, lat2, lng2) result(rads)
    real, intent(in) :: lat1, lng1
    real, intent(in) :: lat2, lng2
    rads = acos(cos(lat1 / deg_in_rad) * cos(lat2 / deg_in_rad) &
          & * cos((lng1 - lng2) / deg_in_rad) &
          & + sin(lat1 / deg_in_rad) * sin(lat2 / deg_in_rad))
  end function

  elemental real function raydist(lat1, lng1, h1, lat2, lng2, h2) result(l)
    real, intent(in) :: lat1, lng1, h1
    real, intent(in) :: lat2, lng2, h2
    real, dimension(3) :: x1, x2
    call geo2xyz(lat1, lng1, x1(1), x1(2), x1(3))
    call geo2xyz(lat2, lng2, x2(1), x2(2), x2(3))
    l = sqrt(sum((x1 * (radearth + h1) - x2 * (radearth + h2))**2))
  end function

  !-----------------------------------------------------------------------------


  pure subroutine georay(lat1, lng1, h1, lat2, lng2, h2, t, latx, lngx, hx)
    real, intent(in) :: lat1, lng1, lat2, lng2, h1, h2
    real, dimension(:), intent(in) :: t
    real, dimension(:), intent(out) :: latx, lngx, hx
    real, dimension(3) :: x1, x2, xi
    integer :: i

    call geo2xyz(lat1, lng1, x1(1), x1(2), x1(3))
    call geo2xyz(lat2, lng2, x2(1), x2(2), x2(3))

    x1 = x1 * (radearth + h1)
    x2 = x2 * (radearth + h2)

    do concurrent (i = 1:size(t))
      xi = (1 - t(i)) * x1 + t(i) * x2
      hx(i) = sqrt(xi(1)**2 + xi(2)**2 + xi(3)**2) - radearth
      call xyz2geo(xi(1), xi(2), xi(3), latx(i), lngx(i))
    end do

  end subroutine

  !-----------------------------------------------------------------------------

  pure subroutine interpolmap(map, gt, lat, lng, mout)
    real, intent(in), dimension(:,:) :: map
    real, intent(in) :: lat, lng
    type(geotransform), intent(in) :: gt
    real, intent(out) :: mout
    integer :: ki, kj
    real :: ri,rj,xi,xj

    call geo2arr(gt, lat, lng, xi, xj)

    ri = xi - floor(xi)
    rj = xj - floor(xj)

    ki = nint(xi - ri)
    kj = nint(xj - rj)

    if ( ki >= 1 .and. ki <= size(map,1) - 1 &
        & .and. kj >= 1 .and. kj <= size(map,2) - 1 ) then
      mout    = map(ki,   kj  ) * (1 - ri)  * (1 - rj)  &
      &       + map(ki+1, kj  ) * ri        * (1 - rj)  &
      &       + map(ki,   kj+1) * (1 - ri)  * rj        &
      &       + map(ki+1, kj+1) * ri        * rj
    end if

  end subroutine

end module
