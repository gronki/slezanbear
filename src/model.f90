module slezanbear

  use globals
  use geo
  use kernel

  implicit none

  type lpmap
    real(sp) :: lat, lng, h
    real(sp) :: iz, tau
  end type lpmap

  logical :: terrain_attenuation = .true.

contains

  !-----------------------------------------------------------------------------

  pure real function integrate(y, x) result(integral)
    real, intent(in), dimension(:) :: x, y
    integral = sum((y(2:size(y)) + y(1:size(y)-1)) &
              &  * (x(2:size(x)) - x(1:size(x)-1))) / 2
  end function

  !-----------------------------------------------------------------------------

  pure subroutine checkline(maph, gth, lat1, lng1, h1, lat2, lng2, h2, attn)
    ! altitude map and geotransform
    real, dimension(:,:), intent(in) :: maph
    type(geotransform), intent(in) :: gth
    ! two points and their relative heights (h / R)
    real, intent(in) :: lat1, lng1, h1
    real, intent(in) :: lat2, lng2, h2
    ! result: attenuation factor due to terrain (0 - obstruction, 1 - no obst)
    real, intent(out) :: attn
    ! internal variables
    integer, parameter :: nsectmax = 2**10
    integer :: nsect,i
    real, parameter :: smooth_meters = 2
    real, dimension(nsectmax) :: t, latx, lngx, hx, hterr

    nsect = int(raydist(lat1, lng1, h1, lat2, lng2, h2) / 10)
    if (nsect > nsectmax) nsect = nsectmax
    if (nsect < 4) nsect = 4

    forall (i = 1:nsect) t(i) = (i - 1) / real(nsect - 1)

    call georay(lat1, lng1, h1, lat2, lng2, h2, t(1:nsect), &
          & latx(1:nsect), lngx(1:nsect), hx(1:nsect))

    do concurrent (i = 1:nsect)
      call interpolmap(maph, gth, latx(i), lngx(i), hterr(i))
    end do

    attn = th(minval(hx(1:nsect) - hterr(1:nsect)), smooth_meters / 2)

  end subroutine checkline

  !-----------------------------------------------------------------------------

  subroutine sbmlist(par, lat, lng, mapi, gti, maph, gth, Iobs, tau)
    ! model parameters
    type(modelpar), intent(in) :: par
    ! coordinates of the point
    real, intent(in), dimension(:) :: lat, lng
    ! altitude and flux maps and respective geotransforms
    real, dimension(:,:), intent(in) :: mapi, maph
    type(geotransform), intent(in) :: gti, gth
    ! output
    real, intent(out), dimension(:) :: Iobs, tau
    ! internal variables
    type(modelpar) :: par1
    real, dimension(size(mapi,1),size(mapi,2)) :: asrc, latsrc, lngsrc, hsrc
    real :: hobs
    integer i,j,k,l
    integer, parameter :: nh = 500
    real, dimension(nh) :: hsct, J0, J1, Jatt, JJ

    ! interpolate altitude for all source points and compute their coordinates
    do concurrent (i = 1:size(mapi,1), j = 1:size(mapi,2))
      call arr2geo(gti, real(i), real(j), latsrc(i,j), lngsrc(i,j))

      asrc(i,j) = abs(gti % det()) * (radearth / deg_in_rad)**2 &
            & * cos((latsrc(i,j)) / deg_in_rad)

      call interpolmap(maph, gth, latsrc(i,j), lngsrc(i,j), hsrc(i,j))
    end do

    iterate_list: do k = 1, size(lat)

      ! get the altitude of the observer (in meters)
      call interpolmap(maph, gth, lat(k), lng(k), hobs)

      ! optical depth of the entire atmosphere
      tau(k) = (par % alpha) * (par % relabs) * exp(-hobs / (par % hscale))
      ! background glow diminished by the absorption
      Iobs(k) = (par % skybg) * exp(-tau(k))

      par1 = par

      forall (i = 1:nh)
        hsct(i) = hobs + (6 * (par % hscale) - hobs) * real(i-1) / (nh - 1)
      end forall

      ! iterate through all map points
      do concurrent (i = 1:size(mapi,1), j = 1:size(mapi,2))
        if (mapi(i,j) .eq. 0) cycle

        par1 % alpha = (par % alpha) * exp(-(hsrc(i,j)) / (par % hscale))

        call kern(par1, asrc(i,j), radearth + hsrc(i,j), &
            & radearth * angdist(lat, lng, latsrc(i,j), lngsrc(i,j)), &
            & hsct - hsrc(i,j), hobs - hsrc(i,j), J0, J1)

        if (terrain_attenuation) then
          do concurrent (l = 1:nh)
            call checkline(maph, gth,                       &
                  & lat(k),       lng(k),       hsct(l),    &
                  & latsrc(i,j),  lngsrc(i,j),  hsrc(i,j),  &
                  & Jatt(l))
          end do
          JJ(:) = J0(:) * J1(:) * Jatt(:)
        else
          JJ(:) = J0(:) * J1(:)
        end if

        Iobs(k) = Iobs(k) + INTEGRATE(JJ,hsct) * mapi(i,j) * 1e-9
      end do

    end do iterate_list


  end subroutine

  subroutine onepoint_c(lat, lng, n, &
          & mapi, nxi, nyi, gi, &
          & maph, nxh, nyh, gh, &
          & Iobs, tau) bind(C, name = 'onepoint')
    use iso_c_binding
    integer(c_int), intent(in) :: n
    real(c_float), intent(in), dimension(n) :: lat, lng
    integer(c_int), intent(in) :: nxi, nyi
    real(c_float), dimension(nxi,nyi), intent(in) :: mapi
    integer(c_int), intent(in) :: nxh, nyh
    real(c_float), dimension(nxh,nyh), intent(in) :: maph
    real(c_float), dimension(6), intent(in) :: gi, gh
    real(c_float), dimension(n), intent(out) :: Iobs, tau
    type(modelpar) :: par
    type(geotransform) :: gti, gth
    call gti % import(gi)
    call gth % import(gh)
    call onepoint(par, lat, lng, mapi, gti, maph, gth, Iobs, tau)
  end subroutine

end module
