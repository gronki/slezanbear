module slezanbear

  use globals
  use geo
  use kernel

  implicit none

  type lpmap
    real(sp) :: lat, lng, h
    real(sp) :: iz, tau
  end type lpmap

  type source
    real(fp) :: lat, lng
    real(fp) :: h, area
    real(fp) :: em
  end type source

  logical :: terrain_attenuation = .true.

contains

  !-----------------------------------------------------------------------------

  pure real(fp) function integrate(y, x) result(integral)
    real(fp), intent(in), dimension(:) :: x, y
    integral = sum((y(2:size(y)) + y(1:size(y)-1)) &
              &  * (x(2:size(x)) - x(1:size(x)-1))) / 2
  end function

  !-----------------------------------------------------------------------------

  pure subroutine checkray(maph, gth, lat1, lng1, h1, lat2, lng2, h2, attn)
    ! altitude map and geotransform
    real, dimension(:,:), intent(in) :: maph
    type(geotransform), intent(in) :: gth
    ! two points and their relative heights (h / R)
    real(fp), intent(in) :: lat1, lng1, h1
    real(fp), intent(in) :: lat2, lng2, h2
    ! result: attenuation factor due to terrain (0 - obstruction, 1 - no obst)
    real(fp), intent(out) :: attn
    ! internal variables
    integer, parameter :: nnmax = 2**10
    integer :: nn,i
    real(fp) :: raylength
    real(fp), dimension(nnmax) :: t, lat_ray, lng_ray, h_ray, hterr, Q

    raylength = raydist(lat1, lng1, h1 / radearth, &
      & lat2, lng2, h2 / radearth) * radearth

    if ( abs(h1 - h2) / raylength > 0.707 ) then
      attn = 1
      return
    end if

    nn = nint(raylength / 250)
    if (nn > nnmax) nn = nnmax
    if (nn < 3) nn = 3

    forall (i = 1:nn) t(i) = real(i,fp) / nn

    call georay(lat1, lng1, h1 / radearth, lat2, lng2, h2 / radearth, &
      & t(1:nn), lat_ray(1:nn), lng_ray(1:nn), h_ray(1:nn))

    do concurrent (i = 1:nn)
      call interpolmap(maph, gth, lat_ray(i), lng_ray(i), hterr(i))
    end do

    Q(1:nn) = (h_ray(1:nn) * radearth - hterr(1:nn)) &
          & / (raylength * t(1:nn))
    attn = th(minval(Q(1:nn)) / real(0.001,fp))
  end subroutine checkray

  !-----------------------------------------------------------------------------

  pure subroutine mksources(mapi, gti, maph, gth, src)
    ! altitude and flux maps and respective geotransforms
    real, dimension(:,:), intent(in) :: mapi, maph
    type(geotransform), intent(in) :: gti, gth
    ! source map for computation
    type(source), dimension(:,:), intent(out) :: src
    integer i,j
    ! interpolate altitude for all source points and compute their coordinates
    do concurrent (i = 1:size(mapi,1), j = 1:size(mapi,2))
      call arr2geo(gti, real(i,fp), real(j,fp), src(i,j) % lat, src(i,j) % lng)

      src(i,j) % area = abs(gti % det()) * (radearth / deg_in_rad)**2 &
            & * cos((src(i,j) % lat) / deg_in_rad)

      call interpolmap(maph, gth, src(i,j) % lat, src(i,j) % lng, src(i,j) % h)

      src(i,j) % em = mapi(i,j)
    end do
  end subroutine

  !-----------------------------------------------------------------------------

  pure subroutine onepoint(par, lat, lng, src, maph, gth, I1, I2)
    type(modelpar), intent(in) :: par
    ! observer's latitude, longitude and
    real(fp), intent(in) :: lat, lng
    ! source map
    type(source), dimension(:,:), intent(in) :: src
    ! altitude map and geotransform for the map
    real, dimension(:,:), intent(in) :: maph
    type(geotransform), intent(in) :: gth
    real(fp), intent(out) :: I1, I2

    type(modelpar) :: par1
    real(fp) :: tau, hobs

    integer i, j, k
    integer, parameter :: nh = 2**6
    real(fp) :: JJP(nh)
    real(fp) :: hsct(nh), JJ(nh,4)

    ! get the altitude of the observer (in meters)
    call interpolmap(maph, gth, lat, lng, hobs)

    ! optical depth of the entire atmosphere
    tau = (par % alpha) * (par % relabs) * exp(-hobs / (par % hscale))
    ! background glow diminished by the absorption
    I1 = (par % skybg) * exp(-tau)
    I2 = (par % skybg) * exp(-tau)

    par1 = par

    forall (i = 1:nh)
      hsct(i) = genh(i, nh, hobs, 7 * par % hscale)
    end forall

    ! iterate through all map points
    do concurrent (i = 1:size(src,1), j = 1:size(src,2), &
          & (src(i,j) % em .ne. 0) &
          & .and. (angdist(lat, lng, src(i,j) % lat, src(i,j) % lng) &
          &     .le. (dmax / radearth)))

      par1 % alpha = (par % alpha) * exp(-(src(i,j) % h) / (par % hscale))

      call kern(par1, src(i,j) % area, radearth + src(i,j) % h, &
          & (radearth + src(i,j) % h) &
          &   * angdist(lat, lng, src(i,j) % lat, src(i,j) % lng), &
          & hsct - src(i,j) % h, hobs - src(i,j) % h, &
          & JJ(:,1), JJ(:,2), JJ(:,3))

      check_lines_all_heights: do k = 1,nh
        call checkray(maph, gth,                                &
              & src(i,j) % lat, src(i,j) % lng, src(i,j) % h,   &
              & lat,            lng,            hsct(k),        &
              & JJ(k,4))
      end do check_lines_all_heights

      JJP(:) = JJ(:,1) * JJ(:,2) * JJ(:,3)
      I1 = I1 + INTEGRATE(JJP,hsct) * (src(i,j) % em) * 1e-9

      JJP(:) = JJ(:,1) * JJ(:,2) * JJ(:,3) * JJ(:,4)
      I2 = I2 + INTEGRATE(JJP,hsct) * (src(i,j) % em) * 1e-9

    end do

  contains

    elemental real(fp) function genh(i, n, h0, h1) result(h)
      integer, intent(in) :: i,n
      real(fp), intent(in) :: h0, h1
      real(fp) :: t,p
      real(fp), parameter :: k = 1.5
      t = real(i - 1, fp) / (n - 1)
      p = 1 - exp(-(h1 - h0) / (k * par % hscale))
      h = h0 - k * (par % hscale) * log(1 - t * p)
    end function

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine sbmap(par, mapi, gti, maph, gth, I1, I2, gt)
    ! model parameters
    type(modelpar), intent(in) :: par
    ! coordinates of the point
    real(fp) :: lat, lng
    ! altitude and flux maps and respective geotransforms
    real(sp), dimension(:,:), intent(in) :: mapi, maph
    type(geotransform), intent(in) :: gti, gth, gt
    ! output
    real(fp), intent(out), dimension(:,:) :: I1, I2
    ! source map for computation
    type(source), dimension(:,:), allocatable :: src
    integer :: i,j

    allocate( src(size(mapi,1), size(mapi,2)) )
    call mksources(mapi, gti, maph, gth, src)

    !$omp parallel do private(j,lat,lng)
    do i = 1,size(I1,1)
      do j = 1,size(I1,2)
        call arr2geo(gt, real(i,fp), real(j,fp), lat, lng)
        write (*,'(I4,1X,I4)')  i, j
        call onepoint(par, lat, lng, src, maph, gth, I1(i,j), I2(i,j))
      end do
    end do
    !$omp end parallel do

    deallocate( src )

  end subroutine

end module

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

module c_interfaces

  use slezanbear
  implicit none

contains

  subroutine sbmap_c(  &
          & mapi, nxi, nyi, gi, &
          & maph, nxh, nyh, gh, &
          & I1, I2, nx, ny, g) bind(C)
    use iso_c_binding
    integer(c_int), intent(in), value :: nxi, nyi
    real(c_float), dimension(nxi,nyi), intent(in) :: mapi
    integer(c_int), intent(in), value :: nxh, nyh
    real(c_float), dimension(nxh,nyh), intent(in) :: maph
    real(c_float), dimension(6), intent(in) :: gi, gh, g
    integer(c_int), intent(in), value :: nx, ny
    real(c_double), dimension(nx,ny), intent(out) :: I1,I2
    type(modelpar) :: par
    type(geotransform) :: gti, gth, gt
    call gti % import(gi)
    call gth % import(gh)
    call gt % import(g)
    call sbmap(par, mapi, gti, maph, gth, I1, I2, gt)
  end subroutine
end module
