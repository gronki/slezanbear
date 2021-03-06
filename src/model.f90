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
    real(dp) :: lat, lng
    real(dp) :: h, area
    real(dp) :: em
  end type source

contains

  !-----------------------------------------------------------------------------

  pure real(dp) function integrate(y, x) result(integral)
    real(dp), intent(in), dimension(:) :: x, y
    integral = sum((y(2:size(y)) + y(1:size(y)-1)) &
              &  * (x(2:size(x)) - x(1:size(x)-1))) / 2
  end function

  !-----------------------------------------------------------------------------

  subroutine regr(x, y, a, b)
    real(dp), dimension(:), intent(in) :: x, y
    real(dp) :: a, b, mx, my
    mx = sum(x) / size(x)
    my = sum(y) / size(y)
    a = sum((x - mx) * (y - my)) / sum((x - mx)**2)
    b = my - a * mx
  end subroutine

  !-----------------------------------------------------------------------------

  pure subroutine checkray(maph, gth, src, lat, lng, hsct, attn)
    ! altitude map and geotransform
    real, dimension(:,:), intent(in) :: maph
    type(geotransform), intent(in) :: gth
    ! two points and their relative heights (h / R)
    type(source), intent(in) :: src
    real(dp), intent(in) :: lat, lng, hsct(:)
    ! result: attenuation factor due to terrain (0 - obstruction, 1 - no obst)
    real(dp), intent(out) :: attn(:)
    ! internal variables
    integer :: nn
    real(dp) :: adist

    adist = angdist(src % lat, src % lng, lat, lng)

    if (adist * radearth < chkray_min_dist) then
      attn(:) = 1
      return
    end if

    nn = ceiling(radearth * adist / chkray_sect_meters)
    if (nn > chkray_sect_num) nn = chkray_sect_num
    if (nn < 3) nn = 3

    trace_ray: block

      real(dp), dimension(nn) :: t, lat_ray, lng_ray, h_ray, hterr
      real(dp), dimension(3) :: x1, x2, xi, xw, xr
      integer :: i, j

      forall (i = 1:nn) t(i) = f(real(i,dp) / (nn + 1), 2.0_dp)

      call geo2xyz(src % lat, src % lng, x1(1), x1(2), x1(3))
      call geo2xyz(      lat,       lng, x2(1), x2(2), x2(3))

      xi = cross(x1, x2)
      xi = xi / sqrt(sum(xi**2))
      xw = cross(xi, x1)

      interpolate_terrain: do i = 1, nn
        xi(:) = x1 * cos(t(i) * adist) + xw * sin(t(i) * adist)
        call xyz2geo(xi(1), xi(2), xi(3), lat_ray(i), lng_ray(i))
        call interpolmap(maph, gth, lat_ray(i), lng_ray(i), hterr(i))
      end do interpolate_terrain

      along_heights: do j = 1, size(hsct)

        if ( hsct(j) > radearth * adist ) then
          attn(j) = 1
          cycle along_heights
        end if

        calc_hray: do i = 1, nn
          xr(:) = x1 * (1 - t(i)) * (radearth + src % h) &
          &     + x2 * t(i)       * (radearth + hsct(j))
          h_ray(i) = sqrt(sum(xr**2)) - radearth
        end do calc_hray

        attn(j) = minval(h_ray - hterr) / chkray_toler_meters + 1
        if (attn(j) < 0) attn(j) = 0
        if (attn(j) > 1) attn(j) = 1

      end do along_heights
    end block trace_ray

  contains

    pure function cross(x,y) result(z)
      real(dp), dimension(3), intent(in) :: x,y
      real(dp), dimension(3) :: z
      z(1) =   x(2)*y(3) - x(3)*y(2)
      z(2) = - x(1)*y(3) + x(3)*y(1)
      z(3) =   x(1)*y(2) - x(2)*y(1)
    end function

    elemental real(dp) function f(t, k)
      real(dp), intent(in) :: t, k
      f = 2 * t + (k - 1) * t**2
      f = f / (k + 1)
    end function

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
      call arr2geo(gti, real(i,dp), real(j,dp), src(i,j) % lat, src(i,j) % lng)

      src(i,j) % area = abs(gti % det()) * (radearth / deg_in_rad)**2 &
            & * cos((src(i,j) % lat) / deg_in_rad)

      call interpolmap(maph, gth, src(i,j) % lat, src(i,j) % lng, src(i,j) % h)

      src(i,j) % em = mapi(i,j) * 1e-5
    end do
  end subroutine

  !-----------------------------------------------------------------------------

  pure subroutine onepoint(par, lat, lng, src, maph, gth, sky, hobs)
    type(modelpar), intent(in) :: par
    ! observer's latitude, longitude and
    real(dp), intent(in) :: lat, lng
    ! source map
    type(source), dimension(:,:), intent(in) :: src
    ! altitude map and geotransform for the map
    real, dimension(:,:), intent(in) :: maph
    type(geotransform), intent(in) :: gth
    real(dp), intent(out) :: sky(:), hobs

    type(modelpar) :: par1
    real(dp) :: tau

    integer i, j, k
    real(dp) :: JJP(height_sect_num)
    real(dp) :: hsct(height_sect_num), JJ(height_sect_num,4)

    ! get the altitude of the observer (in meters)
    call interpolmap(maph, gth, lat, lng, hobs)

    ! optical depth of the entire atmosphere
    tau = (par % alpha) * (par % relabs) * (par % hscale) &
    &                        * exp(-hobs / (par % hscale))
    ! background glow diminished by the absorption
    sky(:) = (par % skybg) ! * exp(-tau)

    par1 = par

    do concurrent (i = 1:height_sect_num)
      hsct(i) = genh(i, height_sect_num, hobs, 12 * par % hscale)
    end do

    ! iterate through all map points
    iter_src_rows: do j = 1,size(src,2)
     iter_src_cols: do i = 1,size(src,1)

        if (src(i,j) % em .eq. 0) cycle
        if (angdist(lat, lng, src(i,j) % lat, src(i,j) % lng) &
        &     .gt. (dmax / radearth)) cycle

        par1 % alpha = (par % alpha) * exp(-(src(i,j) % h) / (par % hscale))

        call kern(par1, src(i,j) % area, radearth + src(i,j) % h, &
            & (radearth + src(i,j) % h) &
            &   * angdist(lat, lng, src(i,j) % lat, src(i,j) % lng), &
            & hsct - src(i,j) % h, hobs - src(i,j) % h, &
            & JJ(:,1), JJ(:,2), JJ(:,3))

        ! najprostszy model: wszystko izotropowe + absorpcja
        JJP(:) = JJ(:,1) * JJ(:,3)
        sky(2) = sky(2) + (src(i,j) % em) * integrate(JJP, hsct)

        ! dodajemy krzywą światłości i funkcję fazową rozpraszania
        JJP(:) = JJP * JJ(:,2)
        sky(3) = sky(3) + (src(i,j) % em) * integrate(JJP, hsct)

        ! dodajemy przesłanianie terenem
        if ( terrain_attenuation .and. size(sky) .ge. 4 ) then
          call checkray(maph, gth, src(i,j), lat, lng, hsct, JJ(:,4))
          JJP(:) = JJP * JJ(:,4)
          sky(4) = sky(4) + (src(i,j) % em) * integrate(JJP, hsct)
        end if

      end do iter_src_cols
    end do iter_src_rows

    if (.not. terrain_attenuation) sky(4) = sky(3)

  contains

    elemental real(dp) function genh(i, n, h0, h1) result(h)
      integer, intent(in) :: i,n
      real(dp), intent(in) :: h0, h1
      real(dp) :: t,p
      real(dp), parameter :: k = 1.0
      t = real(i - 1, dp) / (n - 1)
      p = 1 - exp(-(h1 - h0) / (k * par % hscale))
      h = h0 - k * (par % hscale) * log(1 - t * p)
    end function

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine sbmap(par, mapi, gti, maph, gth, sky, hobs, gt)
    ! model parameters
    type(modelpar), intent(in) :: par
    ! coordinates of the point
    real(dp) :: lat, lng
    ! altitude and flux maps and respective geotransforms
    real(sp), dimension(:,:), intent(in) :: mapi, maph
    type(geotransform), intent(in) :: gti, gth, gt
    ! output
    real(dp), intent(out) :: hobs(:,:), sky(:,:,:)
    ! source map for computation
    type(source), dimension(:,:), allocatable :: src
    integer :: i,j

    allocate( src(size(mapi,1), size(mapi,2)) )
    call mksources(mapi, gti, maph, gth, src)

    write (*,'(2A10,A8,5A8)') 'LAT', 'LNG', 'HOBS', 'SKYBG', &
    &                        'SKY1', 'SKY2', 'SKY3'
    !$omp parallel do private(i,j,lat,lng)
    do j = 1,size(sky,2)
      do i = 1,size(sky,1)
        call arr2geo(gt, real(i,dp), real(j,dp), lat, lng)
        call onepoint(par, lat, lng, src, maph, gth, sky(i,j,:), hobs(i,j))
        write (*,'(2F10.4,F8.1,4F8.3)') lat, lng, hobs(i,j), ity2mag(sky(i,j,:))
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
          & sky, hobs, nx, ny, g) bind(C)
    use iso_c_binding
    integer(c_int), intent(in), value :: nxi, nyi
    real(c_float), dimension(nxi,nyi), intent(in) :: mapi
    integer(c_int), intent(in), value :: nxh, nyh
    real(c_float), dimension(nxh,nyh), intent(in) :: maph
    real(c_float), dimension(6), intent(in) :: gi, gh, g
    integer(c_int), intent(in), value :: nx, ny
    real(c_double), dimension(nx,ny,4), intent(out) :: sky
    real(c_double), dimension(nx,ny), intent(out) :: hobs
    type(modelpar) :: par
    type(geotransform) :: gti, gth, gt
    call gti % import(gi)
    call gth % import(gh)
    call gt % import(g)
    call sbmap(par, mapi, gti, maph, gth, sky, hobs, gt)
  end subroutine
end module
