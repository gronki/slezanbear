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

contains

  !-----------------------------------------------------------------------------

  pure real(fp) function integrate(y, x) result(integral)
    real(fp), intent(in), dimension(:) :: x, y
    integral = sum((y(2:size(y)) + y(1:size(y)-1)) &
              &  * (x(2:size(x)) - x(1:size(x)-1))) / 2
  end function

  !-----------------------------------------------------------------------------

  pure subroutine checkray(maph, gth, src, lat, lng, hsct, attn)
    ! altitude map and geotransform
    real, dimension(:,:), intent(in) :: maph
    type(geotransform), intent(in) :: gth
    ! two points and their relative heights (h / R)
    type(source), intent(in) :: src
    real(fp), intent(in) :: lat, lng, hsct(:)
    ! result: attenuation factor due to terrain (0 - obstruction, 1 - no obst)
    real(fp), intent(out) :: attn(:)
    ! internal variables
    integer :: nn
    real(fp) :: adist

    adist = angdist(src % lat, src % lng, lat, lng)

    nn = ceiling(radearth * adist / chkray_sect_meters)
    if (nn > chkray_sect_num) nn = chkray_sect_num
    if (nn < 3) nn = 3

    trace_ray: block

      real(fp), dimension(nn) :: t, lat_ray, lng_ray, h_ray, hterr, Q
      real(fp), dimension(3) :: x1, x2, xi, xw, xr
      integer :: i, j

      forall (i = 1:nn) t(i) = real(i,fp) / (nn + 1)

      call geo2xyz(src % lat, src % lng, x1(1), x1(2), x1(3))
      call geo2xyz(      lat,       lng, x2(1), x2(2), x2(3))

      call cross(x1, x2, xi)
      call cross(xi / sqrt(sum(xi**2)), x1, xw)

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
          xr(:) = x1 * (1 - t(i)) * (radearth + src % h ) &
          &     + x2 * t(i)       * (radearth + hsct(j) )
          h_ray(i) = sqrt(sum(xr**2)) - radearth
        end do calc_hray

        Q(:) = (h_ray - hterr) / (adist * radearth * t)
        attn(j) = th(minval(Q) / real(0.01,fp) + 0.5)
      end do along_heights

    end block trace_ray

  contains

    pure subroutine cross(x,y,z)
      real(fp), dimension(3), intent(in) :: x,y
      real(fp), dimension(3), intent(out) :: z
      z(1) =   x(2)*y(3) - x(3)*y(2)
      z(2) = - x(1)*y(3) + x(3)*y(1)
      z(3) =   x(1)*y(2) - x(2)*y(1)
    end subroutine

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

      src(i,j) % em = mapi(i,j) * 1e-5
    end do
  end subroutine

  !-----------------------------------------------------------------------------

  pure subroutine onepoint(par, lat, lng, src, maph, gth, I1, I2, hobs)
    type(modelpar), intent(in) :: par
    ! observer's latitude, longitude and
    real(fp), intent(in) :: lat, lng
    ! source map
    type(source), dimension(:,:), intent(in) :: src
    ! altitude map and geotransform for the map
    real, dimension(:,:), intent(in) :: maph
    type(geotransform), intent(in) :: gth
    real(fp), intent(out) :: I1, I2, hobs

    type(modelpar) :: par1
    real(fp) :: tau

    integer i, j, k
    real(fp) :: JJP(height_sect_num)
    real(fp) :: hsct(height_sect_num), JJ(height_sect_num,4)

    ! get the altitude of the observer (in meters)
    call interpolmap(maph, gth, lat, lng, hobs)

    ! optical depth of the entire atmosphere
    tau = (par % alpha) * (par % relabs) * exp(-hobs / (par % hscale))
    ! background glow diminished by the absorption
    I1 = (par % skybg) * exp(-tau)
    I2 = (par % skybg) * exp(-tau)

    par1 = par

    forall (i = 1:height_sect_num)
      hsct(i) = genh(i, height_sect_num, hobs, 7 * par % hscale)
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

      JJP(:) = JJ(:,1) * JJ(:,2) * JJ(:,3)
      I1 = I1 + (src(i,j) % em) * integrate(JJP, hsct)

      if ( terrain_attenuation ) then
        call checkray(maph, gth, src(i,j), lat, lng, hsct, JJ(:,4))
        JJP(:) = JJP * JJ(:,4)
        I2 = I2 + (src(i,j) % em) * integrate(JJP, hsct)
      end if

    end do

    if (.not. terrain_attenuation) I2 = I1

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

  subroutine sbmap(par, mapi, gti, maph, gth, I1, I2, hobs, gt)
    ! model parameters
    type(modelpar), intent(in) :: par
    ! coordinates of the point
    real(fp) :: lat, lng
    ! altitude and flux maps and respective geotransforms
    real(sp), dimension(:,:), intent(in) :: mapi, maph
    type(geotransform), intent(in) :: gti, gth, gt
    ! output
    real(fp), intent(out), dimension(:,:) :: I1, I2, hobs
    ! source map for computation
    type(source), dimension(:,:), allocatable :: src
    integer :: i,j

    allocate( src(size(mapi,1), size(mapi,2)) )
    call mksources(mapi, gti, maph, gth, src)

    !$omp parallel do private(j,lat,lng)
    do i = 1,size(I1,1)
      do j = 1,size(I1,2)
        call arr2geo(gt, real(i,fp), real(j,fp), lat, lng)
        call onepoint(par, lat, lng, src, maph, gth, &
        &             I1(i,j), I2(i,j), hobs(i,j))
        write (*,'(I4,1X,I4)')  i, j
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
          & I1, I2, hobs, nx, ny, g) bind(C)
    use iso_c_binding
    integer(c_int), intent(in), value :: nxi, nyi
    real(c_float), dimension(nxi,nyi), intent(in) :: mapi
    integer(c_int), intent(in), value :: nxh, nyh
    real(c_float), dimension(nxh,nyh), intent(in) :: maph
    real(c_float), dimension(6), intent(in) :: gi, gh, g
    integer(c_int), intent(in), value :: nx, ny
    real(c_double), dimension(nx,ny), intent(out) :: I1,I2,hobs
    type(modelpar) :: par
    type(geotransform) :: gti, gth, gt
    call gti % import(gi)
    call gth % import(gh)
    call gt % import(g)
    call sbmap(par, mapi, gti, maph, gth, I1, I2, hobs, gt)
  end subroutine
end module
