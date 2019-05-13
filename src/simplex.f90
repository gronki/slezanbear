module simplex_m

  use globals
  implicit none
  private
  public :: simplex

contains

  subroutine simplex(x,xlo,xhi,y,f)
    interface
      pure real(dp) function simplexfun(x)
        import sp,dp
        real(dp), dimension(:), intent(in) :: x
      end function
    end interface
    real(dp), dimension(:,:), intent(inout) :: x
    real(dp), dimension(:), intent(in) :: xlo, xhi
    real(dp), dimension(:), intent(out) :: y
    procedure(simplexfun) :: f

    real(dp), dimension(*), parameter :: a = [ 1.0, 2.0, 0.5, 0.5 ]
    integer,  dimension(size(x,1) + 1) :: s
    real(dp), dimension(size(x,1)) :: x_0, x_r, x_e, x_c
    real(dp) :: y_r
    integer :: iter, i, n

    character(*), parameter :: actionfmt = '("simplex ", I2, ": ", &
          & A12, "(", I2, ")")'
    character(*), parameter :: statfmt = &
          & '(8X, "miny = ", Es10.3, " (", I2, ")", &
          &   3X, "maxy = ", Es10.3, " (", I2, ")")'

    n = size(x,1)

    if (size(x,2) /= n + 1) error stop "simplex: x must have shape (n,n+1)"
    if (size(y)   /= n + 1) error stop "simplex: size(y) /= size(x,1) + 1"
    if (size(xlo) /= n) error stop "simplex: size(xlo) /= size(x,1)"
    if (size(xhi) /= n) error stop "simplex: size(xhi) /= size(x,1)"

    !$omp parallel do
    do i = 1,n+1
      y(i) = f(x(:,i))
    end do
    !$omp end parallel do

    simplex_loop: do iter = 1, 72

      call isort(y,s)
      x_0 = sum(x(:,s(1:n)),2) / n

      write (0, statfmt) y(s(1)), s(1), y(s(n+1)), s(n+1)

      ! reflection
      x_r = x_0 + a(1) * (x_0 - x(:,s(n+1)))
      y_r = f(x_r)
      if (y_r < y(s(n))) then
        if (y_r > y(s(1))) then
          write (0,actionfmt) iter, "REFLECT", s(n+1)
          x(:,s(n+1)) = x_r
        else
          ! expansion
          x_e = x_0 + a(2) * (x_r - x_0)
          if (f(x_e) < y_r) then
            write (0,actionfmt) iter, "EXPAND", s(n+1)
            x(:,s(n+1)) = x_e
          else
            write (0,actionfmt) iter, "REFLECT", s(n+1)
            x(:,s(n+1)) = x_r
          end if
        end if
        y(s(n+1)) = f(x(:,s(n+1)))
      else
        ! contraction
        x_c = x_0 + a(3) * (x(:,s(n+1)) - x_0)
        if (f(x_c) < y(s(n+1))) then
          write (0,actionfmt) iter, "CONTRACT", s(n+1)
          x(:,s(n+1)) = x_c
          y(s(n+1)) = f(x(:,s(n+1)))
        else
          ! shrink
          write (0,actionfmt) iter, "SHRINK", s(1)
          !$omp parallel do
          do i = 2,n+1
            x(:,s(i)) = x(:,s(1)) + a(4) * (x(:,s(i)) - x(:,s(1)))
            y(s(i)) = f(x(:,s(i)))
          end do
          !$omp end parallel do
          ! do concurrent (i = 2:n+1)
          !   x(:,s(i)) = x(:,s(1)) + a(4) * (x(:,s(i)) - x(:,s(1)))
          !   y(s(i)) = f(x(:,s(i)))
          ! end do
        end if
      end if
    end do simplex_loop

    call isort(y,s)
    write (0, statfmt) y(s(1)), s(1), y(s(n+1)), s(n+1)

    x = x(:,s)
    y = y(s)

  end subroutine

  pure subroutine isort(x,s)
    real(dp), intent(in)  :: x(:)
    integer,  intent(out) :: s(:)
    integer :: i,j,sref

    if (size(x) /= size(s)) error stop "isort: size(x) /= size(s)"

    forall (i = 1:size(x)) s(i) = i

    do i = 2, size(x)
      sref = s(i)
      j = i - 1
      do while ( j > 0 .and. x(s(j)) > x(sref) )
        s(j + 1) = s(j)
        j = j - 1
      end do
      s(j + 1) = sref
    end do
  end subroutine

end module
