module globals

  implicit none

  integer, parameter :: sp = selected_real_kind(6,37)
  integer, parameter :: fp = selected_real_kind(15,307)

  real(fp), parameter :: radearth = 6371 * 1000
  real(fp), parameter :: pi = 4 * atan(real(1,fp))
  real(fp), parameter :: deg_in_rad = 180 / pi

  real(fp) :: llat = 50.78, llng = 16.52
  real(fp) :: ulat = 51.50, ulng = 17.49
  real(fp) :: dmax = 150e3
  

end module
