module globals

  implicit none

  integer, parameter :: sp = selected_real_kind(6,37)
  integer, parameter :: fp = selected_real_kind(15,307)

  real(fp), parameter :: radearth = 6371 * 1000
  real(fp), parameter :: pi = 4 * atan(real(1,fp))
  real(fp), parameter :: deg_in_rad = 180 / pi

end module
