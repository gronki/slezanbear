module globals

  implicit none

  integer, parameter :: sp = selected_real_kind(6,37)
  integer, parameter :: dp = selected_real_kind(15,307)
  real, parameter :: radearth = 6371 * 1e3
  real, parameter :: pi = 4 * atan(real(1))
  real, parameter :: deg_in_rad = 180 / pi

end module
