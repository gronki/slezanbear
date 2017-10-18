module globals

  implicit none

  integer, parameter :: sp = selected_real_kind(6,37)
  integer, parameter :: fp = selected_real_kind(15,307)

  real(fp), parameter :: radearth = 6371 * 1000
  real(fp), parameter :: pi = 4 * atan(real(1,fp))
  real(fp), parameter :: deg_in_rad = 180 / pi

  real(fp) :: llat = 50.78, llng = 16.52
  real(fp) :: ulat = 51.50, ulng = 17.49
  real(fp) :: dmax = 100e3

  real(fp)           :: model_grid_meters   = 1000
  real(fp)           :: elev_grid_meters    = 125

  integer, parameter :: height_sect_num     = 60

  real(fp)           :: chkray_sect_meters  = 125
  real(fp)           :: chkray_maxangle     = 30
  integer, parameter :: chkray_sect_num     = 2**10


end module
