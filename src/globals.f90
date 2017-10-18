module globals

  implicit none

  integer, parameter :: sp = selected_real_kind(6,37)
  integer, parameter :: fp = selected_real_kind(15,307)

  real(fp), parameter :: radearth = 6371 * 1000
  real(fp), parameter :: pi = 4 * atan(real(1,fp))
  real(fp), parameter :: deg_in_rad = 180 / pi

  character(*), parameter :: outfn = 'model4.tif'

  real(fp) :: llat = 50.70, llng = 16.50
  real(fp) :: ulat = 51.50, ulng = 17.50
  real(fp) :: dmax = 200e3

  real(fp)           :: model_grid_meters   = 500
  real(fp)           :: elev_grid_meters    = 125

  integer, parameter :: height_sect_num     = 80

  logical :: terrain_attenuation = .false.
  real(fp)           :: chkray_sect_meters  = 125
  real(fp)           :: chkray_maxangle     = 30
  integer, parameter :: chkray_sect_num     = 1024

end module
