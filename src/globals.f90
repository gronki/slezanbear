module globals

  implicit none

  integer, parameter :: sp = selected_real_kind(6,37)
  integer, parameter :: fp = selected_real_kind(15,307)

  real(fp), parameter :: radearth = 6371 * 1000
  real(fp), parameter :: pi = 4 * atan(real(1,fp))
  real(fp), parameter :: deg_in_rad = 180 / pi

  character(*), parameter :: outfn = 'model_sil_hires.tif'

  real(fp), parameter :: magconst = 5.5

  ! real(fp) :: llat = 50.60, llng = 14.80
  ! real(fp) :: ulat = 51.20, ulng = 15.90
  real(fp) :: llat = 50.50, llng = 14.80
  real(fp) :: ulat = 51.50, ulng = 17.70

  real(fp) :: dmax = 200e3

  real(fp)           :: model_grid_meters   = 250
  real(fp)           :: elev_grid_meters    = 125
  real(fp)           :: elev_grid_meters    = 50

  integer, parameter :: height_sect_num     = 90

  logical :: terrain_attenuation = .true.
  real(fp)           :: chkray_sect_meters  = 250
  real(fp)           :: chkray_sect_meters  = 100
  integer, parameter :: chkray_sect_num     = 2048

end module
