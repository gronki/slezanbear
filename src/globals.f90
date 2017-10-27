module globals

  implicit none

  integer, parameter :: sp = selected_real_kind(6,37)
  integer, parameter :: fp = selected_real_kind(15,307)

  real(fp), parameter :: radearth = 6371 * 1000
  real(fp), parameter :: pi = 4 * atan(real(1,fp))
  real(fp), parameter :: deg_in_rad = 180 / pi

  character(*), parameter :: outfn = 'lowersilesia'

  ! Wroclaw
  ! real(fp) :: llat = 50.80, llng = 16.55
  ! real(fp) :: ulat = 51.40, ulng = 17.40

  ! Izery
  ! real(fp) :: llat = 50.60, llng = 14.75
  ! real(fp) :: ulat = 51.20, ulng = 15.90

  ! Wroclaw center
  ! real(fp) :: llat = 51.00, llng = 16.87
  ! real(fp) :: ulat = 51.15, ulng = 17.09

  ! Lower Silesia
  real(fp) :: llat = 50.40, llng = 14.70
  real(fp) :: ulat = 51.60, ulng = 17.85

  real(fp) :: dmax = 180e3

  real(fp)           :: model_grid_meters   = 1000
  real(fp)           :: elev_grid_meters    = 50

  integer, parameter :: height_sect_num     = 96

  logical :: terrain_attenuation = .false.
  real(fp)           :: chkray_sect_meters  = 100
  real(fp)           :: chkray_min_dist     = 1500
  integer, parameter :: chkray_sect_num     = 2048

end module
