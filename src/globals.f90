module globals

  implicit none

  integer, parameter :: sp = selected_real_kind(6,37)
  integer, parameter :: dp = selected_real_kind(15,307)

  real(dp), parameter :: radearth = 6371 * 1000
  real(dp), parameter :: pi = 4 * atan(real(1,dp))
  real(dp), parameter :: deg_in_rad = 180 / pi

  character(*), parameter :: outfn = 'lowersilesia'

  ! Wroclaw
  ! real(dp) :: llat = 50.80, llng = 16.55
  ! real(dp) :: ulat = 51.40, ulng = 17.40

  ! Izery
  ! real(dp) :: llat = 50.60, llng = 14.75
  ! real(dp) :: ulat = 51.20, ulng = 15.90

  ! Wroclaw center
  ! real(dp) :: llat = 51.00, llng = 16.87
  ! real(dp) :: ulat = 51.15, ulng = 17.09

  ! Lower Silesia
  real(dp) :: llat = 50.40, llng = 14.70
  real(dp) :: ulat = 51.60, ulng = 17.85

  real(dp) :: dmax = 180e3

  real(dp)           :: model_grid_meters   = 2000
  real(dp)           :: elev_grid_meters    = 125

  integer, parameter :: height_sect_num     = 96

  logical :: terrain_attenuation = .true.
  real(dp)           :: chkray_sect_meters  = 250
  real(dp)           :: chkray_min_dist     = 500
  integer, parameter :: chkray_sect_num     = 2048
  real(dp)           :: chkray_toler_meters = 2

end module
