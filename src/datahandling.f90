module datahandling

  use globals

  implicit none

  type dataentry
    real(dp) :: lat, lng
    real(dp) :: sky, dsky
  end type

  integer, parameter, private :: ndata_max = 1024

contains

  subroutine read_data(fn, dat)
    character(len = *), intent(in) :: fn
    type(dataentry), intent(inout), allocatable :: dat(:)
    type(dataentry) :: dat_all(ndata_max)
    integer :: un, errno, ndata, i

    open(newunit = un, file = fn, action = 'read')

    ndata = 0
    do i = 1, size(dat_all)
      read(un, *, iostat = errno) dat_all(i) % lat, dat_all(i) % lng, dat_all(i) % sky
      if (errno .ne. 0) exit
      ndata = ndata + 1
      if (ndata >= ndata_max) error stop 'increase ndata_max in datahandling.f90!'
    end do

    close(un)

    dat = dat_all(1:ndata)
  end subroutine

end module
