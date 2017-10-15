module gdal_test_common

  use geo
  implicit none

contains

  subroutine checkdata(gt, dat, i, j)
    real(8), intent(in) :: dat(:,:,:)
    integer, intent(in) :: i,j
    type(geotransform) :: gt
    real(8) :: lat, lng
    call arr2geo(gt, real(i,8), real(j,8), lat, lng)
    write (*,'(I5,I5,1X,4F9.5,Es11.4)') i, j, lat, lng, dat(i,j,1), dat(i,j,2), dat(i,j,3)
  end subroutine

end module
