subroutine normal2(pf,ocov)
!from physical to normalized ([0-1]) units

use share, only: dp, ndim, n_p, llimits, steps, ntie, indtie
 
use grid_geometry, only: sparse_active, grid_upper
implicit none

!params
real(dp), intent(inout) :: pf(ndim) ! full vector of parameters
real(dp), intent(inout) :: ocov(ndim,ndim)  !full covariance matrix

!locals
real(dp)                :: ulimiti,ulimitj
integer                 :: i,j

do i=1,ndim
  ulimiti=grid_upper(i)
  do j=1,ndim
    ulimitj=grid_upper(j)
    if (sparse_active .and. (n_p(i)==1 .or. n_p(j)==1)) then
        ocov(j,i)=0._dp
    else
        ocov(j,i)=ocov(j,i)/(ulimiti-llimits(i))/(ulimitj-llimits(j))
    endif
  enddo
enddo

do i=1,ndim
  ulimiti=grid_upper(i)
  !pf(i)=0.5_dp+(2._dp*pf(i)-llimits(i)-ulimiti)/(2._dp*(ulimit-llimits(i)))
  if (sparse_active .and. n_p(i)==1) then
      if (pf(i)/=llimits(i)) then
          pf(i)=2._dp
      else
          pf(i)=0._dp
      endif
  else
      pf(i)=(pf(i)-llimits(i)) / (ulimiti-llimits(i))
  endif
enddo

if (indtie(1) > 0 .and. .not.sparse_active) then 
  do i=1,ntie
    if (pf(indtie(i)) < 0.0_dp) then
      pf(indtie(i))=0.0_dp
    else
        if (pf(indtie(i)) >= 1.0_dp) pf(indtie(i))=0.999999_dp
    endif
  enddo
endif

end subroutine normal2

