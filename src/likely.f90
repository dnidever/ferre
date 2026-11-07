subroutine likely(wr,fname,w,chiscale,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,pos)

!
!	Average estimates
!	pos(i)=sum pos*exp(-chi/min_chi*10) through the whole grid (no interpolation)
!
!	If wr=1 then the chi**2 surface is printed out on fname+'.chi'
!

use share, only: 	dp,flen,liobuffer,n_p,ndim,npix,nlambda1,nov,indv, &
                            llimits,steps,aa,mlsf,nlsf
use fun, only: 	objfun

use grid_geometry
implicit none

!input/output
integer, intent(in)	::	wr			!write chi**2 file?
character(len=flen),intent(in):: 	fname			   !spectrum id
real(dp), intent(in)    :: w(nlambda1)          ! weights
real(dp), intent(in)    :: chiscale		!chi**2/sum(w_i(f_i-obs_i))^2)
real(dp), intent(inout) :: pf(ndim)             ! vector of fixed parameters
real(dp), intent(in)    :: pf0(ndim) ! pars read from pfile (in physical units)
real(dp), intent(in)    :: obs(nlambda1)        ! vector of observations
real(dp), intent(in)    :: lambda_obs(nlambda1) ! vector of wavelengths for observations
real(dp), intent(in)    :: e_obs(nlambda1) ! vector of uncertainties for observations
real(dp), intent(in)    :: mobs             ! mean or median of obs array
real(dp), intent(in)    :: lsfarr(mlsf,nlsf)    ! lsfarray
real(dp), intent(out)	:: pos(nov)		   !expected values

!locals
integer			::	i,j,k
integer			::	nt			!number of discrete pts
real(dp)		::	chi
real(dp)		::	min_chi
real(dp)		::	p(nov)			!local vector of pars.
real(dp)		::	ph(nov)			!local vect of phys. par
real(dp)		::	den(nov)		!sum exp(-chi/min_chi*10)
real(dp)		::	ulimit			!temp variable
integer			::	stlen			!tmp counter
character(len=flen)	:: 	chifile			!output file 


if (wr > 0) then 
 	stlen=len_trim(fname)
	chifile(1:stlen)=fname(1:stlen)
	chifile(stlen+1:stlen+4)='.chi'
	write(*,*)fname
	write(*,*)chifile
	open(5,file=chifile,status='unknown',recl=liobuffer)
endif

if (sparse_active) then
    call sparse_average()
    if (wr > 0) close(5)
    return
endif

!calculate how many discrete points at which to evaluate fitness
nt=n_p(indv(nov))
do j=2,nov
	nt=nt*n_p(indv(nov-j+1))
enddo


!initialize pos and min_chi at center
do j=1,nov
  p(j)=0.5_dp
  pos(j)=p(j)
enddo
call objfun(w,chiscale,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,p,min_chi)

!run over the grid and find the absolute minimum location
do i=1,nt
   do j=1,nov
   	p(j)=(real(aa(j,i))-1._dp)/(real(n_p(indv(j)))-1._dp) !nodes 
   	ulimit=llimits(indv(j))+steps(indv(j))*(n_p(indv(j))-1)
   	ph(j)=p(j)*(ulimit-llimits(indv(j)))+llimits(indv(j)) !physical at nodes
   enddo
   call objfun(w,chiscale,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,p,chi)
   if (wr > 0) write(5,'(1x,21(1x,f8.5))')ph,chi
   if (chi < min_chi) min_chi=chi
enddo

!initialize pos and den 
do j=1,nov
  p(j)=0.0_dp
  pos(j)=0.0_dp
  den(j)=0.0_dp
enddo

!run over the grid and calculate expected values
do i=1,nt
   do j=1,nov
   	p(j)=(real(aa(j,i))-1._dp)/(real(n_p(indv(j)))-1._dp) !nodes 
   enddo
   call objfun(w,chiscale,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,p,chi)
   chi=chi/min_chi
   do j=1,nov
   	pos(j)=pos(j)+p(j)*exp(-chi*10._dp)
   	den(j)=den(j)+exp(-chi*10._dp)
   enddo   
enddo

do j=1,nov
   pos(j)=pos(j)/den(j)
enddo  


if (wr > 0) close(5)

contains

subroutine sparse_average()
integer :: record, axis
real(dp) :: total_weight, weight, scale
min_chi=huge(1._dp)
pos=0.5_dp
do record=1,size(node_coordinates,2)
    do j=1,nov
        axis=indv(j)
        p(j)=(node_coordinates(axis,record)-llimits(axis))/(grid_upper(axis)-llimits(axis))
    enddo
    call objfun(w,chiscale,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,p,chi)
    if (chi<fit_penalty/50._dp) min_chi=min(min_chi,chi)
    if (wr>0) write(5,'(1x,21(1x,f8.5))') node_coordinates(indv(1:nov),record),chi
enddo
if (min_chi==huge(1._dp)) return
scale=max(min_chi,epsilon(1._dp))
pos=0._dp
total_weight=0._dp
do record=1,size(node_coordinates,2)
    do j=1,nov
        axis=indv(j)
        p(j)=(node_coordinates(axis,record)-llimits(axis))/(grid_upper(axis)-llimits(axis))
    enddo
    call objfun(w,chiscale,pf,pf0,obs,lambda_obs,e_obs,mobs,lsfarr,p,chi)
    if (chi>=fit_penalty/50._dp) cycle
    weight=exp(-10._dp*(chi-min_chi)/scale)
    pos=pos+weight*p
    total_weight=total_weight+weight
enddo
pos=pos/total_weight
end subroutine

end subroutine likely
