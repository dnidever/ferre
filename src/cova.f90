subroutine	cova(w,lambda_obs,obs,mobs,lsfarr,pf,e,cov)

!
!	Calculation of the covariance matrix (linear term) assuming 
!	Gaussian errors in the observations. See Eq. (14.4.11) of
!	the Numerical Recipes in fortran, and cov=[alpha]**-1. 
!	A solution needs to be determined first and stored in pf, 
!	which enters through the share data module.
!

use booklib, only	: mat_inv
use share, only		: dp,n_p,ndim,nov,indv,npix,nlambda1,badflux,winter,lambda_syn, & 
                                    mlsf,nlsf
                                    
use grid_geometry
implicit none

!in/out
real(dp), intent(in)       :: w(nlambda1)           ! weights
real(dp), intent(in)       :: lambda_obs(nlambda1)  ! wavelengths for observations
real(dp), intent(in)       :: obs(nlambda1)  ! observations
real(dp), intent(in)       :: mobs                 ! mean or median of obs array
real(dp), intent(in)       :: lsfarr(mlsf,nlsf)    ! lsfarray
real(dp), intent(in)       :: pf(ndim)          ! vector of fixed parameters
real(dp), intent(in)  	   :: e(nlambda1)		!errors in the observations
real(dp), intent(out)  	   :: cov(nov,nov)	!covariance matrix

!locals
real(dp), parameter	   :: delta=1.1					!factor to compute derivatives
integer		  	   :: i,j,k							!counters
real(dp)		   :: p(ndim)	        	   		!parameters
real(dp)		   :: flux(nlambda1)				!fluxes
real(dp)		   :: flux1(nlambda1)				!fluxes
real(dp)		   :: flux2(nlambda1)				!fluxes
real(dp)		   :: par1,par2	        			!tmp vars. for partial derivat.
real(dp)		   :: partial(nov,nlambda1)				!partial derivatives 
real(dp)		   :: alpha(nov,nov)				!curvature matrix
integer			   :: error = 0						!error code from mat_inv

if (sparse_active) then
    call sparse_covariance()
    return
endif
p=pf
do i=1,nov

	call flx(pf,lambda_obs,obs,e,mobs,lsfarr,flux)
	
	k=indv(i)
	
	!initialize
	flux1(1:nlambda1)=badflux
	flux2(1:nlambda1)=badflux
		
	!two estimates of the partials, increasing
	p(k)=pf(k)*delta
	if (p(k) <= 1.0_dp) then
		call flx(p,lambda_obs,obs,e,mobs,lsfarr,flux1)
	endif

	!and decreasing the parameters
	p(k)=pf(k)/delta
	if (p(k) >= 0.0_dp) then
		call flx(p,lambda_obs,obs,e,mobs,lsfarr,flux2)
	endif
	
		
	do j=1,nlambda1
	
	   par1=(flux1(j)-flux(j))/pf(k)/(delta-1._dp)
	   par2=(flux2(j)-flux(j))/pf(k)/(1._dp/delta-1._dp)

	   if (flux1(j) >= 0._dp .and. flux2(j) >=0._dp) then
	      	partial(i,j)=(par1+par2)/2._dp
	   else
		partial(i,j)=0.0_dp
		if (flux1(j) >= 0.0_dp)	partial(i,j)=par1
		if (flux2(j) >= 0.0_dp) partial(i,j)=par2
		
	   endif
	   
	enddo
enddo

do i=1,nov
	do j=1,nov
		alpha(j,i)=sum(partial(i,1:nlambda1)*partial(j,1:nlambda1)/e(1:nlambda1)**2)
	enddo
enddo

!booklib lib
call mat_inv(alpha,cov,nov,nov,error)
!Pang's lib
!call migs(alpha,nov,cov,flux1)

if (error /= 0) then
	write(*,*) 'cova: ERROR'
	write(*,*) 'mat_inv returned an non-zero error=',error
	write(*,*) 'probably indicating a singular matrix'
	cov(:,:)=-1._dp
endif

contains

subroutine sparse_covariance()
integer :: axis, interval, side, d, n
real(dp) :: span, x, h, left_model(nlambda1), right_model(nlambda1)
cov=-1._dp
call flx(pf,lambda_obs,obs,e,mobs,lsfarr,flux)
if (invalid_model(flux)) return
do d=1,nov
    axis=indv(d)
    n=n_p(axis)
    span=axis_values(n,axis)-axis_values(1,axis)
    x=axis_values(1,axis)+pf(axis)*span
    interval=bracket(axis,x)
    h=0.001_dp*(axis_values(interval+1,axis)-axis_values(interval,axis))/span
    if (pf(axis)<h .or. pf(axis)>1._dp-h) return
    p=pf
    p(axis)=pf(axis)-h
    call flx(p,lambda_obs,obs,e,mobs,lsfarr,left_model)
    p(axis)=pf(axis)+h
    call flx(p,lambda_obs,obs,e,mobs,lsfarr,right_model)
    if (invalid_model(left_model) .or. invalid_model(right_model)) return
    partial(d,:)=(right_model-left_model)/(2._dp*h)
enddo
do d=1,nov
    do side=1,nov
        alpha(side,d)=sum(partial(d,:)*partial(side,:)/e**2)
    enddo
enddo
call mat_inv(alpha,cov,nov,nov,error)
if (error/=0) cov=-1._dp
end subroutine

end subroutine cova



