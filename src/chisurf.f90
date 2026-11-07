
subroutine chisurf(w,obs,lambda_obs,fname,e_obs,mobs,lsfarr)

!
!	Output the chi**2 surface
!
!	

use share, only: dp,flen,liobuffer,n_p,ntot,ndim,nov,npix,nlambda1,f,winter,lambda_syn,mlsf,nlsf
use grid_geometry
implicit none

!input/output
real(dp), intent(in)        :: w(nlambda1)       ! weights
real(dp), intent(in)        :: obs(nlambda1)     ! vector of observations
real(dp), intent(in)        :: lambda_obs(nlambda1) ! vector of wavelengths for observations
character(len=flen),intent(in):: 	fname		 !spectrum id

real(dp), intent(in) :: e_obs(nlambda1),mobs,lsfarr(mlsf,nlsf)
!locals
integer			::	i,stlen
real(dp)		::	lchi,flux(nlambda1),p(ndim)
character(len=flen)	:: 	chifile		!output file with chi**2 surf

!write(*,*)'in chisurf'


stlen=len_trim(fname)
chifile(1:stlen)=fname(1:stlen)
chifile(stlen+1:stlen+4)='.chi'
chifile(stlen+5:flen)=' '
chifile=trim(chifile)
!write(*,*)fname
write(*,*)'printing out ',chifile
open(15,file=chifile,status='unknown',recl=liobuffer)
!write(15,*)fname
!write(5,*)obs
if (sparse_active) then
    do i=1,ntot
        p=node_coordinates(:,i)
        call normal(p)
        call flx(p,lambda_obs,obs,e_obs,mobs,lsfarr,flux)
        if (invalid_model(flux)) then
            lchi=log10(fit_penalty)
        else
            lchi=log10(max(sum(w*(obs-flux)**2),tiny(1._dp)))
        endif
        write(15,*) node_coordinates(:,i),lchi
    enddo
    close(15)
    return
endif
do i=1,ntot
	lchi=0.0_dp
	if (winter == 2) then 
		call wresample(lambda_syn,f(1:npix,i),npix,lambda_obs,flux,nlambda1)
		lchi=sum((obs(1:nlambda1)-flux(1:nlambda1))**2)
	else
		lchi=sum((obs(1:nlambda1)-f(1:nlambda1,i))**2)
	endif

	lchi=log10(lchi)
        write(*,*)'i,lchi=',i,lchi
	write(15,*)lchi
enddo
close(15)
!write(*,*)'exiting chisurf'

end subroutine chisurf

