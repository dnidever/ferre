subroutine decompress(influx,outflux)

use share, only: dp,scalef,scaled,npix,                                      &
   npca,nelnpca,nvar,totalnpca,f_access,constant,pcaproject,wpca,meanspca

implicit none

!in/out
real(dp),intent(in)	   :: influx(npix)		!synthetic fluxes
real(dp),intent(out)   :: outflux(totalnpca)		!npca expanded flux vector

!locals
integer				   :: i,ii,jj,kk,offset
real(dp)               :: flux(npix)			!tmp data holder 


flux=influx

if (scaled == 1) flux=flux*scalef
flux = flux - constant

! Reconstruct each block while walking contiguous wavelength elements.
! Each pixel retains the original component accumulation order.
outflux=0._dp
offset=0
do ii=1,nelnpca
    do jj=1,nvar
        do i=1,npca(ii)
            kk=i+offset
            outflux(kk)=outflux(kk)+flux((ii-1)*nvar+jj)*wpca(kk,jj)
        enddo
    enddo
    offset=offset+npca(ii)
enddo
outflux = outflux + meanspca
	

end subroutine decompress
