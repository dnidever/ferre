! Experimental cache-friendly PCA reconstruction.
! Same per-pixel summation order as src/decompress.f90; vectorizes across pixels.
subroutine decompress(influx,outflux)
use share, only: dp,scalef,scaled,npix,npca,nelnpca,nvar,totalnpca, &
                 constant,wpca,meanspca
implicit none
real(dp),intent(in) :: influx(npix)
real(dp),intent(out) :: outflux(totalnpca)
real(dp) :: flux(npix)
integer :: ii,jj,i,kk,offset
flux=influx
if (scaled == 1) flux=flux*scalef
flux=flux-constant
outflux=0.0_dp
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
outflux=outflux+meanspca
end subroutine decompress
