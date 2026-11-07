! Three-label interpolation (8575 pixels) and 12-block PCA (7514 pixels).
program benchmark
use share
use grid_geometry, only: sparse_active
implicit none
integer :: i,j,method,reps,iteration
real(dp) :: begin,finish,checksum,p(3)
real(dp),allocatable :: output(:)
character(20) :: arg
call get_command_argument(1,arg);read(arg,*) method
call get_command_argument(2,arg);read(arg,*) reps
checksum=0
if(method==5) then
    nelnpca=12;nvar=75;npix=900;totalnpca=7514;npca=627;npca(12)=617
    allocate(wpca(totalnpca,nvar),meanspca(totalnpca),output(totalnpca))
    meanspca=.8_dp
    do j=1,nvar
        do i=1,totalnpca
            wpca(i,j)=sin(real(i*j,dp)*.001_dp)*.01_dp
        enddo
    enddo
    allocate(f(npix,1));f(:,1)=.1_dp
    call cpu_time(begin)
    do iteration=1,reps
        f(1,1)=.1_dp+real(mod(iteration,101),dp)*.0001_dp
        call decompress(f(:,1),output)
        checksum=checksum+output(1+mod(iteration,totalnpca))
    enddo
else
    ndim=3;npix=8575;ntot=125;n_p(1:3)=5;inter=method;f_access=0
    allocate(ntimes(3),f(npix,ntot),output(npix))
    indi(1:3)=[3,2,1];ntimes=[25,5,1]
    do j=1,ntot
        do i=1,npix
            f(i,j)=1._dp+.05_dp*sin(real(i*j,dp)*.001_dp)
        enddo
    enddo
    call getee
    call cpu_time(begin)
    do iteration=1,reps
        p=[.1_dp+.79_dp*real(mod(iteration*7,101),dp)/101._dp, &
           .12_dp+.75_dp*real(mod(iteration*13,103),dp)/103._dp, &
           .08_dp+.8_dp*real(mod(iteration*17,107),dp)/107._dp]
        select case(method)
        case(1);call lin(p,output)
        case(2);call qua(p,output)
        case(3);call cub(p,output)
        case(4);call spl(p,output)
        end select
        checksum=checksum+output(1+mod(iteration,npix))
    enddo
endif
call cpu_time(finish)
write(*,'(i2,1x,i7,1x,f12.6,1x,es16.8)') method,reps,finish-begin,checksum
end program
