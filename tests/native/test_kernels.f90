! Numerical regression tests for PCA traversal and reusable spline weights.
! Compare against per-pixel reference calculations, not timing thresholds.
program test_kernels
use share
use omp_lib, only: omp_set_num_threads
implicit none
integer :: components, index, trial, d, order, mode, point
integer :: pixel_count, state_size
integer, allocatable :: seed(:)
real(dp), allocatable :: coefficients(:), output(:), expected(:), points(:,:), reference(:,:), actual(:,:)
real(dp) :: error, worst
call random_seed(size=state_size)
allocate(seed(state_size));seed=271828;call random_seed(put=seed)
call omp_set_num_threads(3)
worst=0._dp

! PCA: scalar accumulation is an independent reference, including scaling,
! means, constant subtraction, unequal block sizes, and multiple components.
nelnpca=3
npca(1:3)=[17,31,9]
totalnpca=sum(npca(1:3))
allocate(meanspca(totalnpca))
call random_number(meanspca)
do components=1,75,37
    nvar=components;npix=nelnpca*nvar
    allocate(wpca(totalnpca,nvar),coefficients(npix),output(totalnpca),expected(totalnpca))
    call random_number(wpca);wpca=(wpca-.5_dp)*.2_dp
    do mode=0,1
        scaled=mode;scalef=1.37_dp;constant=.17_dp
        do trial=1,16
            call random_number(coefficients);coefficients=2._dp*coefficients-1._dp
            if (trial==1) coefficients=0._dp
            call reference_pca(coefficients,expected)
            call decompress(coefficients,output)
            call check(output,expected,'PCA reconstruction')
        enddo
    enddo
    ! Shared PCA basis with distinct inputs on simultaneous fitting threads.
    allocate(points(npix,24),reference(totalnpca,24),actual(totalnpca,24))
    call random_number(points)
    do point=1,24
        call reference_pca(points(:,point),reference(:,point))
    enddo
    !$omp parallel do private(point)
    do point=1,24
        call decompress(points(:,point),actual(:,point))
    enddo
    !$omp end parallel do
    call check(reshape(actual,[size(actual)]),reshape(reference,[size(reference)]),'parallel PCA')
    deallocate(points,reference,actual,wpca,coefficients,output,expected)
enddo
deallocate(meanspca)

! Splines: compare the old per-pixel spline solve at randomized locations,
! exact nodes, both boundaries, and differing interpolation-axis orders.
! Exercise in-memory and binary direct-access libraries with three threads.
npix=137
allocate(output(npix),expected(npix))
do d=1,5
    ndim=d;n_p(1:ndim)=5;ntot=5**ndim
    allocate(ntimes(ndim),f(npix,ntot),points(ndim,32),reference(npix,32),actual(npix,32))
    do index=1,ndim
        ntimes(index)=5**(ndim-index)
    enddo
    call random_number(f);f=.5_dp+f
    call random_number(points)
    points(:,1)=0._dp
    points(:,2)=1._dp-2._dp*epsilon(1._dp)
    points(:,3)=.25_dp
    points(:,4)=.5_dp
    points(:,5)=.75_dp
    do order=1,2
        do index=1,ndim
            if (order==1) then
                indi(index)=ndim-index+1
            else
                indi(index)=index
            endif
        enddo
        inter=4
        call getee
        do mode=0,1
            f_access=mode;f_format=1
            if (mode==1) then
                inquire(iolength=pixel_count) f(:,1)
                open(10,status='scratch',access='direct',form='unformatted',recl=pixel_count)
                do index=1,ntot
                    write(10,rec=index) f(:,index)
                enddo
            endif
            do point=1,32
                call reference_spl(points(:,point),reference(:,point))
            enddo
            !$omp parallel do private(point)
            do point=1,32
                call spl(points(:,point),actual(:,point))
            enddo
            !$omp end parallel do
            call check(reshape(actual,[size(actual)]),reshape(reference,[size(reference)]),'spline interpolation')
            if (mode==1) close(10)
        enddo
    enddo
    deallocate(ntimes,f,points,reference,actual)
enddo
write(*,'(a,es12.4)') 'PCA and spline regression checks passed; largest absolute difference: ',worst

contains
subroutine check(values,truth,label)
real(dp), intent(in) :: values(:),truth(:)
character(*), intent(in) :: label
! Both paths use single precision by default. Permit a few ulps through
! multiple tensor-product stages without accepting NaN or infinite outputs.
use_ieee: block
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    if (.not.all(ieee_is_finite(values))) then
        write(*,*) 'Nonfinite output: ',label
        stop 1
    endif
end block use_ieee
error=maxval(abs(values-truth))
worst=max(worst,error)
if (any(abs(values-truth)>32._dp*epsilon(1._dp)*max(1._dp,abs(truth)))) then
    write(*,*) 'Numerical mismatch: ',label,error
    stop 1
endif
end subroutine

subroutine reference_pca(input,result)
real(dp), intent(in) :: input(npix)
real(dp), intent(out) :: result(totalnpca)
real(dp) :: work(npix)
integer :: b,i,j,offset,k
work=input
if (scaled==1) work=work*scalef
work=work-constant
offset=0
do b=1,nelnpca
    do i=1,npca(b)
        k=offset+i;result(k)=0._dp
        do j=1,nvar
            result(k)=result(k)+work((b-1)*nvar+j)*wpca(k,j)
        enddo
    enddo
    offset=offset+npca(b)
enddo
result=result+meanspca
end subroutine
end program
