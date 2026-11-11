module sparse_interpolation
use, intrinsic :: ieee_arithmetic
use share, only: dp, ndim, npix, n_p, inter, f, f_access, f_format, fmtformat, npca, transposed, ntot
use grid_geometry
use quantized_grid, only: quantized_active, quantized_spectrum
use booklib, only: spline_fit, spline_int
implicit none
private
public :: interpolate_sparse

contains

subroutine interpolate_sparse(p, flux, valid)
real(dp), intent(in) :: p(ndim)
real(dp), intent(out) :: flux(npix)
logical, intent(out) :: valid
integer :: indices(4,ndim), counts(ndim), key(ndim), d
real(dp) :: weights(4,ndim), coordinate
flux = 0._dp
valid = all(ieee_is_finite(p))
if (.not.valid) return
do d=1,ndim
    coordinate = axis_values(1,d)+p(d)*(axis_values(n_p(d),d)-axis_values(1,d))
    if (p(d) < 0._dp .or. p(d) > 1._dp) then
        valid = .false.
        return
    endif
    call axis_weights(d, coordinate, indices(:,d), weights(:,d), counts(d))
enddo
call accumulate(1, 1._dp)
if (.not.all(ieee_is_finite(flux))) valid = .false.

contains

recursive subroutine accumulate(dimension, weight)
integer, intent(in) :: dimension
real(dp), intent(in) :: weight
real(dp) :: spectrum(npix)
real(dp), allocatable :: pixel_row(:)
integer :: i, record, ios, pixel
if (.not.valid) return
if (dimension <= ndim) then
    do i=1,counts(dimension)
        if (weights(i,dimension) == 0._dp) cycle
        key(dimension) = indices(i,dimension)
        call accumulate(dimension+1,weight*weights(i,dimension))
        if (.not.valid) return
    enddo
else
    record = node_record(key)
    if (record == 0) then
        valid = .false.
        return
    endif
    if (quantized_active) then
        call quantized_spectrum(record,spectrum)
    else if (f_access == 0) then
        spectrum = f(:,record)
    else
        !$omp critical(sparse_file_read)
        if (transposed/=0) then
            allocate(pixel_row(ntot))
            do pixel=1,npix
                if (f_format==0) then
                    read(10,fmtformat,rec=pixel,iostat=ios) pixel_row
                else
                    read(10,rec=pixel,iostat=ios) pixel_row
                endif
                if (ios/=0) exit
                spectrum(pixel)=pixel_row(record)
            enddo
            deallocate(pixel_row)
        else if (f_format == 0) then
            read(10,fmtformat,rec=record,iostat=ios) spectrum
        else
            read(10,rec=record,iostat=ios) spectrum
        endif
        !$omp end critical(sparse_file_read)
        if (ios /= 0) call grid_error('Cannot read sparse spectrum record')
    endif
    if (.not.all(ieee_is_finite(spectrum))) then
        valid = .false.
        return
    endif
    ! Zero PCA coefficients can describe a valid mean spectrum.
    if (npca(1) == 0 .and. all(spectrum == 0._dp)) then
        valid = .false.
        return
    endif
    flux = flux+weight*spectrum
endif
end subroutine
end subroutine

subroutine axis_weights(d, x, indices, weights, count)
integer, intent(in) :: d
real(dp), intent(in) :: x
integer, intent(out) :: indices(4), count
real(dp), intent(out) :: weights(4)
integer :: lower, start, i, selected, error
real(dp) :: xx(4), yy(4), second(4), t, value, tolerance
indices = 1
weights = 0._dp
lower = bracket(d,x)
! Snap physical coordinates within a few floating-point ulps of a node.
do i=lower,min(lower+1,n_p(d))
    tolerance = 4._dp*epsilon(x)*max(1._dp,abs(x),abs(axis_values(i,d)))
    if (abs(x-axis_values(i,d)) <= tolerance) then
        count = 1
        indices(1) = i
        weights(1) = 1._dp
        return
    endif
enddo
if (inter <= 0) then
    selected = lower
    if (x-axis_values(lower,d) >= axis_values(lower+1,d)-x) selected = lower+1
    count = 1
    indices(1) = selected
    weights(1) = 1._dp
    return
endif
if (inter == 1) then
    count = 2
    indices(1:2) = [lower,lower+1]
    t = (x-axis_values(lower,d))/(axis_values(lower+1,d)-axis_values(lower,d))
    weights(1:2) = [1._dp-t,t]
    return
endif
if (inter == 2) then
    count = 3
    start = max(1,lower-1)
else
    count = 4
    start = min(max(1,lower-1),n_p(d)-3)
endif
do i=1,count
    indices(i) = start+i-1
enddo
! Express distances relative to the enclosing interval for numerical stability.
xx(1:count) = (axis_values(indices(1:count),d)-axis_values(lower,d))/ &
              (axis_values(lower+1,d)-axis_values(lower,d))
t = (x-axis_values(lower,d))/(axis_values(lower+1,d)-axis_values(lower,d))
do i=1,count
    yy = 0._dp
    yy(i) = 1._dp
    if (inter == 4) then
        error=0
        call spline_fit(xx,yy,4,second,error)
        if (error /= 0) call grid_error('Spline construction failed')
        error=0
        call spline_int(xx,yy,4,second,t,value,error)
        if (error /= 0) call grid_error('Spline evaluation failed')
    else if (count == 3) then
        value = quadratic(xx(1:3),yy(1:3),t,lower-start+1)
    else if (lower == 1) then
        value = quadratic(xx(1:3),yy(1:3),t,1)
    else if (lower == n_p(d)-1) then
        value = quadratic(xx(2:4),yy(2:4),t,2)
    else
        value = yy(2)*(1._dp-t)**3+yy(3)*t**3 + &
            3._dp*t*(1._dp-t)**2*(yy(2)+(yy(3)-yy(1))/(xx(3)-xx(1))/3._dp) + &
            3._dp*t**2*(1._dp-t)*(yy(3)-(yy(4)-yy(2))/(xx(4)-xx(2))/3._dp)
    endif
    weights(i) = value
enddo
end subroutine

real(dp) function quadratic(x, y, t, left) result(value)
real(dp), intent(in) :: x(3), y(3), t
integer, intent(in) :: left
real(dp) :: derivative, control
derivative = (y(3)-y(1))/(x(3)-x(1))
if (left == 1) then
    control = y(2)-.5_dp*derivative
    value = y(1)*(1._dp-t)**2+y(2)*t**2+2._dp*t*(1._dp-t)*control
else
    control = y(2)+.5_dp*derivative
    value = y(2)*(1._dp-t)**2+y(3)*t**2+2._dp*t*(1._dp-t)*control
endif
end function

end module sparse_interpolation
