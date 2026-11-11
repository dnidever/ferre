module quantized_grid
! Native PCA coefficient storage. The persistent grid is int16, never real(dp).
! On disk values are unsigned, little-endian 16-bit words. Signed Fortran
! storage preserves their bit patterns; masking into int32 recovers 0..65535.
use iso_fortran_env, only: int8, int16, int32, int64, real64, file_storage_size
use, intrinsic :: ieee_arithmetic
use share, only: dp, npix, ntot, f_access, transposed, f
use grid_geometry, only: grid_error, coordinate_path
implicit none
private
public :: quantized_active, open_quantized, quantized_spectrum, reset_quantized, quantized_bytes
logical :: quantized_active = .false., swap_bytes = .false.
integer :: data_unit = 0
integer(int16), allocatable :: codes(:,:)
real(real64), allocatable :: offset(:), scale(:)

contains

subroutine reset_quantized()
if (data_unit /= 0) close(data_unit)
data_unit = 0
if (allocated(codes)) deallocate(codes)
if (allocated(offset)) deallocate(offset, scale)
quantized_active = .false.
end subroutine

function quantized_bytes() result(bytes)
integer(int64) :: bytes
bytes = 0_int64
if (allocated(codes)) bytes = size(codes,kind=int64)*2_int64
end function

subroutine endian_swap(words)
integer(int16), intent(inout) :: words(:)
integer :: i
integer(int32) :: word
if (.not.swap_bytes) return
do i=1,size(words)
    word = iand(int(words(i),int32),65535_int32)
    word = ior(ishft(iand(word,255_int32),8),ishft(word,-8))
    if (word >= 32768_int32) word = word-65536_int32
    words(i) = int(word,int16)
enddo
end subroutine

subroutine open_quantized(header, table)
character(*), intent(in) :: header, table
character(:), allocatable :: path, data_path
integer :: version, ncoeff, unit, ios, i, stat
integer(int64) :: bytes, expected
integer(int16) :: one
integer(int8) :: representation(2)
integer(int16), allocatable :: row(:)
namelist / quantization / version, ncoeff, offset, scale
call reset_quantized()
if (allocated(f)) call grid_error('uint16 reader must not retain a float32 whole-grid array')
if (storage_size(1._dp) /= 32 .or. file_storage_size /= 8 .or. storage_size(one) /= 16) &
    call grid_error('uint16 grids require 32-bit FERRE real values and 8-bit file storage units')
if (f_access /= 0 .and. f_access /= 1) call grid_error('uint16 requires F_ACCESS=0 or 1')
if (transposed /= 0 .and. transposed /= 1) call grid_error('uint16 requires TRANSPOSED=0 or 1')
if (len_trim(table) == 0) call grid_error('uint16 requires QUANTFILE')
allocate(offset(npix),scale(npix),stat=stat)
if (stat /= 0) call grid_error('Cannot allocate quantization parameter table')
offset = ieee_value(0._real64,ieee_quiet_nan)
scale = offset
version = 0
ncoeff = 0
path = coordinate_path(header,table)
open(newunit=unit,file=path,status='old',action='read',iostat=ios)
if (ios /= 0) call grid_error('Cannot open QUANTFILE: '//path)
read(unit,nml=quantization,iostat=ios)
close(unit)
if (ios /= 0) call grid_error('Cannot read QUANTIZATION table: '//path)
if (version /= 1 .or. ncoeff /= npix) call grid_error('Quantization version or coefficient count mismatch')
if (.not.all(ieee_is_finite(offset)) .or. .not.all(ieee_is_finite(scale))) &
    call grid_error('Quantization offsets/scales must all be finite and present')
if (any(scale <= 0._real64)) call grid_error('Quantization scales must be positive')
if (any(abs(offset) > real(huge(1._dp),real64)) .or. &
    any(abs(offset+scale*65535._real64) > real(huge(1._dp),real64))) &
    call grid_error('Quantization endpoints exceed the working real range')
one = 1_int16
representation = transfer(one,representation)
swap_bytes = representation(1) /= 1_int8
i = index(trim(header),'.',back=.true.)
if (i == 0) call grid_error('uint16 header must have a file extension')
data_path = header(:i)//'u16'
expected = 2_int64*int(npix,int64)*int(ntot,int64)
inquire(file=data_path,size=bytes,iostat=ios)
if (ios /= 0 .or. bytes /= expected) call grid_error('uint16 payload size mismatch: '//data_path)
open(newunit=data_unit,file=data_path,status='old',access='stream',form='unformatted',action='read',iostat=ios)
if (ios /= 0) call grid_error('Cannot open uint16 payload: '//data_path)
if (f_access == 0) then
    allocate(codes(npix,ntot),stat=stat)
    if (stat /= 0) call grid_error('Cannot allocate uint16 coefficient grid')
    if (transposed == 0) then
        read(data_unit,iostat=ios) codes
        if (ios /= 0) call grid_error('Cannot load uint16 coefficient grid')
    else
        allocate(row(ntot),stat=stat)
        if (stat /= 0) call grid_error('Cannot allocate transposed uint16 row')
        do i=1,npix
            read(data_unit,iostat=ios) row
            if (ios /= 0) call grid_error('Cannot load transposed uint16 coefficient row')
            codes(i,:) = row
        enddo
        deallocate(row)
    endif
    if (swap_bytes) then
        do i=1,ntot
            call endian_swap(codes(:,i))
        enddo
    endif
    close(data_unit)
    data_unit = 0
endif
quantized_active = .true.
write(*,'(a,i0,a)') 'read_f       -> uint16 coefficient array: ',quantized_bytes(),' RAM bytes'
write(*,'(a)') 'read_f       -> only requested neighboring coefficient records are decoded'
end subroutine

subroutine quantized_spectrum(record, spectrum)
integer, intent(in) :: record
real(dp), intent(out) :: spectrum(npix)
integer(int16) :: buffer(npix)
integer :: ios, pixel
integer(int64) :: position
if (record < 1 .or. record > ntot) call grid_error('uint16 record outside grid')
if (allocated(codes)) then
    buffer = codes(:,record)
else
    ios = 0
    !$omp critical(quantized_file_read)
    if (transposed == 0) then
        position = 1_int64+2_int64*int(record-1,int64)*int(npix,int64)
        read(data_unit,pos=position,iostat=ios) buffer
    else
        do pixel=1,npix
            position = 1_int64+2_int64*(int(pixel-1,int64)*int(ntot,int64)+int(record-1,int64))
            read(data_unit,pos=position,iostat=ios) buffer(pixel)
            if (ios /= 0) exit
        enddo
    endif
    !$omp end critical(quantized_file_read)
    if (ios /= 0) call grid_error('Cannot read uint16 coefficient record')
    call endian_swap(buffer)
endif
! Offsets/scales are float64, then the small restored vector is rounded to dp.
! There is deliberately no decoded whole-grid allocation or cache.
spectrum = real(offset+scale*real(iand(int(buffer,int32),65535_int32),real64),dp)
end subroutine

end module quantized_grid
