module grid_geometry
! Geometry is constructed serially and is read-only during parallel fits.
use iso_fortran_env, only: int64
use, intrinsic :: ieee_arithmetic
use share, only: dp, ndim, n_p, ntot, llimits, steps, badflux, flen
implicit none
private
public :: sparse_active, axis_values, node_coordinates, load_nodes, clear_nodes, &
          node_record, bracket, grid_error, invalid_model, fit_penalty, &
          read_spectrum, validate_binary_size, coordinate_path, grid_upper
logical :: sparse_active = .false.
real(dp), allocatable :: axis_values(:,:), node_coordinates(:,:)
integer, allocatable :: node_keys(:,:), slots(:)
real(dp), parameter :: fit_penalty = 1.e20_dp

contains

subroutine grid_error(message)
character(*), intent(in) :: message
write(*,'(a)') 'FERRE grid error: '//trim(message)
error stop 1
end subroutine

subroutine clear_nodes()
if (allocated(axis_values)) deallocate(axis_values, node_coordinates, node_keys, slots)
sparse_active = .false.
end subroutine

function coordinate_path(library, filename) result(path)
character(*), intent(in) :: library, filename
character(:), allocatable :: path
integer :: slash
path = trim(filename)
if (filename(1:1) == '/') return
slash = index(trim(library), '/', back=.true.)
if (slash > 0) path = library(:slash)//trim(filename)
end function

recursive subroutine sort_values(values, first, last)
real(dp), intent(inout) :: values(:)
integer, intent(in) :: first, last
integer :: left, right
real(dp) :: pivot, swap
if (first >= last) return
left = first
right = last
pivot = values(first + (last-first)/2)
do
    do while (values(left) < pivot)
        left = left+1
    enddo
    do while (values(right) > pivot)
        right = right-1
    enddo
    if (left > right) exit
    swap = values(left)
    values(left) = values(right)
    values(right) = swap
    left = left+1
    right = right-1
    if (left > right) exit
enddo
if (first < right) call sort_values(values, first, right)
if (left < last) call sort_values(values, left, last)
end subroutine

subroutine load_nodes(filename, library)
character(*), intent(in) :: filename, library
character(:), allocatable :: path
character(len=32768) :: line
real(dp), allocatable :: sorted(:,:), row(:)
integer :: unit, ios, records, i, j, k, count, position, words
logical :: in_word
call clear_nodes()
if (len_trim(filename) == 0) return
path = coordinate_path(library, filename)
open(newunit=unit, file=path, status='old', action='read', iostat=ios)
if (ios /= 0) call grid_error('Cannot open coordinate file '//path)
records = 0
! Count and validate data rows. Comments and blank lines do not consume records.
do
    call coordinate_line(unit, line, ios)
    if (ios < 0) exit
    if (ios > 0) call grid_error('Error reading '//path)
    if (len_trim(line) == 0) cycle
    words = 0
    in_word = .false.
    do i=1,len_trim(line)
        if (line(i:i) == ' ' .or. line(i:i) == achar(9) .or. line(i:i) == ',') then
            in_word = .false.
        else if (.not.in_word) then
            words = words+1
            in_word = .true.
        endif
    enddo
    if (words /= ndim) call grid_error('Coordinate row must have exactly NDIM values in '//path)
    records = records+1
enddo
if (records == 0) call grid_error('Empty coordinate file '//path)
allocate(node_coordinates(ndim,records), sorted(records,ndim), row(ndim))
rewind(unit)
j = 0
do
    call coordinate_line(unit, line, ios)
    if (ios < 0) exit
    if (len_trim(line) == 0) cycle
    read(line,*,iostat=ios) row
    if (ios /= 0) call grid_error('Invalid coordinate row in '//path)
    if (.not.all(ieee_is_finite(row))) call grid_error('Nonfinite coordinate in '//path)
    j = j+1
    node_coordinates(:,j) = row
enddo
close(unit)
do i=1,ndim
    sorted(:,i) = node_coordinates(i,:)
    call sort_values(sorted(:,i), 1, records)
    count = 1
    do j=2,records
        if (sorted(j,i) /= sorted(count,i)) then
            count = count+1
            sorted(count,i) = sorted(j,i)
        endif
    enddo
    n_p(i) = count
    llimits(i) = sorted(1,i)
    steps(i) = 0._dp
    if (count > 1) steps(i) = (sorted(count,i)-sorted(1,i))/real(count-1,dp)
enddo
allocate(axis_values(maxval(n_p(1:ndim)),ndim), node_keys(ndim,records), slots(2*records+1))
axis_values = 0._dp
slots = 0
do i=1,ndim
    axis_values(1:n_p(i),i) = sorted(1:n_p(i),i)
    do j=1,records
        k = bracket(i,node_coordinates(i,j))
        if (axis_values(k,i) /= node_coordinates(i,j)) k = k+1
        node_keys(i,j) = k
    enddo
enddo
do j=1,records
    position = hash_key(node_keys(:,j))
    do while (slots(position) /= 0)
        if (all(node_keys(:,slots(position)) == node_keys(:,j))) &
            call grid_error('Duplicate coordinate tuple in '//path)
        position = mod(position,size(slots))+1
    enddo
    slots(position) = j
enddo
ntot = records
sparse_active = .true.
write(*,'(a,i0,a)') 'Loaded sparse coordinate table with ',records,' spectra'
end subroutine

subroutine coordinate_line(unit, line, ios)
integer, intent(in) :: unit
character(*), intent(out) :: line
integer, intent(out) :: ios
integer :: comment
read(unit,'(a)',iostat=ios) line
if (ios /= 0) return
comment = scan(line,'#!')
if (comment > 0) line(comment:) = ' '
end subroutine

integer function bracket(dimension, value) result(lower)
integer, intent(in) :: dimension
real(dp), intent(in) :: value
integer :: left, right, middle
left = 1
right = n_p(dimension)
do while (left < right-1)
    middle = left+(right-left)/2
    if (axis_values(middle,dimension) <= value) then
        left = middle
    else
        right = middle
    endif
enddo
lower = left
end function

integer function hash_key(key) result(position)
integer, intent(in) :: key(:)
integer(int64) :: hash, capacity
integer :: i
hash = 0_int64
capacity = int(size(slots),int64)
do i=1,size(key)
    hash = modulo(hash*1000003_int64 + int(key(i),int64),capacity)
enddo
position = int(hash)+1
end function

integer function node_record(key) result(record)
integer, intent(in) :: key(:)
integer :: position, start
record = 0
position = hash_key(key)
start = position
do while (slots(position) /= 0)
    if (all(node_keys(:,slots(position)) == key)) then
        record = slots(position)
        return
    endif
    position = mod(position,size(slots))+1
    if (position == start) exit
enddo
end function

real(dp) function grid_upper(d) result(upper)
integer, intent(in) :: d
if (sparse_active) then
    upper=axis_values(n_p(d),d)
else
    upper=llimits(d)+steps(d)*(n_p(d)-1)
endif
end function

logical function invalid_model(flux)
real(dp), intent(in) :: flux(:)
invalid_model = all(flux == badflux) .or. .not.all(ieee_is_finite(flux))
end function

subroutine read_spectrum(unit, values)
! Retain wrapped list-directed spectra; accept a one-token zero marker.
integer, intent(in) :: unit
real(dp), intent(out) :: values(:)
character(len=32768) :: line
real(dp) :: zero
integer :: ios, words, i
logical :: in_word
read(unit,'(a)',iostat=ios) line
if (ios /= 0) call grid_error('Too few spectrum records')
words = 0
in_word = .false.
do i=1,len_trim(line)
    if (line(i:i) == ' ' .or. line(i:i) == achar(9)) then
        in_word = .false.
    else if (.not.in_word) then
        words = words+1
        in_word = .true.
    endif
enddo
if (words == 1) then
    read(line,*,iostat=ios) zero
    if (ios == 0) then
        if (zero == 0._dp) then
            values = 0._dp
            return
        endif
    endif
endif
backspace(unit)
read(unit,*,iostat=ios) values
if (ios /= 0) call grid_error('Incomplete or invalid spectrum record')
end subroutine

subroutine validate_binary_size(filename, pixels, records)
character(*), intent(in) :: filename
integer, intent(in) :: pixels, records
integer(int64) :: bytes, expected
integer :: ios
inquire(file=filename,size=bytes,iostat=ios)
expected = int(pixels,int64)*int(records,int64)*int(storage_size(1._dp)/8,int64)
if (ios /= 0 .or. bytes /= expected) call grid_error('Binary record count does not match NODEFILE: '//filename)
end subroutine

end module grid_geometry
