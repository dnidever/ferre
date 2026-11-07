MODULE booklib
IMPLICIT NONE

! Set access--private by default, and public for specific
! procedures.
PRIVATE
PUBLIC cross_prod                 ! Cross product
PUBLIC deriv                      ! Derivative
PUBLIC dft                        ! DFT: single / double
PUBLIC fft                        ! FFT: single / double
PUBLIC heapsort                   ! Heap sort
PUBLIC heapsort_2                 ! Heap sort w/carry
PUBLIC histogram                  ! Histogram
PUBLIC idft                       ! Inv. DFT: single / double
PUBLIC ifft                       ! Inv. FFT: single / double
PUBLIC integ                      ! Integration of a function
PUBLIC integ_d                    ! Integration of discrete points
PUBLIC interp                     ! Interpolation
PUBLIC lcase                      ! Shift string to lower case
PUBLIC lsq_fit                    ! Least squares fit to poly
PUBLIC lstmul                     ! Last power of n
PUBLIC mat_inv                    ! Matrix inversion
PUBLIC nxtmul                     ! Next power of n
PUBLIC plot                       ! Plot a data set
PUBLIC plotxy                     ! Plot (x,y) data
PUBLIC random_n                   ! Gaussian random variable
PUBLIC random_r                   ! Rayleigh random variable
PUBLIC random_u                   ! Uniform random variable
PUBLIC simul                      ! Solve simultaneous eqns
PUBLIC sinc                       ! SINC function
PUBLIC spline_fit                 ! Cubic spline fit
PUBLIC spline_int                 ! Cubic spline interpolation
PUBLIC ssort                      ! Selection sort
PUBLIC statistics                 ! Calculate statistics
PUBLIC ucase                      ! Shift string to upper case

! Declare generic procedures.
INTERFACE cross_prod
   MODULE PROCEDURE cross_prod_sgl
   MODULE PROCEDURE cross_prod_dbl
END INTERFACE

INTERFACE deriv
   MODULE PROCEDURE deriv
   MODULE PROCEDURE dderiv
END INTERFACE

INTERFACE dft
   MODULE PROCEDURE dft_sgl
   MODULE PROCEDURE dft_dbl
END INTERFACE

INTERFACE fft
   MODULE PROCEDURE fft_sgl
   MODULE PROCEDURE fft_dbl
END INTERFACE

INTERFACE heapsort
   MODULE PROCEDURE heapsort_real_sgl
   MODULE PROCEDURE heapsort_real_dbl
   MODULE PROCEDURE heapsort_int
   MODULE PROCEDURE heapsort_char
END INTERFACE

INTERFACE heapsort_2
   MODULE PROCEDURE heapsort_2_real_sgl
   MODULE PROCEDURE heapsort_2_real_dbl
   MODULE PROCEDURE heapsort_2_int
   MODULE PROCEDURE heapsort_2_char
END INTERFACE

INTERFACE idft
   MODULE PROCEDURE idft_sgl
   MODULE PROCEDURE idft_dbl
END INTERFACE

INTERFACE ifft
   MODULE PROCEDURE ifft_sgl
   MODULE PROCEDURE ifft_dbl
END INTERFACE

INTERFACE integ
   MODULE PROCEDURE integrate_sgl
   MODULE PROCEDURE integrate_dbl
END INTERFACE

INTERFACE integ_d
   MODULE PROCEDURE integrate_d_sgl
   MODULE PROCEDURE integrate_d_dbl
END INTERFACE

INTERFACE interp
   MODULE PROCEDURE interp_sgl
   MODULE PROCEDURE interp_dbl
END INTERFACE

INTERFACE lsq_fit
   MODULE PROCEDURE lsq_fit_sgl
   MODULE PROCEDURE lsq_fit_dbl
END INTERFACE

INTERFACE mat_inv
   MODULE PROCEDURE mat_inv_sgl
   MODULE PROCEDURE mat_inv_dbl
   MODULE PROCEDURE mat_inv_sgl_cmplx
   MODULE PROCEDURE mat_inv_dbl_cmplx
END INTERFACE

INTERFACE plot
   MODULE PROCEDURE plot_sgl
   MODULE PROCEDURE plot_dbl
END INTERFACE

INTERFACE plotxy
   MODULE PROCEDURE plotxy_sgl
   MODULE PROCEDURE plotxy_dbl
END INTERFACE

INTERFACE simul
   MODULE PROCEDURE simul_sgl
   MODULE PROCEDURE simul_dbl
   MODULE PROCEDURE simul_sgl_cmplx
   MODULE PROCEDURE simul_dbl_cmplx
END INTERFACE

INTERFACE sinc
   MODULE PROCEDURE sinc_sgl
   MODULE PROCEDURE sinc_dbl
END INTERFACE

INTERFACE spline_fit
    MODULE PROCEDURE spline_fit_sgl
    MODULE PROCEDURE spline_fit_dbl
END INTERFACE

INTERFACE spline_int
    MODULE PROCEDURE spline_int_sgl
    MODULE PROCEDURE spline_int_dbl
END INTERFACE

INTERFACE ssort
   MODULE PROCEDURE ssort_sgl
   MODULE PROCEDURE ssort_dbl
   MODULE PROCEDURE ssort_int
   MODULE PROCEDURE ssort_char
END INTERFACE

INTERFACE statistics
   MODULE PROCEDURE statistics_sgl
   MODULE PROCEDURE statistics_dbl
END INTERFACE

! Declare module procedures.
CONTAINS

SUBROUTINE fft_sgl ( array_in, array_out, n, error )
!
!  Purpose:
!    A Cooley-Tukey Radix 2, in place, decimation-in-frequency FFT
!    subroutine.  Although the routine is in place, the input
!    array is not destroyed because it is first copied to the
!    output array, and the FFT is performed on the copy in the
!    output array.
!
!  Record of revisions:
!       Date       Programmer          Description of change
!       ====       ==========          =====================
!     03/18/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=6) ! Precision
REAL(KIND=kind), PARAMETER :: twopi = 6.283185307179586_kind ! 2*pi

! Declare calling arguments:
INTEGER, INTENT(IN) :: n             ! No. of elements in FFT--must
                                     ! be a power of 2:  2, 4, 8, etc.
COMPLEX(KIND=kind),INTENT(IN),DIMENSION(n) :: array_in
                                     ! Input data for FFT
COMPLEX(KIND=kind),INTENT(INOUT),DIMENSION(n) :: array_out
                                     ! Output data for FFT
INTEGER, INTENT(OUT) :: error        ! Error flag:
                                     ! 0 -- No errors
                                     ! 1 -- n not a power of 2.

! List of local variables:
REAL(KIND=kind) :: angle      ! Twiddle factor angle for butterfly
REAL(KIND=kind) :: anginc     ! Angle increment for twiddle factor
INTEGER :: i, j, k, l, n1, n2 ! Index variables
INTEGER :: m                  ! Number such that n = 2**m
COMPLEX(KIND=kind) :: temp    ! Temporary storage
COMPLEX(KIND=kind) :: w       ! Twiddle factor

! Copy the input array
array_out = array_in

! Check for proper array size.  The length of the input data
! array must be a power of 2.
j = 1
m = 0
1 IF ( j < n ) THEN
     M = M + 1
     J = J * 2
     GO TO 1
  END IF

! Make sure N is a power of 2.
IF ( J /= N ) THEN

   ! Not a power of 2--this is a error!
   error = 1

ELSE

   ! This is a power of 2.  Process the FFT!
   N2 = N

   ! Loop over the number of stages in the FFT.
   DO 10 K = 1, M
       N1     = N2
       N2     = N2 / 2
       ANGINC = TWOPI / N1
       ANGLE  = 0.

       ! Loop over all angles at the current stage.
       DO 20 J = 1, N2
          W = CMPLX ( COS(ANGLE), -SIN(ANGLE), KIND=kind )
          ANGLE = J * ANGINC

          ! Perform butterflies.
          DO 30 I = J, N, N1
             L = I + N2
             TEMP         = ARRAY_out(I) - ARRAY_out(L)
             ARRAY_out(I) = ARRAY_out(I) + ARRAY_out(L)
             ARRAY_out(L) = TEMP * W
          30 CONTINUE
       20 CONTINUE
    10 CONTINUE

   ! Unscramble data with a digit reverse counter.
   J = 1
   N1 = N - 1
   DO 104 I = 1, N1
      IF ( I .GE. J ) GO TO 101
         TEMP   = ARRAY_out(J)
         ARRAY_out(J) = ARRAY_out(I)
         ARRAY_out(I) = TEMP
         101 K = N / 2
         102 IF ( K .GE. J ) GO TO 103
              J = J - K
              K = K / 2
              GO TO 102
         103 J = J + K
   104 CONTINUE
END IF

END SUBROUTINE fft_sgl

SUBROUTINE fft_dbl ( array_in, array_out, n, error )
!
!  Purpose:
!    A Cooley-Tukey Radix 2, in place, decimation-in-frequency FFT
!    subroutine.  Although the routine is in place, the input
!    array is not destroyed because it is first copied to the
!    output array, and the FFT is performed on the copy in the
!    output array.
!
!  Record of revisions:
!       Date       Programmer          Description of change
!       ====       ==========          =====================
!     03/18/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=12) ! Precision
REAL(KIND=kind), PARAMETER :: twopi = 6.283185307179586_kind ! 2*pi

! Declare calling arguments:
INTEGER, INTENT(IN) :: n             ! No. of elements in FFT--must
                                     ! be a power of 2:  2, 4, 8, etc.
COMPLEX(KIND=kind),INTENT(IN),DIMENSION(n) :: array_in
                                     ! Input data for FFT
COMPLEX(KIND=kind),INTENT(INOUT),DIMENSION(n) :: array_out
                                     ! Output data for FFT
INTEGER, INTENT(OUT) :: error        ! Error flag:
                                     ! 0 -- No errors
                                     ! 1 -- n not a power of 2.

! List of local variables:
REAL(KIND=kind) :: angle      ! Twiddle factor angle for butterfly
REAL(KIND=kind) :: anginc     ! Angle increment for twiddle factor
INTEGER :: i, j, k, l, n1, n2 ! Index variables
INTEGER :: m                  ! Number such that n = 2**m
COMPLEX(KIND=kind) :: temp    ! Temporary storage
COMPLEX(KIND=kind) :: w       ! Twiddle factor

! Copy the input array
array_out = array_in

! Check for proper array size.  The length of the input data
! array must be a power of 2.
j = 1
m = 0
1 IF ( j < n ) THEN
     M = M + 1
     J = J * 2
     GO TO 1
  END IF

! Make sure N is a power of 2.
IF ( J /= N ) THEN

   ! Not a power of 2--this is a error!
   error = 1

ELSE

   ! This is a power of 2.  Process the FFT!
   N2 = N

   ! Loop over the number of stages in the FFT.
   DO 10 K = 1, M
       N1     = N2
       N2     = N2 / 2
       ANGINC = TWOPI / N1
       ANGLE  = 0.

       ! Loop over all angles at the current stage.
       DO 20 J = 1, N2
          W = CMPLX ( COS(ANGLE), -SIN(ANGLE), KIND=kind )
          ANGLE = J * ANGINC

          ! Perform butterflies.
          DO 30 I = J, N, N1
             L = I + N2
             TEMP         = ARRAY_out(I) - ARRAY_out(L)
             ARRAY_out(I) = ARRAY_out(I) + ARRAY_out(L)
             ARRAY_out(L) = TEMP * W
          30 CONTINUE
       20 CONTINUE
    10 CONTINUE

   ! Unscramble data with a digit reverse counter.
   J = 1
   N1 = N - 1
   DO 104 I = 1, N1
      IF ( I .GE. J ) GO TO 101
         TEMP   = ARRAY_out(J)
         ARRAY_out(J) = ARRAY_out(I)
         ARRAY_out(I) = TEMP
         101 K = N / 2
         102 IF ( K .GE. J ) GO TO 103
              J = J - K
              K = K / 2
              GO TO 102
         103 J = J + K
   104 CONTINUE
END IF

END SUBROUTINE fft_dbl

SUBROUTINE ifft_sgl ( array_in, array_out, n, error )
!
!  Purpose:
!    A Cooley-Tukey Radix 2, in place, decimation-in-frequency
!    inverse FFT subroutine.  Although the routine is in place,
!    the input array is not destroyed because it is first copied
!    to the output array, and the FFT is performed on the copy
!    in the output array.
!
!  Record of revisions:
!       Date       Programmer          Description of change
!       ====       ==========          =====================
!     05/17/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=6) ! Precision
REAL(KIND=kind), PARAMETER :: twopi = 6.283185307179586_kind ! 2*pi

! Declare calling arguments:
INTEGER, INTENT(IN) :: n             ! No. of elements in IFFT--must
                                     ! be a pwer of 2:  2, 4, 8, etc.
COMPLEX(KIND=kind),INTENT(IN),DIMENSION(n) :: array_in
                                     ! Input data for IFFT
COMPLEX(KIND=kind),INTENT(INOUT),DIMENSION(n) :: array_out
                                     ! Output data for IFFT
INTEGER, INTENT(OUT) :: error        ! Error flag:
                                     ! 0 -- No errors
                                     ! 1 -- n not a power of 2.

! List of local variables:
REAL(KIND=kind) :: angle      ! Twiddle factor angle for butterfly
REAL(KIND=kind) :: anginc     ! Angle increment for twiddle factor
INTEGER :: i, j, k, l, n1, n2 ! Index variables
INTEGER :: m                  ! Number such that n = 2**m
COMPLEX(KIND=kind) :: temp    ! Temporary storage
COMPLEX(KIND=kind) :: w       ! Twiddle factor

! Copy the input array
array_out = array_in

! Check for proper array size.  The length of the input data
! array must be a power of 2.
j = 1
m = 0
1 IF ( j < n ) THEN
     M = M + 1
     J = J * 2
     GO TO 1
  END IF

! Make sure N is a power of 2.
IF ( J /= N ) THEN

   ! Not a power of 2--this is a error!
   error = 1

ELSE

   ! This is a power of 2.  Process the FFT!
   N2 = N

   ! Loop over the number of stages in the FFT.
   DO 10 K = 1, M
       N1     = N2
       N2     = N2 / 2
       ANGINC = TWOPI / N1
       ANGLE  = 0.

       ! Loop over all angles at the current stage.
       DO 20 J = 1, N2
          W = CMPLX ( COS(ANGLE), SIN(ANGLE), KIND=kind )
          ANGLE = J * ANGINC

          ! Perform butterflies.
          DO 30 I = J, N, N1
             L = I + N2
             TEMP         = ARRAY_out(I) - ARRAY_out(L)
             ARRAY_out(I) = ARRAY_out(I) + ARRAY_out(L)
             ARRAY_out(L) = TEMP * W
          30 CONTINUE
       20 CONTINUE
    10 CONTINUE

   ! Unscramble data with a digit reverse counter.
   J = 1
   N1 = N - 1
   DO 104 I = 1, N1
      IF ( I .GE. J ) GO TO 101
         TEMP   = ARRAY_out(J)
         ARRAY_out(J) = ARRAY_out(I)
         ARRAY_out(I) = TEMP
         101 K = N / 2
         102 IF ( K .GE. J ) GO TO 103
              J = J - K
              K = K / 2
              GO TO 102
         103 J = J + K
   104 CONTINUE
END IF

! Divide inverse FFT output Coefficients by 1/N
ARRAY_out = ARRAY_out / CMPLX(REAL(N), 0., KIND=kind)

END SUBROUTINE ifft_sgl

SUBROUTINE ifft_dbl ( array_in, array_out, n, error )
!
!  Purpose:
!    A Cooley-Tukey Radix 2, in place, decimation-in-frequency
!    inverse FFT subroutine.  Although the routine is in place,
!    the input array is not destroyed because it is first copied
!    to the output array, and the FFT is performed on the copy
!    in the output array.
!
!  Record of revisions:
!       Date       Programmer          Description of change
!       ====       ==========          =====================
!     05/17/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=12) ! Precision
REAL(KIND=kind), PARAMETER :: twopi = 6.283185307179586_kind ! 2*pi

! Declare calling arguments:
INTEGER, INTENT(IN) :: n             ! No. of elements in IFFT--must
                                     ! be a pwer of 2:  2, 4, 8, etc.
COMPLEX(KIND=kind),INTENT(IN),DIMENSION(n) :: array_in
                                     ! Input data for IFFT
COMPLEX(KIND=kind),INTENT(INOUT),DIMENSION(n) :: array_out
                                     ! Output data for IFFT
INTEGER, INTENT(OUT) :: error        ! Error flag:
                                     ! 0 -- No errors
                                     ! 1 -- n not a power of 2.

! List of local variables:
REAL(KIND=kind) :: angle      ! Twiddle factor angle for butterfly
REAL(KIND=kind) :: anginc     ! Angle increment for twiddle factor
INTEGER :: i, j, k, l, n1, n2 ! Index variables
INTEGER :: m                  ! Number such that n = 2**m
COMPLEX(KIND=kind) :: temp    ! Temporary storage
COMPLEX(KIND=kind) :: w       ! Twiddle factor

! Copy the input array
array_out = array_in

! Check for proper array size.  The length of the input data
! array must be a power of 2.
j = 1
m = 0
1 IF ( j < n ) THEN
     M = M + 1
     J = J * 2
     GO TO 1
  END IF

! Make sure N is a power of 2.
IF ( J /= N ) THEN

   ! Not a power of 2--this is a error!
   error = 1

ELSE

   ! This is a power of 2.  Process the FFT!
   N2 = N

   ! Loop over the number of stages in the FFT.
   DO 10 K = 1, M
       N1     = N2
       N2     = N2 / 2
       ANGINC = TWOPI / N1
       ANGLE  = 0.

       ! Loop over all angles at the current stage.
       DO 20 J = 1, N2
          W = CMPLX ( COS(ANGLE), SIN(ANGLE), KIND=kind )
          ANGLE = J * ANGINC

          ! Perform butterflies.
          DO 30 I = J, N, N1
             L = I + N2
             TEMP         = ARRAY_out(I) - ARRAY_out(L)
             ARRAY_out(I) = ARRAY_out(I) + ARRAY_out(L)
             ARRAY_out(L) = TEMP * W
          30 CONTINUE
       20 CONTINUE
    10 CONTINUE

   ! Unscramble data with a digit reverse counter.
   J = 1
   N1 = N - 1
   DO 104 I = 1, N1
      IF ( I .GE. J ) GO TO 101
         TEMP   = ARRAY_out(J)
         ARRAY_out(J) = ARRAY_out(I)
         ARRAY_out(I) = TEMP
         101 K = N / 2
         102 IF ( K .GE. J ) GO TO 103
              J = J - K
              K = K / 2
              GO TO 102
         103 J = J + K
   104 CONTINUE
END IF

! Divide inverse FFT output Coefficients by 1/N
ARRAY_out = ARRAY_out / CMPLX(REAL(N), 0., KIND=kind)

END SUBROUTINE ifft_dbl

REAL FUNCTION random_u( )
!
!  Purpose:
!    Function to generate a uniform distribution on the range
!    0 <= random1 < 1.0 from a function.  It is based on the
!    intrinsic subroutine random_number.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    03/12/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Local variable:
REAL :: temp                   ! Temporary variable
CALL random_number( temp )
random_u = temp
END FUNCTION random_u

REAL FUNCTION random_n( )
!
!  Purpose:
!    Function to generate a gaussian normal distribution with
!    a mean of 0.0 and a standard deviation of 1.0.  It is
!    based on the Box-Muller method as documented in Chapter 7
!    of "Numerical Recipes".
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    03/12/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local variables.
LOGICAL,SAVE :: available = .FALSE. ! Value available?
REAL :: r                 ! SQRT(ran(1)**2+ran(2)**2)
REAL, DIMENSION(2) :: ran ! Unif random vars [0,1)
REAL, SAVE :: saved       ! Saved value
REAL, DIMENSION(2) :: v   ! Unif random vars [-1,1)
REAL, DIMENSION(2) :: y   ! Gaussian random vars

! Check to see if a saved value is available.
IF ( .NOT. available ) THEN

   ! No variable is available, so we must create new ones.
   ! Get 2 uniform random variables in the range [0.,1.) such
   ! that the square root of the sum of their squares < 1.
   ! Keep trying until we come up with such a combination.
   DO
      CALL random_number( ran )
      v = 2. * ran - 1.
      r = SUM( v**2 )
      IF ( r < 1. ) EXIT
   END DO

   ! Calculate the two Gaussian random from the uniform
   ! variables by Box-Muller method as:
   !   y = SQRT( -2.*LOG(r) / r ) * v
   y = SQRT( -2. * LOG(r) / r ) * v

   ! Return one value and save the other for next time.
   random_n = y(1)
   saved    = y(2)
   available = .TRUE.

ELSE ! A saved value was available.

   random_n  = saved
   available = .FALSE.

END IF

END FUNCTION random_n

REAL FUNCTION random_r( )
!
!  Purpose:
!    Function to generate a Rayleigh distribution.  It is
!    based on the Box-Muller method for generating normal
!    distributions.  This function generates a Rayleigh
!    distribution with a mean of 1.25 and a standard deviation
!    of:
!               sd = SQRT( 4./pi - 1. ) * mean
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    03/12/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local variables.
REAL :: r                 ! SQRT(ran(1)**2+ran(2)**2)
REAL, DIMENSION(2) :: ran ! Unif random vars [0,1)
REAL, DIMENSION(2) :: v   ! Unif random vars [-1,1)
REAL, DIMENSION(2) :: y   ! Gaussian random vars

! Get 2 uniform random variables in the range [0.,1.) such
! that the square root of the sum of their squares < 1.
! Keep trying until we come up with such a combination.
DO
   CALL random_number( ran )
   v = 2. * ran - 1.
   r = SUM( v**2 )
   IF ( r < 1. ) EXIT
END DO

! Calculate the two Gaussian random from the uniform
! variables by Box-Muller method as:
!   y = SQRT( -2.*LOG(r) / r ) * v
y = SQRT( -2. * LOG(r) / r ) * v

! Calculate a Rayleigh-distributed variable from them.
random_r = SQRT( y(1)**2 + y(2)**2 )

END FUNCTION random_r

SUBROUTINE interp_sgl ( x, y, npts, x0, y0, error )
!
!  Purpose:
!    To linearly interpolate the value y0 at position x0, given a
!    set of (x,y) measurements organized in increasing order of x.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    03/17/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=6) ! Precision

! Declare calling arguments:
INTEGER, INTENT(IN) :: npts           ! Dimension of arrays x and y
REAL(KIND=kind),DIMENSION(npts), INTENT(IN) :: x
                                      ! Independent variable x.
REAL(KIND=kind),DIMENSION(npts), INTENT(IN) :: y
                                      ! Dependent variable y.
REAL(KIND=kind),INTENT(IN) :: x0      ! Point to interpolate.
REAL(KIND=kind),INTENT(OUT) :: y0     ! Interpolated value.
INTEGER, INTENT(OUT) :: error         ! Error flag: 0 -- No error
                                      !            -1 -- x0 < x(1)
                                      !             1 -- x0 > x(npts)

! Declare local variables:
INTEGER :: i                     ! Index variable
INTEGER :: ibase                 ! Index for interpolation.
REAL :: slope                    ! Slope for interpolation.


! Assume that the input data set is in ascending order of x.
! See if the measurement position x0 is smaller or larger
! than any value in the data set.
IF ( x0 < x(1) ) THEN
   error = -1
ELSE IF ( x0 > x(npts) ) THEN
   error = 1
ELSE

   ! Point is between x(1) and x(npts).  Find the two points
   ! that it is between.
    DO i = 1, npts-1
       IF ( (x0 >= x(i)) .AND. (x0 <= x(i+1)) ) THEN
            ibase = i        ! Found the point
            EXIT
       END IF
    END DO

    ! Now linearly interpolate point.
    slope = ( y(ibase+1)-y(ibase) ) / ( x(ibase+1)-X(ibase) )
    y0 = slope * ( x0 - x(ibase) ) + y(ibase)

    ! Set error flag to 0.
    error = 0
END IF

END SUBROUTINE interp_sgl

SUBROUTINE interp_dbl ( x, y, npts, x0, y0, error )
!
!  Purpose:
!    To linearly interpolate the value y0 at position x0, given a
!    set of (x,y) measurements organized in increasing order of x.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    03/17/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=12) ! Precision

! Declare calling arguments:
INTEGER, INTENT(IN) :: npts           ! Dimension of arrays x and y
REAL(KIND=kind),DIMENSION(npts), INTENT(IN) :: x
                                      ! Independent variable x.
REAL(KIND=kind),DIMENSION(npts), INTENT(IN) :: y
                                      ! Dependent variable y.
REAL(KIND=kind),INTENT(IN) :: x0      ! Point to interpolate.
REAL(KIND=kind),INTENT(OUT) :: y0     ! Interpolated value.
INTEGER, INTENT(OUT) :: error         ! Error flag: 0 -- No error
                                      !            -1 -- x0 < x(1)
                                      !             1 -- x0 > x(npts)

! Declare local variables:
INTEGER :: i                     ! Index variable
INTEGER :: ibase                 ! Index for interpolation.
REAL :: slope                    ! Slope for interpolation.


! Assume that the input data set is in ascending order of x.
! See if the measurement position x0 is smaller or larger
! than any value in the data set.
IF ( x0 < x(1) ) THEN
   error = -1
ELSE IF ( x0 > x(npts) ) THEN
   error = 1
ELSE

   ! Point is between x(1) and x(npts).  Find the two points
   ! that it is between.
    DO i = 1, npts-1
       IF ( (x0 >= x(i)) .AND. (x0 <= x(i+1)) ) THEN
            ibase = i        ! Found the point
            EXIT
       END IF
    END DO

    ! Now linearly interpolate point.
    slope = ( y(ibase+1)-y(ibase) ) / ( x(ibase+1)-X(ibase) )
    y0 = slope * ( x0 - x(ibase) ) + y(ibase)

    ! Set error flag to 0.
    error = 0
END IF

END SUBROUTINE interp_dbl

SUBROUTINE lstmul ( value, base, exp, mul )
!
!  Purpose:
!    Subroutine to calculate the largest exponent EXP such that
!    value => mul ( = base**exp ).  This calculation is useful
!    for sizing FFT's, etc.
!
IMPLICIT NONE

! List of dummy arguments:
INTEGER,INTENT(IN) :: value ! Value to compare.
INTEGER,INTENT(IN) :: base  ! Base to which power EXP is raised.
INTEGER,INTENT(OUT) :: exp  ! The exponent.
INTEGER,INTENT(OUT) :: mul  ! The next power of 2: MUL = BASE**EXP.

! Calculate "mul" and power of "base" required.
exp = 0
mul = base ** exp

! Find the power of base which exceeds or equals value.
DO
   exp = exp + 1
   mul = base ** exp
   IF ( mul >= value ) EXIT
END DO

! If the pwer exceeds value, back off by one power of base.
IF ( mul > value ) THEN
   exp = exp - 1
   mul = base ** exp
END IF

END SUBROUTINE lstmul

SUBROUTINE mat_inv_sgl ( a, b, ndim, n, error )
!
!  Purpose:
!    Subroutine to n x n matrix using Gaussian elimination
!    and the maximum pivot technique.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    03/16/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=6) ! Precision
REAL(KIND=kind), PARAMETER :: eps = 10. * EPSILON(0.0_kind)
                                     ! A "small" number for comparison
                                     ! when determining singular eqns

! Declare calling arguments:
INTEGER, INTENT(IN) :: ndim          ! Dimension of arrays a and b
REAL(KIND=kind),INTENT(IN),DIMENSION(ndim,ndim) :: a
                                     ! Input matrix (N x N).  This
                                     ! array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
                                     ! The declared dimension ndim
                                     ! must be passed to the sub, or
                                     ! it won't be able to interpret
                                     ! subscripts correctly.
REAL(KIND=kind),INTENT(OUT),DIMENSION(ndim,ndim) :: b
                                     ! Inverse of matrix a.  This
                                     ! array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
                                     ! The declared dimension ndim
                                     ! must be passed to the sub, or
                                     ! it won't be able to interpret
                                     ! subscripts correctly.
INTEGER, INTENT(IN) :: n             ! Number of equations to solve.
INTEGER, INTENT(OUT) :: error        ! Error flag:
                                     !   0 -- No error
                                     !   1 -- Singular equations

! Declare local variables:
REAL(KIND=kind),DIMENSION(n,n) :: a1 ! Copy of a to destroy while
                                     ! building the inverse.  Only
                                     ! actual elements in use are
                                     ! duplicated.
REAL(KIND=kind),DIMENSION(n,n) :: b1 ! Array in which to build
                                     ! inverse.  Only acutal elements
                                     ! in use are duplicated.
REAL(KIND=kind) :: factor            ! Factor to multiply eqn irow by
                                     ! before adding to eqn jrow
INTEGER :: irow                      ! Number of the equation currently
                                     ! currently being processed
INTEGER :: ipeak                     ! Pointer to equation containing
                                     ! maximum pivot value
INTEGER :: jrow                      ! Number of the equation compared
                                     ! to the current equation
REAL(KIND=kind),DIMENSION(n) :: temp ! Scratch array

! Make a copy of the input array.
a1 = a(1:n,1:n)

! Initialize the inverse array.
b1 = 0.
DO irow = 1, n
   b1(irow,irow) = 1.
END DO

! Process n times to get all equations...
mainloop: DO irow = 1, n

   ! Find peak pivot for column irow in rows irow to n
   ipeak = irow
   max_pivot: DO jrow = irow+1, n
      IF (ABS(a1(jrow,irow)) > ABS(a1(ipeak,irow))) THEN
         ipeak = jrow
      END IF
   END DO max_pivot

   ! Check for singular equations.
   singular: IF ( ABS(a1(ipeak,irow)) < eps ) THEN
      error = 1
      RETURN
   END IF singular

   ! Otherwise, if ipeak /= irow, swap equations irow & ipeak
   swap_eqn: IF ( ipeak /= irow ) THEN
      temp = a1(ipeak,:)
      a1(ipeak,:) = a1(irow,:)     ! Swap rows in a1
      a1(irow,:) = temp
      temp = b1(ipeak,:)
      b1(ipeak,:) = b1(irow,:)     ! Swap rows in b1
      b1(irow,:) = temp
   END IF swap_eqn

   ! Multiply equation irow by -a1(jrow,irow)/a1(irow,irow), and
   ! add it to Eqn jrow (for all eqns except irow itself).
   eliminate: DO jrow = 1, n
      IF ( jrow /= irow ) THEN
         factor = -a1(jrow,irow)/a1(irow,irow)
         a1(jrow,:) = a1(irow,:)*factor + a1(jrow,:)
         b1(jrow,:) = b1(irow,:)*factor + b1(jrow,:)
      END IF
   END DO eliminate
END DO mainloop

! End of main loop over all equations.  All off-diagonal
! terms are now zero.  To get the final answer, we must
! divide each equation by the coefficient of its on-diagonal
! term.
divide: DO irow = 1, n
   b1(irow,:) = b1(irow,:) / a1(irow,irow)
   a1(irow,irow) = 1.
END DO divide

! Copy the answer to the output, set error flag
! to 0 and return.
b = 0.
b(1:n,1:n) = b1
error = 0

END SUBROUTINE mat_inv_sgl

SUBROUTINE mat_inv_dbl ( a, b, ndim, n, error )
!
!  Purpose:
!    Subroutine to n x n matrix using Gaussian elimination
!    and the maximum pivot technique.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    03/16/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=12) ! Precision
REAL(KIND=kind), PARAMETER :: eps = 10. * EPSILON(0.0_kind)
                                     ! A "small" number for comparison
                                     ! when determining singular eqns

! Declare calling arguments:
INTEGER, INTENT(IN) :: ndim          ! Dimension of arrays a and b
REAL(KIND=kind),INTENT(IN),DIMENSION(ndim,ndim) :: a
                                     ! Input matrix (N x N).  This
                                     ! array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
                                     ! The declared dimension ndim
                                     ! must be passed to the sub, or
                                     ! it won't be able to interpret
                                     ! subscripts correctly.
REAL(KIND=kind),INTENT(OUT),DIMENSION(ndim,ndim) :: b
                                     ! Inverse of matrix a.  This
                                     ! array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
                                     ! The declared dimension ndim
                                     ! must be passed to the sub, or
                                     ! it won't be able to interpret
                                     ! subscripts correctly.
INTEGER, INTENT(IN) :: n             ! Number of equations to solve.
INTEGER, INTENT(OUT) :: error        ! Error flag:
                                     !   0 -- No error
                                     !   1 -- Singular equations

! Declare local variables:
REAL(KIND=kind),DIMENSION(n,n) :: a1 ! Copy of a to destroy while
                                     ! building the inverse.  Only
                                     ! actual elements in use are
                                     ! duplicated.
REAL(KIND=kind),DIMENSION(n,n) :: b1 ! Array in which to build
                                     ! inverse.  Only acutal elements
                                     ! in use are duplicated.
REAL(KIND=kind) :: factor            ! Factor to multiply eqn irow by
                                     ! before adding to eqn jrow
INTEGER :: irow                      ! Number of the equation currently
                                     ! currently being processed
INTEGER :: ipeak                     ! Pointer to equation containing
                                     ! maximum pivot value
INTEGER :: jrow                      ! Number of the equation compared
                                     ! to the current equation
REAL(KIND=kind),DIMENSION(n) :: temp ! Scratch array

! Make a copy of the input array.
a1 = a(1:n,1:n)

! Initialize the inverse array.
b1 = 0.
DO irow = 1, n
   b1(irow,irow) = 1.
END DO

! Process n times to get all equations...
mainloop: DO irow = 1, n

   ! Find peak pivot for column irow in rows irow to n
   ipeak = irow
   max_pivot: DO jrow = irow+1, n
      IF (ABS(a1(jrow,irow)) > ABS(a1(ipeak,irow))) THEN
         ipeak = jrow
      END IF
   END DO max_pivot

   ! Check for singular equations.
   singular: IF ( ABS(a1(ipeak,irow)) < eps ) THEN
      error = 1
      RETURN
   END IF singular

   ! Otherwise, if ipeak /= irow, swap equations irow & ipeak
   swap_eqn: IF ( ipeak /= irow ) THEN
      temp = a1(ipeak,:)
      a1(ipeak,:) = a1(irow,:)     ! Swap rows in a1
      a1(irow,:) = temp
      temp = b1(ipeak,:)
      b1(ipeak,:) = b1(irow,:)     ! Swap rows in b1
      b1(irow,:) = temp
   END IF swap_eqn

   ! Multiply equation irow by -a1(jrow,irow)/a1(irow,irow), and
   ! add it to Eqn jrow (for all eqns except irow itself).
   eliminate: DO jrow = 1, n
      IF ( jrow /= irow ) THEN
         factor = -a1(jrow,irow)/a1(irow,irow)
         a1(jrow,:) = a1(irow,:)*factor + a1(jrow,:)
         b1(jrow,:) = b1(irow,:)*factor + b1(jrow,:)
      END IF
   END DO eliminate
END DO mainloop

! End of main loop over all equations.  All off-diagonal
! terms are now zero.  To get the final answer, we must
! divide each equation by the coefficient of its on-diagonal
! term.
divide: DO irow = 1, n
   b1(irow,:) = b1(irow,:) / a1(irow,irow)
   a1(irow,irow) = 1.
END DO divide

! Copy the answer to the output, set error flag
! to 0 and return.
b = 0.
b(1:n,1:n) = b1
error = 0

END SUBROUTINE mat_inv_dbl

SUBROUTINE mat_inv_sgl_cmplx ( a, b, ndim, n, error )
!
!  Purpose:
!    Subroutine to n x n matrix using Gaussian elimination
!    and the maximum pivot technique.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    03/16/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=6) ! Precision
REAL(KIND=kind), PARAMETER :: eps = 10. * EPSILON(0.0_kind)
                                     ! A "small" number for comparison
                                     ! when determining singular eqns

! Declare calling arguments:
INTEGER, INTENT(IN) :: ndim          ! Dimension of arrays a and b
COMPLEX(KIND=kind),INTENT(IN),DIMENSION(ndim,ndim) :: a
                                     ! Input matrix (N x N).  This
                                     ! array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
                                     ! The declared dimension ndim
                                     ! must be passed to the sub, or
                                     ! it won't be able to interpret
                                     ! subscripts correctly.
COMPLEX(KIND=kind),INTENT(OUT),DIMENSION(ndim,ndim) :: b
                                     ! Inverse of matrix a.  This
                                     ! array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
                                     ! The declared dimension ndim
                                     ! must be passed to the sub, or
                                     ! it won't be able to interpret
                                     ! subscripts correctly.
INTEGER, INTENT(IN) :: n             ! Number of equations to solve.
INTEGER, INTENT(OUT) :: error        ! Error flag:
                                     !   0 -- No error
                                     !   1 -- Singular equations

! Declare local variables:
COMPLEX(KIND=kind),DIMENSION(n,n) :: a1 ! Copy of a to destroy while
                                     ! building the inverse.  Only
                                     ! actual elements in use are
                                     ! duplicated.
COMPLEX(KIND=kind),DIMENSION(n,n) :: b1 ! Array in which to build
                                     ! inverse.  Only acutal elements
                                     ! in use are duplicated.
COMPLEX(KIND=kind) :: factor         ! Factor to multiply eqn irow by
                                     ! before adding to eqn jrow
INTEGER :: irow                      ! Number of the equation currently
                                     ! currently being processed
INTEGER :: ipeak                     ! Pointer to equation containing
                                     ! maximum pivot value
INTEGER :: jrow                      ! Number of the equation compared
                                     ! to the current equation
COMPLEX(KIND=kind),DIMENSION(n) :: temp ! Scratch array

! Make a copy of the input array.
a1 = a(1:n,1:n)

! Initialize the inverse array.
b1 = 0.
DO irow = 1, n
   b1(irow,irow) = (1.0_kind,0.0_kind)
END DO

! Process n times to get all equations...
mainloop: DO irow = 1, n

   ! Find peak pivot for column irow in rows irow to n
   ipeak = irow
   max_pivot: DO jrow = irow+1, n
      IF (ABS(a1(jrow,irow)) > ABS(a1(ipeak,irow))) THEN
         ipeak = jrow
      END IF
   END DO max_pivot

   ! Check for singular equations.
   singular: IF ( ABS(a1(ipeak,irow)) < eps ) THEN
      error = 1
      RETURN
   END IF singular

   ! Otherwise, if ipeak /= irow, swap equations irow & ipeak
   swap_eqn: IF ( ipeak /= irow ) THEN
      temp = a1(ipeak,:)
      a1(ipeak,:) = a1(irow,:)     ! Swap rows in a1
      a1(irow,:) = temp
      temp = b1(ipeak,:)
      b1(ipeak,:) = b1(irow,:)     ! Swap rows in b1
      b1(irow,:) = temp
   END IF swap_eqn

   ! Multiply equation irow by -a1(jrow,irow)/a1(irow,irow), and
   ! add it to Eqn jrow (for all eqns except irow itself).
   eliminate: DO jrow = 1, n
      IF ( jrow /= irow ) THEN
         factor = -a1(jrow,irow)/a1(irow,irow)
         a1(jrow,:) = a1(irow,:)*factor + a1(jrow,:)
         b1(jrow,:) = b1(irow,:)*factor + b1(jrow,:)
      END IF
   END DO eliminate
END DO mainloop

! End of main loop over all equations.  All off-diagonal
! terms are now zero.  To get the final answer, we must
! divide each equation by the coefficient of its on-diagonal
! term.
divide: DO irow = 1, n
   b1(irow,:) = b1(irow,:) / a1(irow,irow)
   a1(irow,irow) = 1.
END DO divide

! Copy the answer to the output, set error flag
! to 0 and return.
b = 0.
b(1:n,1:n) = b1
error = 0

END SUBROUTINE mat_inv_sgl_cmplx

SUBROUTINE mat_inv_dbl_cmplx ( a, b, ndim, n, error )
!
!  Purpose:
!    Subroutine to n x n matrix using Gaussian elimination
!    and the maximum pivot technique.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    03/16/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=12) ! Precision
REAL(KIND=kind), PARAMETER :: eps = 10. * EPSILON(0.0_kind)
                                     ! A "small" number for comparison
                                     ! when determining singular eqns

! Declare calling arguments:
INTEGER, INTENT(IN) :: ndim          ! Dimension of arrays a and b
COMPLEX(KIND=kind),INTENT(IN),DIMENSION(ndim,ndim) :: a
                                     ! Input matrix (N x N).  This
                                     ! array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
                                     ! The declared dimension ndim
                                     ! must be passed to the sub, or
                                     ! it won't be able to interpret
                                     ! subscripts correctly.
COMPLEX(KIND=kind),INTENT(OUT),DIMENSION(ndim,ndim) :: b
                                     ! Inverse of matrix a.  This
                                     ! array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
                                     ! The declared dimension ndim
                                     ! must be passed to the sub, or
                                     ! it won't be able to interpret
                                     ! subscripts correctly.
INTEGER, INTENT(IN) :: n             ! Number of equations to solve.
INTEGER, INTENT(OUT) :: error        ! Error flag:
                                     !   0 -- No error
                                     !   1 -- Singular equations

! Declare local variables:
COMPLEX(KIND=kind),DIMENSION(n,n) :: a1 ! Copy of a to destroy while
                                     ! building the inverse.  Only
                                     ! actual elements in use are
                                     ! duplicated.
COMPLEX(KIND=kind),DIMENSION(n,n) :: b1 ! Array in which to build
                                     ! inverse.  Only acutal elements
                                     ! in use are duplicated.
COMPLEX(KIND=kind) :: factor         ! Factor to multiply eqn irow by
                                     ! before adding to eqn jrow
INTEGER :: irow                      ! Number of the equation currently
                                     ! currently being processed
INTEGER :: ipeak                     ! Pointer to equation containing
                                     ! maximum pivot value
INTEGER :: jrow                      ! Number of the equation compared
                                     ! to the current equation
COMPLEX(KIND=kind),DIMENSION(n) :: temp ! Scratch array

! Make a copy of the input array.
a1 = a(1:n,1:n)

! Initialize the inverse array.
b1 = 0.
DO irow = 1, n
   b1(irow,irow) = (1.0_kind,0.0_kind)
END DO

! Process n times to get all equations...
mainloop: DO irow = 1, n

   ! Find peak pivot for column irow in rows irow to n
   ipeak = irow
   max_pivot: DO jrow = irow+1, n
      IF (ABS(a1(jrow,irow)) > ABS(a1(ipeak,irow))) THEN
         ipeak = jrow
      END IF
   END DO max_pivot

   ! Check for singular equations.
   singular: IF ( ABS(a1(ipeak,irow)) < eps ) THEN
      error = 1
      RETURN
   END IF singular

   ! Otherwise, if ipeak /= irow, swap equations irow & ipeak
   swap_eqn: IF ( ipeak /= irow ) THEN
      temp = a1(ipeak,:)
      a1(ipeak,:) = a1(irow,:)     ! Swap rows in a1
      a1(irow,:) = temp
      temp = b1(ipeak,:)
      b1(ipeak,:) = b1(irow,:)     ! Swap rows in b1
      b1(irow,:) = temp
   END IF swap_eqn

   ! Multiply equation irow by -a1(jrow,irow)/a1(irow,irow), and
   ! add it to Eqn jrow (for all eqns except irow itself).
   eliminate: DO jrow = 1, n
      IF ( jrow /= irow ) THEN
         factor = -a1(jrow,irow)/a1(irow,irow)
         a1(jrow,:) = a1(irow,:)*factor + a1(jrow,:)
         b1(jrow,:) = b1(irow,:)*factor + b1(jrow,:)
      END IF
   END DO eliminate
END DO mainloop

! End of main loop over all equations.  All off-diagonal
! terms are now zero.  To get the final answer, we must
! divide each equation by the coefficient of its on-diagonal
! term.
divide: DO irow = 1, n
   b1(irow,:) = b1(irow,:) / a1(irow,irow)
   a1(irow,irow) = 1.
END DO divide

! Copy the answer to the output, set error flag
! to 0 and return.
b = 0.
b(1:n,1:n) = b1
error = 0

END SUBROUTINE mat_inv_dbl_cmplx

SUBROUTINE nxtmul ( value, base, exp, mul )
!
!  Purpose:
!    Subroutine to calculate the smallest exponent EXP such that
!    value <= mul ( = base**exp ).  This calculation is useful
!    for sizing FFT's, etc.
!
IMPLICIT NONE

! List of dummy arguments:
INTEGER,INTENT(IN) :: value ! Value to compare.
INTEGER,INTENT(IN) :: base  ! Base to which power EXP is raised.
INTEGER,INTENT(OUT) :: exp  ! The exponent.
INTEGER,INTENT(OUT) :: mul  ! The next power of 2: MUL = BASE**EXP.

! Calculate "mul" and power of "base" required.
exp = 0
mul = base ** exp

! Find the power of base which exceeds or equals value.
DO
   exp = exp + 1
   mul = base ** exp
   IF ( mul >= value ) EXIT
END DO

END SUBROUTINE nxtmul

SUBROUTINE plot_sgl ( y, npts, lu, minamp, maxamp )
!
!  Purpose:
!    Subroutine to plot the points in array y.  The data in the
!    array is assumed to be at a uniform spacing.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    07/01/95    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=6) ! Precision
INTEGER, PARAMETER :: nbins = 65   ! Number of plotting bins

! List of calling arguments:
INTEGER, INTENT(IN) :: npts        ! Number of points to plot
REAL(KIND=kind), DIMENSION(npts), INTENT(IN) :: y
                                   ! Data to plot.
INTEGER, INTENT(IN) :: lu          ! I/o unit to plot on.
REAL(KIND=kind), INTENT(IN), OPTIONAL :: minamp
                                   ! Smallest value to plot
REAL(KIND=kind), INTENT(IN), OPTIONAL :: maxamp
                                   ! Largest value to plot

! List of local variables:
CHARACTER(len=14) :: annot         ! annotation (y-amplitude)
INTEGER :: i                       ! Loop index
INTEGER :: ibin                    ! Bin no. for current Y value
INTEGER :: ibin0                   ! Bin no. for zero-crossing
REAL(KIND=kind) :: frac            ! frac of plot width to position data
CHARACTER(len=79) :: line          ! Output buffer
REAL(KIND=kind) :: maxvl1          ! Largest value to plot
REAL(KIND=kind) :: minvl1          ! Smallest value to plot
CHARACTER(len=65) :: pltbuf        ! Plotting buffer
CHARACTER(len=65) :: scl           ! Scale on border of plot

! Set border
scl( 1:30) = '+-----------------------------'
scl(31:60) = '------------------------------'
scl(61:65) = '----+'

! If the scales are defaulted, set min and max of y axis.
IF ( PRESENT(minamp) ) THEN
   minvl1 = minamp
ELSE
   minvl1 = MINVAL(y)
END IF

IF ( PRESENT(maxamp) ) THEN
   maxvl1 = maxamp
ELSE
   maxvl1 = MAXVAL(y)
END IF

! We will divide minvl1 to maxvl1 into 65 bins for
! plotting purposes.  Locate the zero bin, if it is
! between minvl1 and maxvl1.
!
IF ( (maxvl1 > 0.) .AND. (minvl1 < 0) ) THEN
   frac = ( 0. - minvl1) / (maxvl1 - minvl1 )
   ibin0 = NINT ( (nbins-1) * frac ) + 1
ELSE
   ibin0 = 0
END IF

! Set zero into border of plot, if it is within
! the limits of the plot.
annot = ' '
IF ( ibin0 > 0 ) THEN
   scl(ibin0:ibin0) = '+'
END IF

! Print upper border.
line = ' '
line(10:21) = real_to_char(minvl1)
line(68:79) = real_to_char(maxvl1)
WRITE (lu,'(A)' ) line
WRITE (lu,'(A,1X,A)') annot, scl

! Plot data points.
DO i = 1, npts

   ! Clear line.
   pltbuf = ' '
   annot  = ' '

   ! Set value of y data point.
   annot(2:13) = real_to_char( y(i) )

   ! Set min and max borders.
   pltbuf(1:1)   = '|'
   pltbuf(65:65) = '|'

   ! Set zero line, if within borders.
   IF ( ibin0 > 0 ) THEN
      pltbuf(ibin0:ibin0) = '|'
   END IF

   ! Plot point on array.
   frac = ( y(i) - minvl1) / (maxvl1 - minvl1 )
   ibin = NINT ( (nbins-1) * frac ) + 1
   IF ( (ibin >= 1) .AND. (ibin <= nbins) ) THEN
      pltbuf(ibin:ibin) = '*'
   END IF

   !Write out line.
   WRITE (lu,'(A,1X,A)') annot, pltbuf

END DO

! Print amplitude scale at bottom of plot.
annot = ' '
WRITE (lu,'(A,1X,A)') annot, scl
WRITE (lu,'(A)' ) line

!     Print out summary info.
WRITE (lu,'(/,10X,A,I12)' ) 'Number of Points = ', npts

CONTAINS
   FUNCTION real_to_char ( value )
   !
   !  Purpose:
   !    To convert a real value into a 12-character string, with the
   !    number printed in as readable a format as possible considering
   !    its range.  This routine prints out the number according to the
   !    following rules:
   !       1.  value > 9999999.                ES12.5
   !       2.  value < -999999.                ES12.5
   !       3.  0.    <  ABS(value) < 0.01      ES12.5
   !       5.  value = 0.0                     F12.4
   !       6.  Otherwise                       F12.4
   !
   !  Record of revisions:
   !      Date       Programmer          Description of change
   !      ====       ==========          =====================
   !    11/26/95    S. J. Chapman        Original code
   !
   IMPLICIT NONE

   ! Declare calling arguments:
   REAL(KIND=kind), INTENT(IN) :: value ! value to convert to char form
   CHARACTER (len=12) :: real_to_char  ! Output character string

   ! Declare local variables:
   CHARACTER(len=9) :: fmt             ! Format descriptor
   CHARACTER(len=12) :: string = ' '   ! Output string

   ! Select proper format
   IF ( value > 9999999. ) THEN
      fmt = '(ES12.5)'
   ELSE IF ( value < -999999. ) THEN
      fmt = '(ES12.5)'
   ELSE IF ( value == 0. ) THEN
      fmt = '(F12.4)'
   ELSE IF ( ABS(value) < 0.01 ) THEN
      fmt = '(ES12.5)'
   ELSE
      fmt = '(F12.4)'
   END IF

   ! Convert to character form.
   WRITE (string,fmt) value
   real_to_char = string

   END FUNCTION real_to_char

END SUBROUTINE plot_sgl

SUBROUTINE plot_dbl ( y, npts, lu, minamp, maxamp )
!
!  Purpose:
!    Subroutine to plot the points in array y.  The data in the
!    array is assumed to be at a uniform spacing.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    07/01/95    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=12) ! Precision
INTEGER, PARAMETER :: nbins = 65   ! Number of plotting bins

! List of calling arguments:
INTEGER, INTENT(IN) :: npts        ! Number of points to plot
REAL(KIND=kind), DIMENSION(npts), INTENT(IN) :: y
                                   ! Data to plot.
INTEGER, INTENT(IN) :: lu          ! I/o unit to plot on.
REAL(KIND=kind), INTENT(IN), OPTIONAL :: minamp
                                   ! Smallest value to plot
REAL(KIND=kind), INTENT(IN), OPTIONAL :: maxamp
                                   ! Largest value to plot

! List of local variables:
CHARACTER(len=14) :: annot         ! annotation (y-amplitude)
INTEGER :: i                       ! Loop index
INTEGER :: ibin                    ! Bin no. for current Y value
INTEGER :: ibin0                   ! Bin no. for zero-crossing
REAL(KIND=kind) :: frac            ! frac of plot width to position data
CHARACTER(len=79) :: line          ! Output buffer
REAL(KIND=kind) :: maxvl1          ! Largest value to plot
REAL(KIND=kind) :: minvl1          ! Smallest value to plot
CHARACTER(len=65) :: pltbuf        ! Plotting buffer
CHARACTER(len=65) :: scl           ! Scale on border of plot

! Set border
scl( 1:30) = '+-----------------------------'
scl(31:60) = '------------------------------'
scl(61:65) = '----+'

! If the scales are defaulted, set min and max of y axis.
IF ( PRESENT(minamp) ) THEN
   minvl1 = minamp
ELSE
   minvl1 = MINVAL(y)
END IF

IF ( PRESENT(maxamp) ) THEN
   maxvl1 = maxamp
ELSE
   maxvl1 = MAXVAL(y)
END IF

! We will divide minvl1 to maxvl1 into 65 bins for
! plotting purposes.  Locate the zero bin, if it is
! between minvl1 and maxvl1.
!
IF ( (maxvl1 > 0.) .AND. (minvl1 < 0) ) THEN
   frac = ( 0. - minvl1) / (maxvl1 - minvl1 )
   ibin0 = NINT ( (nbins-1) * frac ) + 1
ELSE
   ibin0 = 0
END IF

! Set zero into border of plot, if it is within
! the limits of the plot.
annot = ' '
IF ( ibin0 > 0 ) THEN
   scl(ibin0:ibin0) = '+'
END IF

! Print upper border.
line = ' '
line(10:21) = real_to_char(minvl1)
line(68:79) = real_to_char(maxvl1)
WRITE (lu,'(A)' ) line
WRITE (lu,'(A,1X,A)') annot, scl

! Plot data points.
DO i = 1, npts

   ! Clear line.
   pltbuf = ' '
   annot  = ' '

   ! Set value of y data point.
   annot(2:13) = real_to_char( y(i) )

   ! Set min and max borders.
   pltbuf(1:1)   = '|'
   pltbuf(65:65) = '|'

   ! Set zero line, if within borders.
   IF ( ibin0 > 0 ) THEN
      pltbuf(ibin0:ibin0) = '|'
   END IF

   ! Plot point on array.
   frac = ( y(i) - minvl1) / (maxvl1 - minvl1 )
   ibin = NINT ( (nbins-1) * frac ) + 1
   IF ( (ibin >= 1) .AND. (ibin <= nbins) ) THEN
      pltbuf(ibin:ibin) = '*'
   END IF

   !Write out line.
   WRITE (lu,'(A,1X,A)') annot, pltbuf

END DO

! Print amplitude scale at bottom of plot.
annot = ' '
WRITE (lu,'(A,1X,A)') annot, scl
WRITE (lu,'(A)' ) line

!     Print out summary info.
WRITE (lu,'(/,10X,A,I12)' ) 'Number of Points = ', npts

CONTAINS
   FUNCTION real_to_char ( value )
   !
   !  Purpose:
   !    To convert a real value into a 12-character string, with the
   !    number printed in as readable a format as possible considering
   !    its range.  This routine prints out the number according to the
   !    following rules:
   !       1.  value > 9999999.                ES12.5
   !       2.  value < -999999.                ES12.5
   !       3.  0.    <  ABS(value) < 0.01      ES12.5
   !       5.  value = 0.0                     F12.4
   !       6.  Otherwise                       F12.4
   !
   !  Record of revisions:
   !      Date       Programmer          Description of change
   !      ====       ==========          =====================
   !    11/26/95    S. J. Chapman        Original code
   !
   IMPLICIT NONE

   ! Declare calling arguments:
   REAL(KIND=kind), INTENT(IN) :: value ! value to convert to char form
   CHARACTER (len=12) :: real_to_char  ! Output character string

   ! Declare local variables:
   CHARACTER(len=9) :: fmt             ! Format descriptor
   CHARACTER(len=12) :: string = ' '   ! Output string

   ! Select proper format
   IF ( value > 9999999. ) THEN
      fmt = '(ES12.5)'
   ELSE IF ( value < -999999. ) THEN
      fmt = '(ES12.5)'
   ELSE IF ( value == 0. ) THEN
      fmt = '(F12.4)'
   ELSE IF ( ABS(value) < 0.01 ) THEN
      fmt = '(ES12.5)'
   ELSE
      fmt = '(F12.4)'
   END IF

   ! Convert to character form.
   WRITE (string,fmt) value
   real_to_char = string

   END FUNCTION real_to_char

END SUBROUTINE plot_dbl

SUBROUTINE plotxy_sgl ( x, y, npts, lu, minx, maxx, miny, maxy, &
                        nbinx, nbiny )
!
!  Purpose:
!    Subroutine to plot (x,y) pairs of points.  Note that
!    x and y scales may either be set or defaulted on an
!    item-by-item basis.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    03/20/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=6) ! Precision
INTEGER, PARAMETER :: maxbin = 65           ! Max num. of bins

! List of dummy arguments:
INTEGER, INTENT(IN) :: npts                 ! No. points to plot
REAL(KIND=kind),INTENT(IN),DIMENSION(npts) :: x  ! x data
REAL(KIND=kind),INTENT(IN),DIMENSION(npts) :: y  ! y data
INTEGER, INTENT(IN) :: lu                   ! Output unit
REAL(KIND=kind),INTENT(IN),OPTIONAL :: minx ! Smallest x to plot
REAL(KIND=kind),INTENT(IN),OPTIONAL :: maxx ! Largest x to plot
REAL(KIND=kind),INTENT(IN),OPTIONAL :: miny ! Smallest y to plot
REAL(KIND=kind),INTENT(IN),OPTIONAL :: maxy ! Largest y to plot
INTEGER,INTENT(IN),OPTIONAL :: nbinx        ! No. bins for x
INTEGER,INTENT(IN),OPTIONAL :: nbiny        ! No. bins for y

! List of local variables:
CHARACTER(len=14) :: annot    ! Annotation (y-amplitude)
REAL(KIND=kind) :: frac       ! Fraction of plot width
INTEGER :: i                  ! Loop index
!INTEGER :: ibinx0            ! Bin for x zero crossing
INTEGER :: ibiny              ! Bin for current y value
INTEGER :: ibiny0             ! Bin for y zero crossing
INTEGER :: j                  ! Loop index
REAL(KIND=kind) :: maxx1      ! Largest x to plot
REAL(KIND=kind) :: minx1      ! Samllest x to plot
REAL(KIND=kind) :: maxy1      ! Largest y to plot
REAL(KIND=kind) :: miny1      ! Samllest y to plot
INTEGER :: nbinx1             ! No. of x bins to plot
INTEGER :: nbiny1             ! No. of y bins to plot
CHARACTER(len=65) :: pltbuf   ! Plotting buffer
CHARACTER(len=65) :: scl      ! Scale on plot border
CHARACTER(len=65),SAVE:: scl1 ! Template for scale
CHARACTER(len=80) :: sclbuf   ! Plotting scale buffer
REAL(KIND=kind) :: xpm        ! Beginning of x range (cur x pos)
REAL(KIND=kind) :: xpos       ! Center of x range (cur x pos)
REAL(KIND=kind) :: xpp        ! Ending of x range (cur x pos)

! Initialize scale border.
scl1 = '+-----------------------------&
       &-----------------------------------'

!  If any scale values are defaulted, set min and max of x and y
!  axes according to the input data.  Note that we are making
!  local copies of minx1, maxx1, miny1, and maxy1 if the exist,
!  so that we may modify their values.  Check for and set defaults
! on all optional arguments.
!
IF ( PRESENT(minx) ) THEN
   minx1 = minx
ELSE
   minx1 = minval(x)
END IF

IF ( PRESENT(maxx) ) THEN
   maxx1 = maxx
ELSE
   maxx1 = maxval(x)
END IF

IF ( PRESENT(miny) ) THEN
   miny1 = miny
ELSE
   miny1 = minval(y)
END IF

IF ( PRESENT(maxy) ) THEN
   maxy1 = maxy
ELSE
   maxy1 = maxval(y)
END IF

IF ( PRESENT(nbinx) ) THEN
   nbinx1 = MAX(1, MIN(nbinx,maxbin) )
ELSE
   nbinx1 = maxbin
END IF

IF ( PRESENT(nbiny) ) THEN
   nbiny1 = MAX(1, MIN(nbiny,maxbin) )
ELSE
   nbiny1 = maxbin
END IF

! We will divide minx1 to maxx1 into nbinx1 bins for
! plotting purposes.  Locate the zero bin, if it is
! between minx1 and maxx1.
!IF ( (maxx1 > 0.) .AND. (minx1 < 0.) ) THEN
!   frac = ( 0. - minx1) / (maxx1 - minx1 )
!   ibinx0 = NINT ( (nbinx1-1) * frac ) + 1
!ELSE
!   ibinx0 = 0
!END IF

! We will divide miny1 to maxy1 into nbiny1 bins for
! plotting purposes.  Locate the zero bin, if it is
! between miny1 and maxy1.
 IF ( (maxy1 > 0.) .AND. (miny1 < 0.) ) THEN
   frac = ( 0. - miny1) / (maxy1 - miny1 )
   ibiny0 = NINT ( (nbiny1-1) * frac ) + 1
ELSE
   ibiny0 = 0
END IF

! Set upper right hand corner of plot.
annot = ' '
scl   = scl1(1:nbiny1)
scl(nbiny1:nbiny1) = '+'

! Set y axis zero into border of plot, if it is within
! the limits of the plot.
IF ( ibiny0 > 0 ) THEN
   scl(ibiny0:ibiny0) = '+'
END IF

! Generate plot scale label.
sclbuf = ' '
sclbuf(10:22) = real_to_char(miny1)
sclbuf(15+nbiny1-11:15+nbiny1) = real_to_char(maxy1)

! Set y axis zero into scale of plot, if it is within
! the limits of the plot.
IF ( ibiny0 > 0 ) THEN
   sclbuf(15+ibiny0:15+ibiny0) = '0'
END IF

! Print upper border.
WRITE (lu,'(A)') sclbuf
WRITE (lu,'(A,1x,A)') annot, scl

! Plot data points.
plot: DO i = 1, nbinx1

   ! Calculate x position of this data point.
   xpos = ((REAL(I)-1.0)/(nbinx1-1))*(maxx1-minx1) + minx1
   xpm  = ((REAL(I)-1.5)/(nbinx1-1))*(maxx1-minx1) + minx1
   xpp  = ((REAL(I)-0.5)/(nbinx1-1))*(maxx1-minx1) + minx1

   ! Clear line.
   pltbuf = ' '
   annot  = ' '

   ! Set value of x data point.
   annot(2:13) = real_to_char( xpos )

   ! Set min and max borders.
   pltbuf(1:1)   = '|'
   pltbuf(nbiny1:nbiny1) = '|'

   ! Set zero line, if within borders.
   IF ( ibiny0 > 0 ) THEN
      pltbuf(ibiny0:ibiny0) = '|'
   END IF

   ! Find points that map to this x position.
   xpts: DO j = 1, npts
      IF ( (x(j) > xpm) .AND. (x(j) < xpp) ) THEN

         ! Plot point on array.
         frac  = ( y(J) - miny1) / (maxy1 - miny1 )
         ibiny = NINT ( (nbiny1-1) * frac ) + 1
         IF ( (ibiny >= 1) .AND. (ibiny <= nbiny1) ) THEN
            pltbuf(ibiny:ibiny) = '*'
         END IF
      END IF
   END DO xpts

   ! Output line.
   WRITE (lu,'(A,1X,A)') annot, pltbuf

END DO plot

! Print lower border.
annot = ' '
WRITE (lu,'(A,1X,A)') annot, scl
WRITE (lu,'(A)')      sclbuf

! Print out summary info.
WRITE (lu,'(/,10x,A,I12)' ) 'Number of Points = ', npts

CONTAINS
   FUNCTION real_to_char ( value )
   !
   !  Purpose:
   !    To convert a real value into a 12-character string, with the
   !    number printed in as readable a format as possible considering
   !    its range.  This routine prints out the number according to the
   !    following rules:
   !       1.  value > 9999999.                ES12.5
   !       2.  value < -999999.                ES12.5
   !       3.  0.    <  ABS(value) < 0.01      ES12.5
   !       5.  value = 0.0                     F12.4
   !       6.  Otherwise                       F12.4
   !
   !  Record of revisions:
   !      Date       Programmer          Description of change
   !      ====       ==========          =====================
   !    11/26/95    S. J. Chapman        Original code
   !
   IMPLICIT NONE

   ! Declare calling arguments:
   REAL(KIND=kind), INTENT(IN) :: value ! value to convert to char form
   CHARACTER (len=12) :: real_to_char  ! Output character string

   ! Declare local variables:
   CHARACTER(len=9) :: fmt             ! Format descriptor
   CHARACTER(len=12) :: string = ' '   ! Output string

   ! Select proper format
   IF ( value > 9999999. ) THEN
      fmt = '(ES12.5)'
   ELSE IF ( value < -999999. ) THEN
      fmt = '(ES12.5)'
   ELSE IF ( value == 0. ) THEN
      fmt = '(F12.4)'
   ELSE IF ( ABS(value) < 0.01 ) THEN
      fmt = '(ES12.5)'
   ELSE
      fmt = '(F12.4)'
   END IF

   ! Convert to character form.
   WRITE (string,fmt) value
   real_to_char = string

   END FUNCTION real_to_char

END SUBROUTINE plotxy_sgl

SUBROUTINE plotxy_dbl ( x, y, npts, lu, minx, maxx, miny, maxy, &
                        nbinx, nbiny )
!
!  Purpose:
!    Subroutine to plot (x,y) pairs of points.  Note that
!    x and y scales may either be set or defaulted on an
!    item-by-item basis.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    03/20/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=12) ! Precision
INTEGER, PARAMETER :: maxbin = 65           ! Max num. of bins

! List of dummy arguments:
INTEGER, INTENT(IN) :: npts                 ! No. points to plot
REAL(KIND=kind),INTENT(IN),DIMENSION(npts) :: x  ! x data
REAL(KIND=kind),INTENT(IN),DIMENSION(npts) :: y  ! y data
INTEGER, INTENT(IN) :: lu                   ! Output unit
REAL(KIND=kind),INTENT(IN),OPTIONAL :: minx ! Smallest x to plot
REAL(KIND=kind),INTENT(IN),OPTIONAL :: maxx ! Largest x to plot
REAL(KIND=kind),INTENT(IN),OPTIONAL :: miny ! Smallest y to plot
REAL(KIND=kind),INTENT(IN),OPTIONAL :: maxy ! Largest y to plot
INTEGER,INTENT(IN),OPTIONAL :: nbinx        ! No. bins for x
INTEGER,INTENT(IN),OPTIONAL :: nbiny        ! No. bins for y

! List of local variables:
CHARACTER(len=14) :: annot    ! Annotation (y-amplitude)
REAL(KIND=kind) :: frac       ! Fraction of plot width
INTEGER :: i                  ! Loop index
!INTEGER :: ibinx0            ! Bin for x zero crossing
INTEGER :: ibiny              ! Bin for current y value
INTEGER :: ibiny0             ! Bin for y zero crossing
INTEGER :: j                  ! Loop index
REAL(KIND=kind) :: maxx1      ! Largest x to plot
REAL(KIND=kind) :: minx1      ! Samllest x to plot
REAL(KIND=kind) :: maxy1      ! Largest y to plot
REAL(KIND=kind) :: miny1      ! Samllest y to plot
INTEGER :: nbinx1             ! No. of x bins to plot
INTEGER :: nbiny1             ! No. of y bins to plot
CHARACTER(len=65) :: pltbuf   ! Plotting buffer
CHARACTER(len=65) :: scl      ! Scale on plot border
CHARACTER(len=65),SAVE:: scl1 ! Template for scale
CHARACTER(len=80) :: sclbuf   ! Plotting scale buffer
REAL(KIND=kind) :: xpm        ! Beginning of x range (cur x pos)
REAL(KIND=kind) :: xpos       ! Center of x range (cur x pos)
REAL(KIND=kind) :: xpp        ! Ending of x range (cur x pos)

! Initialize scale border.
scl1 = '+-----------------------------&
       &-----------------------------------'

! If any scale values are defaulted, set min and max of x and y
! axes according to the input data.  Note that we are making
! local copies of minx1, maxx1, miny1, and maxy1 if the exist,
! so that we may modify their values.  Check for and set defaults
! on all optional arguments.
!
IF ( PRESENT(minx) ) THEN
   minx1 = minx
ELSE
   minx1 = minval(x)
END IF

IF ( PRESENT(maxx) ) THEN
   maxx1 = maxx
ELSE
   maxx1 = maxval(x)
END IF

IF ( PRESENT(miny) ) THEN
   miny1 = miny
ELSE
   miny1 = minval(y)
END IF

IF ( PRESENT(maxy) ) THEN
   maxy1 = maxy
ELSE
   maxy1 = maxval(y)
END IF

IF ( PRESENT(nbinx) ) THEN
   nbinx1 = MAX(1, MIN(nbinx,maxbin) )
ELSE
   nbinx1 = maxbin
END IF

IF ( PRESENT(nbiny) ) THEN
   nbiny1 = MAX(1, MIN(nbiny,maxbin) )
ELSE
   nbiny1 = maxbin
END IF

! We will divide minx1 to maxx1 into nbinx1 bins for
! plotting purposes.  Locate the zero bin, if it is
! between minx1 and maxx1.
!IF ( (maxx1 > 0.) .AND. (minx1 < 0.) ) THEN
!   frac = ( 0. - minx1) / (maxx1 - minx1 )
!   ibinx0 = NINT ( (nbinx1-1) * frac ) + 1
!ELSE
!   ibinx0 = 0
!END IF

! We will divide miny1 to maxy1 into nbiny1 bins for
! plotting purposes.  Locate the zero bin, if it is
! between miny1 and maxy1.
 IF ( (maxy1 > 0.) .AND. (miny1 < 0.) ) THEN
   frac = ( 0. - miny1) / (maxy1 - miny1 )
   ibiny0 = NINT ( (nbiny1-1) * frac ) + 1
ELSE
   ibiny0 = 0
END IF

! Set upper right hand corner of plot.
annot = ' '
scl   = scl1(1:nbiny1)
scl(nbiny1:nbiny1) = '+'

! Set y axis zero into border of plot, if it is within
! the limits of the plot.
IF ( ibiny0 > 0 ) THEN
   scl(ibiny0:ibiny0) = '+'
END IF

! Generate plot scale label.
sclbuf = ' '
sclbuf(10:22) = real_to_char(miny1)
sclbuf(15+nbiny1-11:15+nbiny1) = real_to_char(maxy1)

! Set y axis zero into scale of plot, if it is within
! the limits of the plot.
IF ( ibiny0 > 0 ) THEN
   sclbuf(15+ibiny0:15+ibiny0) = '0'
END IF

! Print upper border.
WRITE (lu,'(A)') sclbuf
WRITE (lu,'(A,1x,A)') annot, scl

! Plot data points.
plot: DO i = 1, nbinx1

   ! Calculate x position of this data point.
   xpos = ((REAL(I)-1.0)/(nbinx1-1))*(maxx1-minx1) + minx1
   xpm  = ((REAL(I)-1.5)/(nbinx1-1))*(maxx1-minx1) + minx1
   xpp  = ((REAL(I)-0.5)/(nbinx1-1))*(maxx1-minx1) + minx1

   ! Clear line.
   pltbuf = ' '
   annot  = ' '

   ! Set value of x data point.
   annot(2:13) = real_to_char( xpos )

   ! Set min and max borders.
   pltbuf(1:1)   = '|'
   pltbuf(nbiny1:nbiny1) = '|'

   ! Set zero line, if within borders.
   IF ( ibiny0 > 0 ) THEN
      pltbuf(ibiny0:ibiny0) = '|'
   END IF

   ! Find points that map to this x position.
   xpts: DO j = 1, npts
      IF ( (x(j) > xpm) .AND. (x(j) < xpp) ) THEN

         ! Plot point on array.
         frac  = ( y(J) - miny1) / (maxy1 - miny1 )
         ibiny = NINT ( (nbiny1-1) * frac ) + 1
         IF ( (ibiny >= 1) .AND. (ibiny <= nbiny1) ) THEN
            pltbuf(ibiny:ibiny) = '*'
         END IF
      END IF
   END DO xpts

   ! Output line.
   WRITE (lu,'(A,1X,A)') annot, pltbuf

END DO plot

! Print lower border.
annot = ' '
WRITE (lu,'(A,1X,A)') annot, scl
WRITE (lu,'(A)')      sclbuf

! Print out summary info.
WRITE (lu,'(/,10x,A,I12)' ) 'Number of Points = ', npts

CONTAINS
   FUNCTION real_to_char ( value )
   !
   !  Purpose:
   !    To convert a real value into a 12-character string, with the
   !    number printed in as readable a format as possible considering
   !    its range.  This routine prints out the number according to the
   !    following rules:
   !       1.  value > 9999999.                ES12.5
   !       2.  value < -999999.                ES12.5
   !       3.  0.    <  ABS(value) < 0.01      ES12.5
   !       5.  value = 0.0                     F12.4
   !       6.  Otherwise                       F12.4
   !
   !  Record of revisions:
   !      Date       Programmer          Description of change
   !      ====       ==========          =====================
   !    11/26/95    S. J. Chapman        Original code
   !
   IMPLICIT NONE

   ! Declare calling arguments:
   REAL(KIND=kind), INTENT(IN) :: value ! value to convert to char form
   CHARACTER (len=12) :: real_to_char  ! Output character string

   ! Declare local variables:
   CHARACTER(len=9) :: fmt             ! Format descriptor
   CHARACTER(len=12) :: string = ' '   ! Output string

   ! Select proper format
   IF ( value > 9999999. ) THEN
      fmt = '(ES12.5)'
   ELSE IF ( value < -999999. ) THEN
      fmt = '(ES12.5)'
   ELSE IF ( value == 0. ) THEN
      fmt = '(F12.4)'
   ELSE IF ( ABS(value) < 0.01 ) THEN
      fmt = '(ES12.5)'
   ELSE
      fmt = '(F12.4)'
   END IF

   ! Convert to character form.
   WRITE (string,fmt) value
   real_to_char = string

   END FUNCTION real_to_char

END SUBROUTINE plotxy_dbl

SUBROUTINE simul_sgl ( a, b, soln, ndim, n, error )
!
!  Purpose:
!    Subroutine to solve a set of N linear equations in N
!    unknowns using Gaussian elimination and the maximum
!    pivot technique.  This version of simul has been
!    modified to use array sections and automatic arrays.
!    It DOES NOT DESTROY the original input values.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    10/16/95    S. J. Chapman        Original code
! 1. 10/17/95    S. J. Chapman        Modified to use array ops
! 2. 03/02/96    S. J. Chapman        Add allocatable arrays
!
IMPLICIT NONE

! Declare calling arguments:
INTEGER, INTENT(IN) :: ndim          ! Dimension of arrays a and b
REAL, INTENT(IN), DIMENSION(ndim,ndim) :: a
                                     ! Array of coefficients (N x N).
                                     ! This array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
REAL, INTENT(IN), DIMENSION(ndim) :: b
                                     ! Input: Right-hand side of eqns.
REAL, INTENT(OUT), DIMENSION(ndim) :: soln
                                     ! Output: Solution vector.
INTEGER, INTENT(IN) :: n             ! Number of equations to solve.
INTEGER, INTENT(OUT) :: error        ! Error flag:
                                     !   0 -- No error
                                     !   1 -- Singular equations

! Declare local parameters
REAL, PARAMETER :: epsilon = 1.0E-6  ! A "small" number for comparison
                                     ! when determining singular eqns

! Declare local variables:
REAL, DIMENSION(n,n) :: a1           ! Copy of "a" which will be
                                     ! destroyed during the solution
REAL :: factor                       ! Factor to multiply eqn irow by
                                     ! before adding to eqn jrow
INTEGER :: irow                      ! Number of the equation currently
                                     ! currently being processed
INTEGER :: ipeak                     ! Pointer to equation containing
                                     ! maximum pivot value
INTEGER :: jrow                      ! Number of the equation compared
                                     ! to the current equation
REAL :: temp                         ! Scratch value
REAL, DIMENSION(n) :: temp1          ! Scratch array

! Make copies of arrays "a" and "b" for local use
a1 = a(1:n,1:n)
soln = b(1:n)

! Process N times to get all equations...
mainloop: DO irow = 1, n

   ! Find peak pivot for column irow in rows irow to N
   ipeak = irow
   max_pivot: DO jrow = irow+1, n
      IF (ABS(a1(jrow,irow)) > ABS(a1(ipeak,irow))) THEN
         ipeak = jrow
      END IF
   END DO max_pivot

   ! Check for singular equations.
   singular: IF ( ABS(a1(ipeak,irow)) < epsilon ) THEN
      error = 1
      RETURN
   END IF singular

   ! Otherwise, if ipeak /= irow, swap equations irow & ipeak
   swap_eqn: IF ( ipeak /= irow ) THEN
      temp1 = a1(ipeak,1:n)
      a1(ipeak,1:n) = a1(irow,1:n)   ! Swap rows in a
      a1(irow,1:n) = temp1
      temp = soln(ipeak)
      soln(ipeak) = soln(irow)       ! Swap rows in b
      soln(irow)  = temp
   END IF swap_eqn

   ! Multiply equation irow by -a1(jrow,irow)/a1(irow,irow),
   ! and add it to Eqn jrow (for all eqns except irow itself).
   eliminate: DO jrow = 1, n
      IF ( jrow /= irow ) THEN
         factor = -a1(jrow,irow)/a1(irow,irow)
         a1(jrow,1:n) = a1(irow,1:n)*factor + a1(jrow,1:n)
         soln(jrow) = soln(irow)*factor + soln(jrow)
      END IF
   END DO eliminate
END DO mainloop

! End of main loop over all equations.  All off-diagonal
! terms are now zero.  To get the final answer, we must
! divide each equation by the coefficient of its on-diagonal
! term.
divide: DO irow = 1, n
   soln(irow) = soln(irow) / a1(irow,irow)
END DO divide

! Set error flag to 0 and return.
error = 0

END SUBROUTINE simul_sgl

SUBROUTINE simul_dbl ( a, b, soln, ndim, n, error )
!
!  Purpose:
!    Subroutine to solve a set of N linear equations in N
!    unknowns using Gaussian elimination and the maximum
!    pivot technique.  This version of simul has been
!    modified to use array sections and automatic arrays.
!    It uses double precision arithmetic to avoid
!    cumulative roundoff errors.  It DOES NOT DESTROY the
!    original input values.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    10/16/95    S. J. Chapman        Original code
! 1. 10/17/95    S. J. Chapman        Modified to use array ops
! 2. 03/02/96    S. J. Chapman        Add allocatable arrays
! 3. 05/08/96    S. J. Chapman        Double precision
!
IMPLICIT NONE

! Declare parameters
INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND(p=12)

! Declare calling arguments:
INTEGER, INTENT(IN) :: ndim          ! Dimension of arrays a and b
REAL(KIND=dbl), INTENT(IN), DIMENSION(ndim,ndim) :: a
                                     ! Array of coefficients (N x N).
                                     ! This array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
REAL(KIND=dbl), INTENT(IN), DIMENSION(ndim) :: b
                                     ! Input: Right-hand side of eqns.
REAL(KIND=dbl), INTENT(OUT), DIMENSION(ndim) :: soln
                                     ! Output: Solution vector.
INTEGER, INTENT(IN) :: n             ! Number of equations to solve.
INTEGER, INTENT(OUT) :: error        ! Error flag:
                                     !   0 -- No error
                                     !   1 -- Singular equations

! Declare local parameters
REAL(KIND=dbl), PARAMETER :: epsilon = 1.0E-12
                                     ! A "small" number for comparison
                                     ! when determining singular eqns

! Declare local variables:
REAL(KIND=dbl), DIMENSION(n,n) :: a1 ! Copy of "a" which will be
                                     ! destroyed during the solution
REAL(KIND=dbl) :: factor             ! Factor to multiply eqn irow by
                                     ! before adding to eqn jrow
INTEGER :: irow                      ! Number of the equation currently
                                     ! currently being processed
INTEGER :: ipeak                     ! Pointer to equation containing
                                     ! maximum pivot value
INTEGER :: jrow                      ! Number of the equation compared
                                     ! to the current equation
REAL(KIND=dbl) :: temp               ! Scratch value
REAL(KIND=dbl),DIMENSION(n) :: temp1 ! Scratch array

! Make copies of arrays "a" and "b" for local use
a1 = a(1:n,1:n)
soln = b(1:n)

! Process N times to get all equations...
mainloop: DO irow = 1, n

   ! Find peak pivot for column irow in rows irow to N
   ipeak = irow
   max_pivot: DO jrow = irow+1, n
      IF (ABS(a1(jrow,irow)) > ABS(a1(ipeak,irow))) THEN
         ipeak = jrow
      END IF
   END DO max_pivot

   ! Check for singular equations.
   singular: IF ( ABS(a1(ipeak,irow)) < epsilon ) THEN
      error = 1
      RETURN
   END IF singular

   ! Otherwise, if ipeak /= irow, swap equations irow & ipeak
   swap_eqn: IF ( ipeak /= irow ) THEN
      temp1 = a1(ipeak,1:n)
      a1(ipeak,1:n) = a1(irow,1:n)   ! Swap rows in a
      a1(irow,1:n) = temp1
      temp = soln(ipeak)
      soln(ipeak) = soln(irow)       ! Swap rows in b
      soln(irow)  = temp
   END IF swap_eqn

   ! Multiply equation irow by -a1(jrow,irow)/a1(irow,irow),
   ! and add it to Eqn jrow (for all eqns except irow itself).
   eliminate: DO jrow = 1, n
      IF ( jrow /= irow ) THEN
         factor = -a1(jrow,irow)/a1(irow,irow)
         a1(jrow,1:n) = a1(irow,1:n)*factor + a1(jrow,1:n)
         soln(jrow) = soln(irow)*factor + soln(jrow)
      END IF
   END DO eliminate
END DO mainloop

! End of main loop over all equations.  All off-diagonal
! terms are now zero.  To get the final answer, we must
! divide each equation by the coefficient of its on-diagonal
! term.
divide: DO irow = 1, n
   soln(irow) = soln(irow) / a1(irow,irow)
END DO divide

! Set error flag to 0 and return.
error = 0

END SUBROUTINE simul_dbl

SUBROUTINE simul_sgl_cmplx ( a, b, soln, ndim, n, error )
!
!  Purpose:
!    Subroutine to solve a set of N linear equations in N
!    unknowns using Gaussian elimination and the maximum
!    pivot technique.  This version of simul has been
!    modified to use array sections and automatic arrays.
!    It uses single precision complex arithmetic.  It DOES
!    NOT DESTROY the original input values.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    10/16/95    S. J. Chapman        Original code
! 1. 10/17/95    S. J. Chapman        Modified to use array ops
! 2. 03/02/96    S. J. Chapman        Add allocatable arrays
! 3. 05/08/96    S. J. Chapman        Double precision
! 4. 05/09/96    S. J. Chapman        Complex version
!
IMPLICIT NONE

! Declare parameters
INTEGER, PARAMETER :: sgl = SELECTED_REAL_KIND(p=6)

! Declare calling arguments:
INTEGER, INTENT(IN) :: ndim          ! Dimension of arrays a and b
COMPLEX(KIND=sgl), INTENT(IN), DIMENSION(ndim,ndim) :: a
                                     ! Array of coefficients (N x N).
                                     ! This array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
COMPLEX(KIND=sgl), INTENT(IN), DIMENSION(ndim) :: b
                                     ! Input: Right-hand side of eqns.
COMPLEX(KIND=sgl), INTENT(OUT), DIMENSION(ndim) :: soln
                                     ! Output: Solution vector.
INTEGER, INTENT(IN) :: n             ! Number of equations to solve.
INTEGER, INTENT(OUT) :: error        ! Error flag:
                                     !   0 -- No error
                                     !   1 -- Singular equations

! Declare local parameters
REAL(KIND=sgl), PARAMETER :: epsilon = 1.0E-6
                                     ! A "small" number for comparison
                                     ! when determining singular eqns

! Declare local variables:
REAL(KIND=sgl), DIMENSION(n,n) :: a1 ! Copy of "a" which will be
                                     ! destroyed during the solution
REAL(KIND=sgl) :: factor             ! Factor to multiply eqn irow by
                                     ! before adding to eqn jrow
INTEGER :: irow                      ! Number of the equation currently
                                     ! currently being processed
INTEGER :: ipeak                     ! Pointer to equation containing
                                     ! maximum pivot value
INTEGER :: jrow                      ! Number of the equation compared
                                     ! to the current equation
COMPLEX(KIND=sgl) :: temp            ! Scratch value
COMPLEX(KIND=sgl),DIMENSION(n) :: temp1 ! Scratch array

! Make copies of arrays "a" and "b" for local use
a1 = a(1:n,1:n)
soln = b(1:n)

! Process N times to get all equations...
mainloop: DO irow = 1, n

   ! Find peak pivot for column irow in rows irow to N
   ipeak = irow
   max_pivot: DO jrow = irow+1, n
      IF (ABS(a1(jrow,irow)) > ABS(a1(ipeak,irow))) THEN
         ipeak = jrow
      END IF
   END DO max_pivot

   ! Check for singular equations.
   singular: IF ( ABS(a1(ipeak,irow)) < epsilon ) THEN
      error = 1
      RETURN
   END IF singular

   ! Otherwise, if ipeak /= irow, swap equations irow & ipeak
   swap_eqn: IF ( ipeak /= irow ) THEN
      temp1 = a1(ipeak,1:n)
      a1(ipeak,1:n) = a1(irow,1:n)   ! Swap rows in a
      a1(irow,1:n) = temp1
      temp = soln(ipeak)
      soln(ipeak) = soln(irow)       ! Swap rows in b
      soln(irow)  = temp
   END IF swap_eqn

   ! Multiply equation irow by -a1(jrow,irow)/a1(irow,irow),
   ! and add it to Eqn jrow (for all eqns except irow itself).
   eliminate: DO jrow = 1, n
      IF ( jrow /= irow ) THEN
         factor = -a1(jrow,irow)/a1(irow,irow)
         a1(jrow,1:n) = a1(irow,1:n)*factor + a1(jrow,1:n)
         soln(jrow) = soln(irow)*factor + soln(jrow)
      END IF
   END DO eliminate
END DO mainloop

! End of main loop over all equations.  All off-diagonal
! terms are now zero.  To get the final answer, we must
! divide each equation by the coefficient of its on-diagonal
! term.
divide: DO irow = 1, n
   soln(irow) = soln(irow) / a1(irow,irow)
END DO divide

! Set error flag to 0 and return.
error = 0

END SUBROUTINE simul_sgl_cmplx

SUBROUTINE simul_dbl_cmplx ( a, b, soln, ndim, n, error )
!
!  Purpose:
!    Subroutine to solve a set of N linear equations in N
!    unknowns using Gaussian elimination and the maximum
!    pivot technique.  This version of simul has been
!    modified to use array sections and automatic arrays.
!    It uses double precision complex arithmetic.  It DOES
!    NOT DESTROY the original input values.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    10/16/95    S. J. Chapman        Original code
! 1. 10/17/95    S. J. Chapman        Modified to use array ops
! 2. 03/02/96    S. J. Chapman        Add allocatable arrays
! 3. 05/08/96    S. J. Chapman        Double precision
! 4. 05/09/96    S. J. Chapman        Double complex version
!
IMPLICIT NONE

! Declare parameters
INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND(p=12)

! Declare calling arguments:
INTEGER, INTENT(IN) :: ndim          ! Dimension of arrays a and b
COMPLEX(KIND=dbl), INTENT(IN), DIMENSION(ndim,ndim) :: a
                                     ! Array of coefficients (N x N).
                                     ! This array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
COMPLEX(KIND=dbl), INTENT(IN), DIMENSION(ndim) :: b
                                     ! Input: Right-hand side of eqns.
COMPLEX(KIND=dbl), INTENT(OUT), DIMENSION(ndim) :: soln
                                     ! Output: Solution vector.
INTEGER, INTENT(IN) :: n             ! Number of equations to solve.
INTEGER, INTENT(OUT) :: error        ! Error flag:
                                     !   0 -- No error
                                     !   1 -- Singular equations

! Declare local parameters
REAL(KIND=dbl), PARAMETER :: epsilon = 1.0E-12
                                     ! A "small" number for comparison
                                     ! when determining singular eqns

! Declare local variables:
REAL(KIND=dbl), DIMENSION(n,n) :: a1 ! Copy of "a" which will be
                                     ! destroyed during the solution
REAL(KIND=dbl) :: factor             ! Factor to multiply eqn irow by
                                     ! before adding to eqn jrow
INTEGER :: irow                      ! Number of the equation currently
                                     ! currently being processed
INTEGER :: ipeak                     ! Pointer to equation containing
                                     ! maximum pivot value
INTEGER :: jrow                      ! Number of the equation compared
                                     ! to the current equation
COMPLEX(KIND=dbl) :: temp            ! Scratch value
COMPLEX(KIND=dbl),DIMENSION(n) :: temp1 ! Scratch array

! Make copies of arrays "a" and "b" for local use
a1 = a(1:n,1:n)
soln = b(1:n)

! Process N times to get all equations...
mainloop: DO irow = 1, n

   ! Find peak pivot for column irow in rows irow to N
   ipeak = irow
   max_pivot: DO jrow = irow+1, n
      IF (ABS(a1(jrow,irow)) > ABS(a1(ipeak,irow))) THEN
         ipeak = jrow
      END IF
   END DO max_pivot

   ! Check for singular equations.
   singular: IF ( ABS(a1(ipeak,irow)) < epsilon ) THEN
      error = 1
      RETURN
   END IF singular

   ! Otherwise, if ipeak /= irow, swap equations irow & ipeak
   swap_eqn: IF ( ipeak /= irow ) THEN
      temp1 = a1(ipeak,1:n)
      a1(ipeak,1:n) = a1(irow,1:n)   ! Swap rows in a
      a1(irow,1:n) = temp1
      temp = soln(ipeak)
      soln(ipeak) = soln(irow)       ! Swap rows in b
      soln(irow)  = temp
   END IF swap_eqn

   ! Multiply equation irow by -a1(jrow,irow)/a1(irow,irow),
   ! and add it to Eqn jrow (for all eqns except irow itself).
   eliminate: DO jrow = 1, n
      IF ( jrow /= irow ) THEN
         factor = -a1(jrow,irow)/a1(irow,irow)
         a1(jrow,1:n) = a1(irow,1:n)*factor + a1(jrow,1:n)
         soln(jrow) = soln(irow)*factor + soln(jrow)
      END IF
   END DO eliminate
END DO mainloop

! End of main loop over all equations.  All off-diagonal
! terms are now zero.  To get the final answer, we must
! divide each equation by the coefficient of its on-diagonal
! term.
divide: DO irow = 1, n
   soln(irow) = soln(irow) / a1(irow,irow)
END DO divide

! Set error flag to 0 and return.
error = 0

END SUBROUTINE simul_dbl_cmplx

SUBROUTINE ssort_sgl(array, n)
!
!  Purpose:
!    To sort real array "array" into ascending order using a selection
!    sort.
!
IMPLICIT NONE

! Declare parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=6)

! Declare calling parameters:
INTEGER, INTENT(IN) :: n                             ! No of vals
REAL(KIND=kind),DIMENSION(n),INTENT(INOUT) :: array  ! Array to sort

! Declare local variables:
INTEGER :: i                  ! Loop index
INTEGER :: iptr               ! Pointer to smallest value
INTEGER :: j                  ! Loop index
REAL(KIND=kind) :: temp       ! Temp variable for swaps

! Sort the array
outer: DO i = 1, n-1

   ! Find the minimum value in array(i) through array(n)
   iptr = i
   inner: DO j = i+1, n
      minval: IF ( array(j) < array(iptr) ) THEN
         iptr = j
      END IF minval
   END DO inner

   ! iptr now points to the minimum value, so swap array(iptr)
   ! with array(i) if i /= iptr.
   swap: IF ( i /= iptr ) THEN
      temp        = array(i)
      array(i)    = array(iptr)
      array(iptr) = temp
   END IF swap

END DO outer

END SUBROUTINE ssort_sgl

SUBROUTINE ssort_dbl(array, n)
!
!  Purpose:
!    To sort real array "array" into ascending order using a selection
!    sort.
!
IMPLICIT NONE

! Declare parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=12)

! Declare calling parameters:
INTEGER, INTENT(IN) :: n                             ! No of vals
REAL(KIND=kind),DIMENSION(n),INTENT(INOUT) :: array  ! Array to sort

! Declare local variables:
INTEGER :: i                  ! Loop index
INTEGER :: iptr               ! Pointer to smallest value
INTEGER :: j                  ! Loop index
REAL(KIND=kind) :: temp       ! Temp variable for swaps

! Sort the array
outer: DO i = 1, n-1

   ! Find the minimum value in array(i) through array(n)
   iptr = i
   inner: DO j = i+1, n
      minval: IF ( array(j) < array(iptr) ) THEN
         iptr = j
      END IF minval
   END DO inner

   ! iptr now points to the minimum value, so swap array(iptr)
   ! with array(i) if i /= iptr.
   swap: IF ( i /= iptr ) THEN
      temp        = array(i)
      array(i)    = array(iptr)
      array(iptr) = temp
   END IF swap

END DO outer

END SUBROUTINE ssort_dbl

SUBROUTINE ssort_int(array, n)
!
!  Purpose:
!    To sort real array "array" into ascending order using a selection
!    sort.
!
IMPLICIT NONE

! Declare calling parameters:
INTEGER, INTENT(IN) :: n                     ! No of vals
INTEGER,DIMENSION(n),INTENT(INOUT) :: array  ! Array to sort

! Declare local variables:
INTEGER :: i                  ! Loop index
INTEGER :: iptr               ! Pointer to smallest value
INTEGER :: j                  ! Loop index
INTEGER :: temp               ! Temp variable for swaps

! Sort the array
outer: DO i = 1, n-1

   ! Find the minimum value in array(i) through array(n)
   iptr = i
   inner: DO j = i+1, n
      minval: IF ( array(j) < array(iptr) ) THEN
         iptr = j
      END IF minval
   END DO inner

   ! iptr now points to the minimum value, so swap array(iptr)
   ! with array(i) if i /= iptr.
   swap: IF ( i /= iptr ) THEN
      temp        = array(i)
      array(i)    = array(iptr)
      array(iptr) = temp
   END IF swap

END DO outer

END SUBROUTINE ssort_int

SUBROUTINE ssort_char(array, n)
!
!  Purpose:
!    To sort real array "array" into ascending order using a selection
!    sort.
!
IMPLICIT NONE

! Declare calling parameters:
INTEGER, INTENT(IN) :: n                              ! No of vals
CHARACTER(len=*),DIMENSION(n),INTENT(INOUT) :: array  ! Array to sort

! Declare local variables:
INTEGER :: i                  ! Loop index
INTEGER :: iptr               ! Pointer to smallest value
INTEGER :: j                  ! Loop index
CHARACTER(len=LEN(array)) :: temp  ! Temp variable for swaps

! Sort the array
outer: DO i = 1, n-1

   ! Find the minimum value in array(i) through array(n)
   iptr = i
   inner: DO j = i+1, n
      minval: IF ( array(j) < array(iptr) ) THEN
         iptr = j
      END IF minval
   END DO inner

   ! iptr now points to the minimum value, so swap array(iptr)
   ! with array(i) if i /= iptr.
   swap: IF ( i /= iptr ) THEN
      temp        = array(i)
      array(i)    = array(iptr)
      array(iptr) = temp
   END IF swap

END DO outer

END SUBROUTINE ssort_char

SUBROUTINE deriv ( f, x0, dx, dfdx, error )
!
!  Purpose:
!    To take the derivative of function f(x) at point x0
!    using step size dx.  If dx = 0., then take derivative
!    with as much accuracy as possible.  This subroutine
!    expects the function f(x)) to be passed as a calling
!    argument.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    05/02/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=6) ! Precision
INTEGER, PARAMETER :: nsteps = 10   ! Maximum number of steps
                                    ! to achieve accuracy

! Declare calling arguments
REAL(KIND=kind), EXTERNAL :: f      ! Function to differentiate
REAL(KIND=kind),INTENT(IN) :: x0    ! Point to get derivative
REAL(KIND=kind),INTENT(INOUT) :: dx ! Desired step size.  If dx = 0,
                                    ! actual step size is returned.
REAL(KIND=kind),INTENT(OUT) :: dfdx ! Derivative
INTEGER,INTENT(OUT) :: error        ! Error flag: 0=no error
                                    !             1=dx<0

! List of local variables:
REAL(KIND=kind) :: delta            ! ABS(dfdx - dfdx0)
REAL(KIND=kind) :: delta0           ! Last delta
REAL(KIND=kind) :: dfdx0            ! Derivative at last step size
INTEGER :: i                        ! Index variable

! If dx < 0., this is an error.
IF ( dx < 0. ) THEN
   error = 1
   return

! IF dx is specified, then calculate derivative using the
! central difference method and the specified dx.
ELSE IF ( DX > 0. ) THEN
   dfdx = (f(x0+dx/2.) - f(x0-dx/2.) ) / dx
   error = 0
   RETURN

! IF dx is not specified, then calculate derivative using the
! central difference method, starting with dx = 0.1 and going
! down until the maximum accuracy is achieved.
ELSE
   error = 0
   dx    = 0.1
   dfdx0 = (f(x0+dx/2.) - f(x0-dx/2.) ) / dx
   delta0 = dfdx0

!  Iterate for most precision.  Note that we do a maximum
!  of NSTEPS-1 iterations to avoid the possibility of
!  infinite loops.
   DO i = 2, nsteps
      dx    = dx / 10.
      dfdx  = (f(x0+dx/2.) - f(x0-dx/2.) ) / dx
      delta = ABS ( dfdx - dfdx0 )
      IF ( delta > delta0 ) THEN    ! We have passed optimal point
         dfdx = dfdx0
         dx   = dx * 10.
         RETURN
      ELSE                          ! Keep looking
         delta0 = delta
         dfdx0  = dfdx
      END IF
   END DO
END IF

END SUBROUTINE deriv

SUBROUTINE dderiv ( f, x0, dx, dfdx, error )
!
!  Purpose:
!    To take the derivative of function f(x) at point x0
!    using step size dx.  If dx = 0., then take derivative
!    with as much accuracy as possible.  This subroutine
!    expects the function f(x)) to be passed as a calling
!    argument.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    05/02/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=12) ! Precision
INTEGER, PARAMETER :: nsteps = 20   ! Maximum number of steps
                                    ! to achieve accuracy

! Declare calling arguments
REAL(KIND=kind), EXTERNAL :: f      ! Function to differentiate
REAL(KIND=KIND),INTENT(IN) :: x0    ! Point to get derivative
REAL(KIND=kind),INTENT(INOUT) :: dx ! Desired step size.  If dx = 0,
                                    ! actual step size is returned.
REAL(KIND=kind),INTENT(OUT) :: dfdx ! Derivative
INTEGER,INTENT(OUT) :: error        ! Error flag: 0=no error
                                    !             1=dx<0

! List of local variables:
REAL(KIND=kind) :: delta            ! ABS(dfdx - dfdx0)
REAL(KIND=kind) :: delta0           ! Last delta
REAL(KIND=kind) :: dfdx0            ! Derivative at last step size
INTEGER :: i                        ! Index variable

! If dx < 0., this is an error.
IF ( dx < 0. ) THEN
   error = 1
   return

! IF dx is specified, then calculate derivative using the
! central difference method and the specified dx.
ELSE IF ( dx > 0. ) THEN
   dfdx = (f(x0+dx/2.) - f(x0-dx/2.) ) / dx
   error = 0
   RETURN

! IF dx is not specified, then calculate derivative using the
! central difference method, starting with dx = 0.1 and going
! down until the maximum accuracy is achieved.
ELSE
   error = 0
   dx    = 0.1
   dfdx0 = (f(x0+dx/2.) - f(x0-dx/2.) ) / dx
   delta0 = dfdx0

!  Iterate for most precision.  Note that we do a maximum
!  of NSTEPS-1 iterations to avoid the possibility of
!  infinite loops.
   DO i = 2, nsteps
      dx    = dx / 10.
      dfdx  = (f(x0+dx/2.) - f(x0-dx/2.) ) / dx
      delta = ABS ( dfdx - dfdx0 )
      IF ( delta > delta0 ) THEN    ! We have passed optimal point
         dfdx = dfdx0
         dx   = dx * 10.
         RETURN
      ELSE                          ! Keep looking
         delta0 = delta
         dfdx0  = dfdx
      END IF
   END DO
END IF

END SUBROUTINE dderiv

SUBROUTINE dft_sgl ( array_in, array_out, n )
!
!  Purpose:
!    An implementation of the Discrete Fourier Transform from
!    the defining equation.
!
!  Record of revisions:
!       Date       Programmer          Description of change
!       ====       ==========          =====================
!     04/05/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=6) ! Precision
REAL(KIND=kind), PARAMETER :: twopi = 6.283185307179586_kind ! 2*pi

! Declare calling arguments:
INTEGER, INTENT(IN) :: n             ! No. of elements in DFT
COMPLEX(KIND=kind),INTENT(IN),DIMENSION(0:n-1) :: array_in
                                     ! Input data for DFT
COMPLEX(KIND=kind),INTENT(INOUT),DIMENSION(0:n-1) :: array_out
                                     ! Output data for DFT

! List of local variables:
INTEGER :: i, k                      ! Index variables

! Calculate the DFT from its definition.
DO i = 0, n-1
   array_out(i) = (0.,0.)
   DO k = 0, n-1
      array_out(i) = array_out(i) + array_in(k) &
                   * EXP( -twopi*(0._kind,1._kind)*k*i/n )
   END DO
END DO

END SUBROUTINE dft_sgl

SUBROUTINE dft_dbl ( array_in, array_out, n )
!
!  Purpose:
!    An implementation of the Discrete Fourier Transform from
!    the defining equation.
!
!  Record of revisions:
!       Date       Programmer          Description of change
!       ====       ==========          =====================
!     04/05/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=12) ! Precision
REAL(KIND=kind), PARAMETER :: twopi = 6.283185307179586_kind ! 2*pi

! Declare calling arguments:
INTEGER, INTENT(IN) :: n             ! No. of elements in DFT
COMPLEX(KIND=kind),INTENT(IN),DIMENSION(0:n-1) :: array_in
                                     ! Input data for DFT
COMPLEX(KIND=kind),INTENT(INOUT),DIMENSION(0:n-1) :: array_out
                                     ! Output data for DFT

! List of local variables:
INTEGER :: i, k                      ! Index variables

! Calculate the DFT from its definition.
DO i = 0, n-1
   array_out(i) = (0.,0.)
   DO k = 0, n-1
      array_out(i) = array_out(i) + array_in(k) &
                   * EXP( -twopi*(0._kind,1._kind)*k*i/n )
   END DO
END DO

END SUBROUTINE dft_dbl

SUBROUTINE idft_sgl ( array_in, array_out, n )
!
!  Purpose:
!    An implementation of the inverse Discrete Fourier Transform
!    from the defining equation.
!
!  Record of revisions:
!       Date       Programmer          Description of change
!       ====       ==========          =====================
!     04/05/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=6) ! Precision
REAL(KIND=kind), PARAMETER :: twopi = 6.283185307179586_kind ! 2*pi

! Declare calling arguments:
INTEGER, INTENT(IN) :: n             ! No. of elements in IDFT
COMPLEX(KIND=kind),INTENT(IN),DIMENSION(0:n-1) :: array_in
                                     ! Input data for IDFT
COMPLEX(KIND=kind),INTENT(INOUT),DIMENSION(0:n-1) :: array_out
                                     ! Output data for IDFT

! List of local variables:
INTEGER :: i, k                      ! Index variables

! Calculate the inverse DFT from its definition.
DO i = 0, n-1
   array_out(i) = (0.,0.)
   DO k = 0, n-1
      array_out(i) = array_out(i) + array_in(k) &
                   * EXP( twopi*(0._kind,1._kind)*k*i/n )
   END DO
   array_out(i) = array_out(i) / REAL(n,KIND=kind)
END DO

END SUBROUTINE idft_sgl

SUBROUTINE idft_dbl ( array_in, array_out, n )
!
!  Purpose:
!    An implementation of the inverse Discrete Fourier Transform
!    from the defining equation.
!
!  Record of revisions:
!       Date       Programmer          Description of change
!       ====       ==========          =====================
!     04/05/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=12) ! Precision
REAL(KIND=kind), PARAMETER :: twopi = 6.283185307179586_kind ! 2*pi

! Declare calling arguments:
INTEGER, INTENT(IN) :: n             ! No. of elements in IDFT
COMPLEX(KIND=kind),INTENT(IN),DIMENSION(0:n-1) :: array_in
                                     ! Input data for IDFT
COMPLEX(KIND=kind),INTENT(INOUT),DIMENSION(0:n-1) :: array_out
                                     ! Output data for IDFT

! List of local variables:
INTEGER :: i, k                      ! Index variables

! Calculate the inverse DFT from its definition.
DO i = 0, n-1
   array_out(i) = (0.,0.)
   DO k = 0, n-1
      array_out(i) = array_out(i) + array_in(k) &
                   * EXP( twopi*(0._kind,1._kind)*k*i/n )
   END DO
   array_out(i) = array_out(i) / REAL(n,KIND=kind)
END DO

END SUBROUTINE idft_dbl

SUBROUTINE ucase ( string )
!
!  Purpose:
!    To shift a character string to upper case on any processor,
!    regardless of collating sequence.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    11/25/95    S. J. Chapman        Original code
! 1. 11/25/95    S. J. Chapman        Modified for any collating seq
!
IMPLICIT NONE

! Declare calling parameters:
CHARACTER(len=*), INTENT(INOUT) :: string

! Declare local variables:
INTEGER :: i                 ! Loop index
INTEGER :: length            ! Length of input string

! Get length of string
length = LEN ( string )

! Now shift lower case letters to upper case.
DO i = 1, length
   IF ( ( string(i:i) >= 'a' ) .AND.  &
        ( string(i:i) <= 'z' ) ) THEN
      string(i:i) = ACHAR ( IACHAR ( string(i:i) ) - 32 )
   END IF
END DO

END SUBROUTINE ucase

SUBROUTINE lcase ( string )
!
!  Purpose:
!    To shift a character string to lower case on any processor,
!    regardless of collating sequence.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    11/27/95    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare calling parameters:
CHARACTER(len=*), INTENT(INOUT) :: string

! Declare local variables:
INTEGER :: i                 ! Loop index
INTEGER :: length            ! Length of input string

! Get length of string
length = LEN ( string )

! Now shift upper case letters to lower case.
DO i = 1, length
   IF ( ( string(i:i) >= 'A' ) .AND.  &
        ( string(i:i) <= 'Z' ) ) THEN
      string(i:i) = ACHAR ( IACHAR ( string(i:i) ) + 32 )
   END IF
END DO

END SUBROUTINE lcase

SUBROUTINE integrate_sgl ( f, x1, x2, dx, area, error )
!
!  Purpose:
!    To integrate function f(x) between x1 and x2 using
!    rectangles of width dx to approximate the area
!    under the curve f(x).
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    02/15/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=6) ! Precision

! Declare calling arguments:
REAL(KIND=kind), EXTERNAL :: f        ! Function to integrate
REAL(KIND=kind), INTENT(IN) :: x1     ! Starting point for integral
REAL(KIND=kind), INTENT(IN) :: x2     ! Ending point for integral
REAL(KIND=kind), INTENT(IN) :: dx     ! Step size
REAL(KIND=kind), INTENT(OUT) :: area  ! Area under the curve.
INTEGER, INTENT(OUT) :: error         ! Error flag:
                                      !   0 = No error.
                                      !   1 = x1 > x2.

! Declare local variables:
REAL(KIND=kind) :: height       ! Height of rectangle
INTEGER :: i                    ! Index variable
INTEGER :: n                    ! Number of rectangles to integrate
REAL(KIND=kind) :: width        ! Width of rectangle
REAL(KIND=kind) :: xstart       ! Starting position of rectangle

! First, check to make sure that x1 <= x2.
errchk: IF ( x1 > x2 ) THEN
   error = 1
ELSE
   ! Clear error flag and area.
   error = 0
   area = 0

   ! Calculate the number of intervals to use.
   n = INT( (x2-x1) / dx + 1. )

   ! Calculate and sum the areas of each rectangle.
   sum: DO i = 1, n
      xstart = x1 + REAL(i-1) * dx
      width  = MIN ( dx, x2 - xstart )
      height = f( xstart + width/2. )
      area   = area + width * height
   END DO sum
END IF errchk

END SUBROUTINE integrate_sgl

SUBROUTINE integrate_dbl ( f, x1, x2, dx, area, error )
!
!  Purpose:
!    To integrate function f(x) between x1 and x2 using
!    rectangles of width dx to approximate the area
!    under the curve f(x).
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    02/15/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=12) ! Precision

! Declare calling arguments:
REAL(KIND=kind), EXTERNAL :: f        ! Function to integrate
REAL(KIND=kind), INTENT(IN) :: x1     ! Starting point for integral
REAL(KIND=kind), INTENT(IN) :: x2     ! Ending point for integral
REAL(KIND=kind), INTENT(IN) :: dx     ! Step size
REAL(KIND=kind), INTENT(OUT) :: area  ! Area under the curve.
INTEGER, INTENT(OUT) :: error         ! Error flag:
                                      !   0 = No error.
                                      !   1 = x1 > x2.

! Declare local variables:
REAL(KIND=kind) ::…27878 tokens truncated…s>1
if (nthreads > 1) then
  write(*,*)'sorting output ...'
  if (scratch == '') then
    call msort() !sorting files in memory
  else
    call fsort() !sorting files on disk
  endif
endif


end do !loop over multiple control (input nml) files

call quit(status)

end program ferre

