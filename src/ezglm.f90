!#####################################################################
module program_constants
   ! Programming constants used throughout the ELOGIT program.
   ! Unlike most modules, everything here is public.
   implicit none
   public
   ! Define compiler-specific KIND numbers for integers,
   ! single and double-precision reals to help ensure consistency of
   ! performance across platforms:
   integer, parameter :: our_int = selected_int_kind(9), &
        our_sgle = selected_real_kind(6,37), &
        our_dble = SELECTED_REAL_KIND(12,80)
   ! Define UNIT numbers for Fortran I/O:
   integer, parameter :: ctrl_file_handle = 11, &
        data_file_handle = 12, names_file_handle = 13, &
        out_file_handle = 13
   ! Define maximum lengths for various types of character strings:
   integer, parameter :: file_name_length = 256, &
        var_name_length = 8, case_id_length = 8
   ! Define the maximum line widths for various types of files:
   integer, parameter :: ctrl_line_width = 80, &
        data_line_width = 1024, names_line_width = 80, &
        out_line_width = 70
   ! Common integer values returned by all functions to indicate
   ! success or failure:
   integer(kind=our_int), parameter :: RETURN_SUCCESS = 0, &
        RETURN_FAIL = -1
   ! Character strings describing this program:
   character(len=*), parameter :: &
        program_name = "ELOGIT", &
        program_description = &
        "A simple program for logistic regression analysis", &
        program_version = "Version 1.0", &
        program_version_and_date = "Version 1.0 - June, 2004", &
        program_author = "Written by J.L. Schafer", &
        program_institution_1 = &
        "Department of Statistics and The Methodology Center", &
        program_institution_2 = "The Pennsylvania State University"
end module program_constants
!#####################################################################
function alnorm ( x, upper )

!*****************************************************************************80
!
!! ALNORM computes the cumulative density of the standard normal distribution.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by David Hill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Hill,
!    Algorithm AS 66:
!    The Normal Integral,
!    Applied Statistics,
!    Volume 22, Number 3, 1973, pages 424-427.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, is one endpoint of the semi-infinite interval
!    over which the integration takes place.
!
!    Input, logical UPPER, determines whether the upper or lower
!    interval is to be integrated:
!    .TRUE.  => integrate from X to + Infinity;
!    .FALSE. => integrate from - Infinity to X.
!
!    Output, real ( kind = 8 ) ALNORM, the integral of the standard normal
!    distribution over the desired interval.
!
  implicit none

  real ( kind = 8 ), parameter :: a1 = 5.75885480458D+00
  real ( kind = 8 ), parameter :: a2 = 2.62433121679D+00
  real ( kind = 8 ), parameter :: a3 = 5.92885724438D+00
  real ( kind = 8 ) alnorm
  real ( kind = 8 ), parameter :: b1 = -29.8213557807D+00
  real ( kind = 8 ), parameter :: b2 = 48.6959930692D+00
  real ( kind = 8 ), parameter :: c1 = -0.000000038052D+00
  real ( kind = 8 ), parameter :: c2 = 0.000398064794D+00
  real ( kind = 8 ), parameter :: c3 = -0.151679116635D+00
  real ( kind = 8 ), parameter :: c4 = 4.8385912808D+00
  real ( kind = 8 ), parameter :: c5 = 0.742380924027D+00
  real ( kind = 8 ), parameter :: c6 = 3.99019417011D+00
  real ( kind = 8 ), parameter :: con = 1.28D+00
  real ( kind = 8 ), parameter :: d1 = 1.00000615302D+00
  real ( kind = 8 ), parameter :: d2 = 1.98615381364D+00
  real ( kind = 8 ), parameter :: d3 = 5.29330324926D+00
  real ( kind = 8 ), parameter :: d4 = -15.1508972451D+00
  real ( kind = 8 ), parameter :: d5 = 30.789933034D+00
  real ( kind = 8 ), parameter :: ltone = 7.0D+00
  real ( kind = 8 ), parameter :: p = 0.398942280444D+00
  real ( kind = 8 ), parameter :: q = 0.39990348504D+00
  real ( kind = 8 ), parameter :: r = 0.398942280385D+00
  logical up
  logical upper
  real ( kind = 8 ), parameter :: utzero = 18.66D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  up = upper
  z = x

  if ( z < 0.0D+00 ) then
    up = .not. up
    z = - z
  end if

  if ( ltone < z .and. ( ( .not. up ) .or. utzero < z ) ) then

    if ( up ) then
      alnorm = 0.0D+00
    else
      alnorm = 1.0D+00
    end if

    return

  end if

  y = 0.5D+00 * z * z

  if ( z <= con ) then

    alnorm = 0.5D+00 - z * ( p - q * y &
      / ( y + a1 + b1 &
      / ( y + a2 + b2 &
      / ( y + a3 ))))

  else

    alnorm = r * exp ( - y ) &
      / ( z + c1 + d1 &
      / ( z + c2 + d2 &
      / ( z + c3 + d3 &
      / ( z + c4 + d4 &
      / ( z + c5 + d5 &
      / ( z + c6 ))))))

  end if

  if ( .not. up ) then
    alnorm = 1.0D+00 - alnorm
  end if

  return
end
subroutine normal_01_cdf_values ( n_data, x, fx )

!*****************************************************************************80
!
!! NORMAL_01_CDF_VALUES returns some values of the Normal 01 CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = NormalDistribution [ 0, 1 ]
!      CDF [ dist, x ]
!
!  Modified:
!
!    28 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 17

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.5000000000000000D+00, &
    0.5398278372770290D+00, &
    0.5792597094391030D+00, &
    0.6179114221889526D+00, &
    0.6554217416103242D+00, &
    0.6914624612740131D+00, &
    0.7257468822499270D+00, &
    0.7580363477769270D+00, &
    0.7881446014166033D+00, &
    0.8159398746532405D+00, &
    0.8413447460685429D+00, &
    0.9331927987311419D+00, &
    0.9772498680518208D+00, &
    0.9937903346742239D+00, &
    0.9986501019683699D+00, &
    0.9997673709209645D+00, &
    0.9999683287581669D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0000000000000000D+00, &
    0.1000000000000000D+00, &
    0.2000000000000000D+00, &
    0.3000000000000000D+00, &
    0.4000000000000000D+00, &
    0.5000000000000000D+00, &
    0.6000000000000000D+00, &
    0.7000000000000000D+00, &
    0.8000000000000000D+00, &
    0.9000000000000000D+00, &
    0.1000000000000000D+01, &
    0.1500000000000000D+01, &
    0.2000000000000000D+01, &
    0.2500000000000000D+01, &
    0.3000000000000000D+01, &
    0.3500000000000000D+01, &
    0.4000000000000000D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine normp ( z, p, q, pdf )

!*****************************************************************************80
!
!! NORMP computes the cumulative density of the standard normal distribution.
!
!  Discussion:
!
!    This is algorithm 5666 from Hart, et al.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by Alan Miller.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thacher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) Z, divides the real line into two
!    semi-infinite intervals, over each of which the standard normal
!    distribution is to be integrated.
!
!    Output, real ( kind = 8 ) P, Q, the integrals of the standard normal
!    distribution over the intervals ( - Infinity, Z] and
!    [Z, + Infinity ), respectively.
!
!    Output, real ( kind = 8 ) PDF, the value of the standard normal
!    distribution at Z.
!
  implicit none

  real ( kind = 8 ) :: cutoff = 7.071D+00
  real ( kind = 8 ) expntl
  real ( kind = 8 ) p
  real ( kind = 8 ) :: p0 = 220.2068679123761D+00
  real ( kind = 8 ) :: p1 = 221.2135961699311D+00
  real ( kind = 8 ) :: p2 = 112.0792914978709D+00
  real ( kind = 8 ) :: p3 = 33.91286607838300D+00
  real ( kind = 8 ) :: p4 = 6.373962203531650D+00
  real ( kind = 8 ) :: p5 = 0.7003830644436881D+00
  real ( kind = 8 ) :: p6 = 0.03526249659989109D+00
  real ( kind = 8 ) pdf
  real ( kind = 8 ) q
  real ( kind = 8 ) :: q0 = 440.4137358247522D+00
  real ( kind = 8 ) :: q1 = 793.8265125199484D+00
  real ( kind = 8 ) :: q2 = 637.3336333788311D+00
  real ( kind = 8 ) :: q3 = 296.5642487796737D+00
  real ( kind = 8 ) :: q4 = 86.78073220294608D+00
  real ( kind = 8 ) :: q5 = 16.06417757920695D+00
  real ( kind = 8 ) :: q6 = 1.755667163182642D+00
  real ( kind = 8 ) :: q7 = 0.08838834764831844D+00
  real ( kind = 8 ) :: root2pi = 2.506628274631001D+00
  real ( kind = 8 ) z
  real ( kind = 8 ) zabs

  zabs = abs ( z )
!
!  37 < |Z|.
!
  if ( 37.0D+00 < zabs ) then

    pdf = 0.0D+00
    p = 0.0D+00
!
!  |Z| <= 37.
!
  else

    expntl = exp ( - 0.5D+00 * zabs * zabs )
    pdf = expntl / root2pi
!
!  |Z| < CUTOFF = 10 / sqrt(2).
!
    if ( zabs < cutoff ) then

      p = expntl * (((((( &
          p6   * zabs &
        + p5 ) * zabs &
        + p4 ) * zabs &
        + p3 ) * zabs &
        + p2 ) * zabs &
        + p1 ) * zabs &
        + p0 ) / ((((((( &
          q7   * zabs &
        + q6 ) * zabs &
        + q5 ) * zabs &
        + q4 ) * zabs &
        + q3 ) * zabs &
        + q2 ) * zabs &
        + q1 ) * zabs &
      + q0 )
!
!  CUTOFF <= |Z|.
!
    else

      p = pdf / ( &
        zabs + 1.0D+00 / ( &
        zabs + 2.0D+00 / ( &
        zabs + 3.0D+00 / ( &
        zabs + 4.0D+00 / ( &
        zabs + 0.65D+00 )))))

    end if

  end if

  if ( z < 0.0D+00 ) then
    q = 1.0D+00 - p
  else
    q = p
    p = 1.0D+00 - q
  end if

  return
end
subroutine nprob ( z, p, q, pdf )

!*****************************************************************************80
!
!! NPROB computes the cumulative density of the standard normal distribution.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by AG Adams.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    AG Adams,
!    Algorithm 39:
!    Areas Under the Normal Curve,
!    Computer Journal,
!    Volume 12, Number 2, May 1969, pages 197-198.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) Z, divides the real line into
!    two semi-infinite intervals, over each of which the standard normal
!    distribution is to be integrated.
!
!    Output, real ( kind = 8 ) P, Q, the integrals of the standard normal
!    distribution over the intervals ( - Infinity, Z] and
!    [Z, + Infinity ), respectively.
!
!    Output, real ( kind = 8 ) PDF, the value of the standard normal
!    distribution at Z.
!
  implicit none

  real ( kind = 8 ), parameter :: a0 = 0.5D+00
  real ( kind = 8 ), parameter :: a1 = 0.398942280444D+00
  real ( kind = 8 ), parameter :: a2 = 0.399903438504D+00
  real ( kind = 8 ), parameter :: a3 = 5.75885480458D+00
  real ( kind = 8 ), parameter :: a4 = 29.8213557808D+00
  real ( kind = 8 ), parameter :: a5 = 2.62433121679D+00
  real ( kind = 8 ), parameter :: a6 = 48.6959930692D+00
  real ( kind = 8 ), parameter :: a7 = 5.92885724438D+00
  real ( kind = 8 ), parameter :: b0 = 0.398942280385D+00
  real ( kind = 8 ), parameter :: b1 = 0.000000038052D+00
  real ( kind = 8 ), parameter :: b2 = 1.00000615302D+00
  real ( kind = 8 ), parameter :: b3 = 0.000398064794D+00
  real ( kind = 8 ), parameter :: b4 = 1.98615381364D+00
  real ( kind = 8 ), parameter :: b5 = 0.151679116635D+00
  real ( kind = 8 ), parameter :: b6 = 5.29330324926D+00
  real ( kind = 8 ), parameter :: b7 = 4.8385912808D+00
  real ( kind = 8 ), parameter :: b8 = 15.1508972451D+00
  real ( kind = 8 ), parameter :: b9 = 0.742380924027D+00
  real ( kind = 8 ), parameter :: b10 = 30.789933034D+00
  real ( kind = 8 ), parameter :: b11 = 3.99019417011D+00
  real ( kind = 8 ) p
  real ( kind = 8 ) pdf
  real ( kind = 8 ) q
  real ( kind = 8 ) y
  real ( kind = 8 ) z
  real ( kind = 8 ) zabs

  zabs = abs ( z )
!
!  |Z| between 0 and 1.28
!
  if ( abs ( z ) <= 1.28D+00 ) then

    y = a0 * z * z
    pdf = exp ( - y ) * b0

    q = a0 - zabs * ( a1 - a2 * y &
      / ( y + a3 - a4 &
      / ( y + a5 + a6 &
      / ( y + a7 ))))
!
!  |Z| between 1.28 and 12.7
!
  else if ( abs ( z ) <= 12.7D+00 ) then

    y = a0 * z * z
    pdf = exp ( - y ) * b0

    q = pdf &
      / ( zabs - b1 + b2 &
      / ( zabs + b3 + b4 &
      / ( zabs - b5 + b6 &
      / ( zabs + b7 - b8 &
      / ( zabs + b9 + b10 &
      / ( zabs + b11 ))))))
!
!  Z far out in tail.
!
  else

    q = 0.0D+00
    pdf = 0.0D+00

  end if

  if ( z < 0.0D+00 ) then
    p = q
    q = 1.0D+00 - p
  else
    p = 1.0D+00 - q
  end if

  return
end



!##################################################################### 
module error_handler
   ! Generic error message handler for both console and non-console
   ! applications. The routines in this module do not halt program
   ! execution; they merely store error messages for subsequent
   ! retrieval.
   ! Written by J.L. Schafer, 4/19/02; revised 6/03/03
   implicit none
   private ! by default
   ! declare public types
   public :: error_type
   ! declare public subroutines and functions
   public :: err_reset, err_handle, err_msg_present, err_get_msgs
   ! Parameters private to this module
   integer, parameter :: &
        ! max width of any single error message line
        err_msg_width = 70
   !##################################################################
   type :: msg_line_type
      ! Private type for a single node in the linked list
      sequence
      character (len=err_msg_width) :: line = ""
      type(msg_line_type), pointer :: next => null()
   end type msg_line_type
   !##################################################################
   type :: error_type
      ! Public type for holding a linked list of messages
      sequence
      private ! contents of this type are private to this module
      logical :: msg_present = .false.
      type(msg_line_type), pointer :: head => null(), tail => null()
   end type error_type 
   !##################################################################
contains
   !##################################################################
   subroutine insert_msg_line( text_line, err )
      ! inserts a new message line at the end of the list
      implicit none
      ! declare arguments
      character(len=*), intent(in) :: text_line
      type(error_type), intent(inout) :: err
      ! begin
      if( .not. err%msg_present ) then
         ! begin a linked list
         allocate( err%head )
         err%tail => err%head
         err%head%line = text_line
         err%msg_present = .true.
      else
         ! add a node to the list; point tail to the new node
         allocate( err%tail%next )
         err%tail => err%tail%next
         err%tail%line = text_line
      end if
   end subroutine insert_msg_line
   !##################################################################
   subroutine err_reset( err )
      ! Public: deletes all messages from the list
      implicit none
      type(error_type), intent(inout) :: err
      type(msg_line_type), pointer :: current_line
      if( .not. err%msg_present ) return
      do
         current_line => err%head
         err%head => err%head%next
         deallocate( current_line )
         if( .not.associated(err%head) ) exit
      end do
      nullify(err%tail)
      err%msg_present = .false.
   end subroutine err_reset
   !##################################################################
   subroutine err_handle(err, err_code, called_from, file_name, &
        line_no, object_name, custom_1, custom_2, custom_3)
      ! Public: Stores a message in the error handler
      ! Meaning of err_code
      !     0 = no error
      ! 1-99: I/O errors
      !     1 = file could not be opened for read-access
      !     2 = file could not be opened for write-access
      !     3 = error in reading file
      !     4 = error in writing to file
      !     5 = error in reading file: EOR/EOF encountered
      !     6 = file open for write-access could not be closed
      !     7 = file open for read-access could not be closed
      ! 100-199: numerical errors
      !   100 = matrix apparently singular
      !   101 = matrix not positive definite
      !   102 = attempted division by zero
      !   103 = attempted logarithm of non-positive number
      !   104 = argument to exp function too large
      !   105 = attempted square root of negative number
      ! 200-299: memory errors
      !   200 = unable to dynamically allocate memory for object
      !   201 = unable to deallocate memory for object
      ! 300-399: array dimension errors
      !   300 = non-square matrix encountered where square
      !     matrix expected
      !   301 = dimensions of matrix arguments not conformable
      ! 1000: other error
      !  1000 = reserved for custom error messages
      implicit none
      ! declare required arguments
      type(error_type), intent(inout) :: err
      integer, intent(in) :: err_code
      ! declare optional arguments
      character (len=*), optional :: called_from
      character (len=*), optional :: file_name
      character (len=*), optional :: object_name
      character (len=*), optional :: custom_1
      character (len=*), optional :: custom_2
      character (len=*), optional :: custom_3
      integer, optional :: line_no
      ! local variables
      character(len=12) :: ichar
      ! begin
      select case(err_code)
         case(0)
            call insert_msg_line( &
                 "No errors", err)
         ! I/O errors
         case(1)
            call insert_msg_line( &
                 "File could not be opened for read-access", err)
         case(2)
            call insert_msg_line( &
                 "File could not be opened for write-access", err)
         case(3)
            call insert_msg_line( &
                 "Error in reading file", err)
         case(4)
            call insert_msg_line( &
                 "Error in writing to file", err)
         case(5)
            call insert_msg_line( &
                 "Error in reading file: EOR or EOF encountered", err)
         case(6)
            call insert_msg_line( &
                 "File open for write-access could not be closed", &
                 err)
         case(7)
            call insert_msg_line( &
                 "File open for read-access could not be closed", err)
         ! numerical errors
         case(100)
            call insert_msg_line( &
                 "Matrix apparently singular", err)
         case(101)
            call insert_msg_line(  &
                 "Matrix not positive definite", err)
         case(102)
            call insert_msg_line(  &
                 "Attempted division by zero", err)
         case(103)
            call insert_msg_line( &
                 "Attempted logarithm of non-positive number", err)
         case(104)
            call insert_msg_line( &
                 "Argument to exp function too large", err)
         case(105)
            call insert_msg_line( &
                 "Attempted square root of negative number", err)
         ! memory errors
         case(200)
            call insert_msg_line(  &
                 "Unable to allocate memory for object", err)
         case(201)
            call insert_msg_line( &
                 "Unable to deallocate memory for object", err)
         ! array dimension errors
         case(300)
            call insert_msg_line( &
                 "Non-square matrix encountered; square expected", &
                 err)
         case(301)
            call insert_msg_line( &
                 "Dimensions of matrix arguments not conformable", &
                 err)
         ! custom error message
         case(1000)
            ! don't do anything yet
         ! anything else
         case default
            call insert_msg_line("Unknown error code.", err)
      end select
      ! append other optional information if present
      if( present(custom_1) ) then
         call insert_msg_line(custom_1, err)
      end if
      if( present(custom_2) ) then
         call insert_msg_line(custom_2, err)
      end if
      if( present(custom_3) ) then
         call insert_msg_line(custom_3, err)
      end if
      if(present(file_name)) then
         call insert_msg_line("FILE: " // trim(file_name), err)
      end if
      if(present(line_no)) then
         ichar = adjustl(ichar)
         call insert_msg_line("LINE: " // trim(ichar), err)
      end if
      if(present(object_name)) then
         call insert_msg_line(trim(object_name), err)
      end if
      if(present(called_from)) then
         call insert_msg_line("OCCURRED IN: " // &
              trim(called_from), err)
      end if
   end subroutine err_handle
   !##################################################################
   logical function err_msg_present(err)
      ! Public: Queries the error_type to see if a message is present
      implicit none
      type(error_type), intent(inout) :: err
      err_msg_present = err%msg_present
   end function err_msg_present
   !##################################################################
   subroutine err_get_msgs(err, msg_string, platform)
      ! Public: Retrieves all stored messages as a single character
      ! string, with message lines separated by platform-appropriate
      ! ASCII carriage control characters.
      ! Values for platform may be "UNIX", "MAC" or "PC"
      implicit none
      ! required arguments
      type(error_type), intent(inout) :: err
      character(len=*), intent(out) :: msg_string
      ! optional arguments
      character(len=*), intent(in), optional :: platform
      ! local variables
      character(len=4) :: plat
      integer :: posn
      logical :: first_time
      type(msg_line_type), pointer :: cur_line
      ! determine platform
      if( present(platform) ) then
         plat = platform
      else
         plat = "PC"
      end if
      ! clean out msg_string
      msg_string = ""
      ! step through the linked list, appending the lines
      first_time = .true.
      cur_line => err%head
      do
         if( .not.associated(cur_line) ) exit
         posn = len_trim(msg_string)
         if( (posn+3) >= len(msg_string) ) exit ! out of space
         posn = posn + 1
         if( .not.first_time) then
            select case(plat)
               case("UNIX")
                  ! Separate lines with LF
                  msg_string(posn:) = achar(10)
                  posn = posn + 1
               case("MAC")
                  ! Separate lines with CR
                  msg_string(posn:) = achar(13)
                  posn = posn + 1
               case default
                  msg_string(posn:) = achar(13) // achar(10)
                  posn = posn + 2
            end select
         end if
         msg_string(posn:) = trim(cur_line%line)
         first_time = .false.
         cur_line => cur_line%next
      end do
   end subroutine err_get_msgs
   !##################################################################
end module error_handler
!#####################################################################

!#######################################################################
module matrix_methods
   use program_constants
   use error_handler
   implicit none
   private ! by default
   public :: cholesky_saxpy, invert_lower, premult_lower_by_transpose
   character(len=*), parameter :: modname = "matrix_methods"
   !####################################################################
contains 
   !####################################################################
   subroutine cholesky_saxpy(a, err)
      !### Overwrites lower triangle of a symmetric, pos.-def.
      !### matrix a with its cholesky factor.
      implicit none
      ! declare arguments
      real(kind=our_dble), intent(inout) :: a(:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      character(len=*), parameter :: subname = "cholesky_saxpy"
      integer(kind=our_int) :: p, j, k
      real(kind=our_dble) :: den
      ! begin
      p = size(a,1)
      if( p /= size(a,2) ) goto 700
      do j = 1, p
         do k = 1, j-1
            a(j:p,j) = a(j:p,j) - a(j:p,k) * a(j,k)
         end do
         if( a(j,j) <= 0.D0 ) goto 710
         den = sqrt( a(j,j) )
         a(j:p,j) = a(j:p,j) / den
      end do
      ! normal exit
      return
      ! error traps
700   call err_handle(err, 300, &
           called_from = subname//" in MOD "//modname)
      return
710   call err_handle(err, 101, &
           called_from = subname//" in MOD "//modname)
      return
   end subroutine cholesky_saxpy
   !####################################################################
   subroutine invert_lower(a, err)
      !### Overwrites a lower-triangular matrix a with its inverse
      !### by forward substitution. The upper triangle is untouched.
      implicit none
      ! declare arguments
      real(kind=our_dble), intent(inout) :: a(:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      character(len=*), parameter :: subname = "invert_lower"
      integer(kind=our_int) :: p, i, j, k
      real(kind=our_dble) :: sum
      ! begin
      p = size(a,1)
      if( p /= size(a,2) ) goto 700
      if( a(1,1) == 0.D0 ) goto 710
      a(1,1) = 1.D0 / a(1,1)
      do i = 2, p
         if( a(i,i) == 0.D0 ) goto 710
         a(i,i) = 1.D0 / a(i,i)
         do j = 1, i-1
            sum = 0.D0
            do k = j, i-1
               sum = sum + a(k,j) * a(i,k)
            end do
            a(i,j) = - sum * a(i,i)
         end do
      end do
      ! normal exit
      return
      ! error traps
700   call err_handle(err, 300, &
           called_from = subname//" in MOD "//modname)
      return
710   call err_handle(err, 100, &
           called_from = subname//" in MOD "//modname)
      return
   end subroutine invert_lower
   !####################################################################
   subroutine premult_lower_by_transpose(a, b, err)
      !### Premultiplies a lower-triangular matrix a by its upper-
      !### triangular transpose to produce a symmetric matrix b.
      implicit none
      ! declare arguments
      real(kind=our_dble), intent(in) :: a(:,:)
      real(kind=our_dble), intent(out) :: b(:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      character(len=*), parameter :: subname = &
           "premult_lower_by_transpose"
      integer(kind=our_int) :: p, i, j, k
      ! begin
      p = size(a,1)
      if( p /= size(a,2) ) goto 700
      if( ( p /= size(b,1) ) .or. ( p /= size(b,2) ) ) goto 710
      do i = 1, p
         do j = 1, i
            b(i,j) = 0.D0
            do k = max(i,j), p   ! skip zero elements
               b(i,j) = b(i,j) + a(k,i) * a(k,j)
            end do
            b(j,i) = b(i,j)      ! copy upper triangle from lower
         end do
      end do
      ! normal exit
      return
      ! error traps 
700   call err_handle(err, 300, &
           called_from = subname//" in MOD "//modname)
      return
710   call err_handle(err, 301, &
           called_from = subname//" in MOD "//modname)
      return
  end subroutine premult_lower_by_transpose
   !####################################################################
end module matrix_methods
!#######################################################################
!#####################################################################
module dynamic_allocation
   ! Routines for allocating and deallocating pointers to arrays
   use error_handler
   use program_constants
   implicit none
   private ! by default
   public :: dyn_alloc, dyn_dealloc
   ! For allocating pointers to arrays
   interface dyn_alloc
      module procedure int1_alloc
      module procedure int2_alloc
      module procedure int3_alloc
      module procedure int4_alloc
      module procedure int5_alloc
      module procedure int6_alloc
      module procedure dbl1_alloc
      module procedure dbl2_alloc
      module procedure dbl3_alloc
      module procedure dbl4_alloc
      module procedure dbl5_alloc
      module procedure dbl6_alloc
      module procedure str1_alloc
      module procedure str2_alloc
      module procedure str3_alloc
      module procedure str4_alloc
      module procedure str5_alloc
      module procedure str6_alloc
   end interface
   ! For deallocating pointers to arrays
   interface dyn_dealloc
      module procedure int1_dealloc
      module procedure int2_dealloc
      module procedure int3_dealloc
      module procedure int4_dealloc
      module procedure int5_dealloc
      module procedure int6_dealloc
      module procedure dbl1_dealloc
      module procedure dbl2_dealloc
      module procedure dbl3_dealloc
      module procedure dbl4_dealloc
      module procedure dbl5_dealloc
      module procedure dbl6_dealloc
      module procedure str1_dealloc
      module procedure str2_dealloc
      module procedure str3_dealloc
      module procedure str4_dealloc
      module procedure str5_dealloc
      module procedure str6_dealloc
   end interface
   ! parameters private to this module
   character(len=*), parameter :: modname = "dynalloc"
   !##################################################################
contains
   !##################################################################
   integer(our_int) function int1_alloc(intArray, dim1, err, &
        lbound1) result(answer)
      ! Allocates an integer array of rank 1
      implicit none
      ! declare required arguments
      integer(kind=our_int), pointer :: intArray(:)
      integer, intent(in) :: dim1
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1
      ! declare local variables and parameters
      integer :: lb1, status
      character(len=*), parameter :: subname = "int1_alloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(intArray) ) deallocate(intArray, stat=status)
      if( status /= 0) goto 800
      lb1 = 1
      if( present(lbound1) ) lb1 = lbound1
      allocate( intArray(lb1:dim1), stat=status )
      if( status /= 0) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
810   call err_handle(err, 200, &
              called_from = subname//" in MOD "//modname)
      return
   end function int1_alloc
   !##################################################################
   integer(our_int) function int2_alloc(intArray, dim1, dim2, err, &
        lbound1, lbound2) result(answer)
      ! Allocates an integer array of rank 2
      implicit none
      ! declare required arguments
      integer(kind=our_int), pointer :: intArray(:,:)
      integer, intent(in) :: dim1, dim2
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2
      ! declare local variables and parameters
      integer :: status, lb1, lb2
      character(len=*), parameter :: subname = "int2_alloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(intArray) ) deallocate(intArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if(status == 0) allocate( intArray(lb1:dim1, lb2:dim2), &
           stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
810   call err_handle(err, 200, &
              called_from = subname//" in MOD "//modname)
      return
   end function int2_alloc
   !##################################################################
   integer(our_int) function int3_alloc(intArray, dim1, dim2, dim3, &
        err, lbound1, lbound2, lbound3) result(answer)
      ! Allocates an integer array of rank 3
      implicit none
      ! declare required arguments
      integer(kind=our_int), pointer :: intArray(:,:,:)
      integer, intent(in) :: dim1, dim2, dim3
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2, lbound3
      ! declare local variables and parameters
      integer :: status, lb1, lb2, lb3
      character(len=*), parameter :: subname = "int3_alloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(intArray) ) deallocate(intArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      lb3 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if( present(lbound3) ) lb3 = lbound3
      if(status == 0) allocate( intArray(lb1:dim1, lb2:dim2, &
           lb3:dim3), stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
810   call err_handle(err, 200, &
              called_from = subname//" in MOD "//modname)
      return
   end function int3_alloc
   !##################################################################
   integer(our_int) function int4_alloc(intArray, dim1, dim2, dim3, &
        dim4, err, lbound1, lbound2, lbound3, lbound4) result(answer)
      ! Allocates an integer array of rank 3
      implicit none
      ! declare required arguments
      integer(kind=our_int), pointer :: intArray(:,:,:,:)
      integer, intent(in) :: dim1, dim2, dim3, dim4
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2, lbound3, &
         lbound4
      ! declare local variables and parameters
      integer :: status, lb1, lb2, lb3, lb4
      character(len=*), parameter :: subname = "int4_alloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(intArray) ) deallocate(intArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      lb3 = 1
      lb4 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if( present(lbound3) ) lb3 = lbound3
      if( present(lbound4) ) lb4 = lbound4
      if(status == 0) allocate( intArray(lb1:dim1, lb2:dim2, &
           lb3:dim3, lb4:dim4), stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
810   call err_handle(err, 200, &
              called_from = subname//" in MOD "//modname)
      return
   end function int4_alloc
   !##################################################################
   integer(our_int) function int5_alloc(intArray, dim1, dim2, dim3, &
        dim4, dim5, err, lbound1, lbound2, lbound3, lbound4, &
        lbound5) result(answer)
      ! Allocates an integer array of rank 3
      implicit none
      ! declare required arguments
      integer(kind=our_int), pointer :: intArray(:,:,:,:,:)
      integer, intent(in) :: dim1, dim2, dim3, dim4, dim5
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2, lbound3, &
         lbound4, lbound5
      ! declare local variables and parameters
      integer :: status, lb1, lb2, lb3, lb4, lb5
      character(len=*), parameter :: subname = "int5_alloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(intArray) ) deallocate(intArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      lb3 = 1
      lb4 = 1
      lb5 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if( present(lbound3) ) lb3 = lbound3
      if( present(lbound4) ) lb4 = lbound4
      if( present(lbound5) ) lb5 = lbound5
      if(status == 0) allocate( intArray(lb1:dim1, lb2:dim2, &
           lb3:dim3, lb4:dim4, lb5:dim5), stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
810   call err_handle(err, 200, &
              called_from = subname//" in MOD "//modname)
      return
   end function int5_alloc
   !##################################################################
   integer(our_int) function int6_alloc(intArray, dim1, dim2, dim3, &
        dim4, dim5, dim6, err, lbound1, lbound2, lbound3, lbound4, &
        lbound5, lbound6) result(answer)
      ! Allocates an integer array of rank 3
      implicit none
      ! declare required arguments
      integer(kind=our_int), pointer :: intArray(:,:,:,:,:,:)
      integer, intent(in) :: dim1, dim2, dim3, dim4, dim5, dim6
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2, lbound3, &
         lbound4, lbound5, lbound6
      ! declare local variables and parameters
      integer :: status, lb1, lb2, lb3, lb4, lb5, lb6
      character(len=*), parameter :: subname = "int6_alloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(intArray) ) deallocate(intArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      lb3 = 1
      lb4 = 1
      lb5 = 1
      lb6 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if( present(lbound3) ) lb3 = lbound3
      if( present(lbound4) ) lb4 = lbound4
      if( present(lbound5) ) lb5 = lbound5
      if( present(lbound6) ) lb6 = lbound6
      if(status == 0) allocate( intArray(lb1:dim1, lb2:dim2, &
           lb3:dim3, lb4:dim4, lb5:dim5, lb6:dim6), stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
810   call err_handle(err, 200, &
              called_from = subname//" in MOD "//modname)
      return
   end function int6_alloc
   !##################################################################
   integer(our_int) function dbl1_alloc(dblArray, dim1, err, &
        lbound1) result(answer)
      ! Allocates a double-precision real array of rank 1
      implicit none
      ! declare required arguments
      real(kind=our_dble), pointer :: dblArray(:)
      integer, intent(in) :: dim1
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1
      ! declare local variables and parameters
      integer :: lb1, status
      character(len=*), parameter :: subname = "dbl1_alloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(dblArray) ) deallocate(dblArray, stat=status)
      if( status /= 0) goto 800
      lb1 = 1
      if( present(lbound1) ) lb1 = lbound1
      allocate( dblArray(lb1:dim1), stat=status )
      if( status /= 0) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
810   call err_handle(err, 200, &
              called_from = subname//" in MOD "//modname)
      return
   end function dbl1_alloc
   !##################################################################
   integer(our_int) function dbl2_alloc(dblArray, dim1, dim2, err, &
        lbound1, lbound2) result(answer)
      ! Allocates a double-precision real array of rank 2
      implicit none
      ! declare required arguments
      real(kind=our_dble), pointer :: dblArray(:,:)
      integer, intent(in) :: dim1, dim2
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2
      ! declare local variables and parameters
      integer :: status, lb1, lb2
      character(len=*), parameter :: subname = "dbl2_alloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(dblArray) ) deallocate(dblArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if(status == 0) allocate( dblArray(lb1:dim1, lb2:dim2), &
           stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
810   call err_handle(err, 200, &
              called_from = subname//" in MOD "//modname)
      return
   end function dbl2_alloc
   !##################################################################
   integer(our_int) function dbl3_alloc(dblArray, dim1, dim2, dim3, &
        err, lbound1, lbound2, lbound3) result(answer)
      ! Allocates a double-precision real array of rank 3
      implicit none
      ! declare required arguments
      real(kind=our_dble), pointer :: dblArray(:,:,:)
      integer, intent(in) :: dim1, dim2, dim3
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2, lbound3
      ! declare local variables and parameters
      integer :: status, lb1, lb2, lb3
      character(len=*), parameter :: subname = "dbl3_alloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(dblArray) ) deallocate(dblArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      lb3 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if( present(lbound3) ) lb3 = lbound3
      if(status == 0) allocate( dblArray(lb1:dim1, lb2:dim2, &
           lb3:dim3), stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
810   call err_handle(err, 200, &
              called_from = subname//" in MOD "//modname)
      return
   end function dbl3_alloc
   !##################################################################
   integer(our_int) function dbl4_alloc(dblArray, dim1, dim2, dim3, &
        dim4, err, lbound1, lbound2, lbound3, lbound4) result(answer)
      ! Allocates a double-precision real array of rank 4
      implicit none
      ! declare required arguments
      real(kind=our_dble), pointer :: dblArray(:,:,:,:)
      integer, intent(in) :: dim1, dim2, dim3, dim4
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2, lbound3, &
         lbound4
      ! declare local variables and parameters
      integer :: status, lb1, lb2, lb3, lb4
      character(len=*), parameter :: subname = "dbl4_alloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(dblArray) ) deallocate(dblArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      lb3 = 1
      lb4 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if( present(lbound3) ) lb3 = lbound3
      if( present(lbound4) ) lb4 = lbound4
      if(status == 0) allocate( dblArray(lb1:dim1, lb2:dim2, &
           lb3:dim3, lb4:dim4), stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
810   call err_handle(err, 200, &
              called_from = subname//" in MOD "//modname)
      return
   end function dbl4_alloc
   !##################################################################
   integer(our_int) function dbl5_alloc(dblArray, dim1, dim2, dim3, &
        dim4, dim5, err, lbound1, lbound2, lbound3, lbound4, &
        lbound5) result(answer)
      ! Allocates a double-precision real array of rank 5
      implicit none
      ! declare required arguments
      real(kind=our_dble), pointer :: dblArray(:,:,:,:,:)
      integer, intent(in) :: dim1, dim2, dim3, dim4, dim5
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2, lbound3, &
         lbound4, lbound5
      ! declare local variables and parameters
      integer :: status, lb1, lb2, lb3, lb4, lb5
      character(len=*), parameter :: subname = "dbl5_alloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(dblArray) ) deallocate(dblArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      lb3 = 1
      lb4 = 1
      lb5 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if( present(lbound3) ) lb3 = lbound3
      if( present(lbound4) ) lb4 = lbound4
      if( present(lbound5) ) lb5 = lbound5
      if(status == 0) allocate( dblArray(lb1:dim1, lb2:dim2, &
           lb3:dim3, lb4:dim4, lb5:dim5), stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
810   call err_handle(err, 200, &
              called_from = subname//" in MOD "//modname)
      return
   end function dbl5_alloc
   !##################################################################
   integer(our_int) function dbl6_alloc(dblArray, dim1, dim2, dim3, &
        dim4, dim5, dim6, err, lbound1, lbound2, lbound3, lbound4, &
        lbound5, lbound6) result(answer)
      ! Allocates a double-precision real array of rank 5
      implicit none
      ! declare required arguments
      real(kind=our_dble), pointer :: dblArray(:,:,:,:,:,:)
      integer, intent(in) :: dim1, dim2, dim3, dim4, dim5, dim6
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2, lbound3, &
         lbound4, lbound5, lbound6
      ! declare local variables and parameters
      integer :: status, lb1, lb2, lb3, lb4, lb5, lb6
      character(len=*), parameter :: subname = "dbl6_alloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(dblArray) ) deallocate(dblArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      lb3 = 1
      lb4 = 1
      lb5 = 1
      lb6 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if( present(lbound3) ) lb3 = lbound3
      if( present(lbound4) ) lb4 = lbound4
      if( present(lbound5) ) lb5 = lbound5
      if( present(lbound6) ) lb6 = lbound6
      if(status == 0) allocate( dblArray(lb1:dim1, lb2:dim2, &
           lb3:dim3, lb4:dim4, lb5:dim5, lb6:dim6), stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
810   call err_handle(err, 200, &
              called_from = subname//" in MOD "//modname)
      return
   end function dbl6_alloc
   !##################################################################
   integer(our_int) function str1_alloc(strArray, dim1, err, &
        lbound1) result(answer)
      ! Allocates a character-string array of rank 1
      implicit none
      ! declare required arguments
      character(len=*), pointer :: strArray(:)
      integer, intent(in) :: dim1
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1
      ! declare local variables and parameters
      integer :: lb1, status
      character(len=*), parameter :: subname = "str1_alloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(strArray) ) deallocate(strArray, stat=status)
      if( status /= 0) goto 800
      lb1 = 1
      if( present(lbound1) ) lb1 = lbound1
      allocate( strArray(lb1:dim1), stat=status )
      if( status /= 0) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
810   call err_handle(err, 200, &
              called_from = subname//" in MOD "//modname)
      return
   end function str1_alloc
   !##################################################################
   integer(our_int) function str2_alloc(strArray, dim1, dim2, err, &
        lbound1, lbound2) result(answer)
      ! Allocates a character-string array of rank 2
      implicit none
      ! declare required arguments
      character(len=*), pointer :: strArray(:,:)
      integer, intent(in) :: dim1, dim2
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2
      ! declare local variables and parameters
      integer :: status, lb1, lb2
      character(len=*), parameter :: subname = "str2_alloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(strArray) ) deallocate(strArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if(status == 0) allocate( strArray(lb1:dim1, lb2:dim2), &
           stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
810   call err_handle(err, 200, &
              called_from = subname//" in MOD "//modname)
      return
   end function str2_alloc
   !##################################################################
   integer(our_int) function str3_alloc(strArray, dim1, dim2, dim3, &
        err, lbound1, lbound2, lbound3) result(answer)
      ! Allocates a character-string array of rank 3
      implicit none
      ! declare required arguments
      character(len=*), pointer :: strArray(:,:,:)
      integer, intent(in) :: dim1, dim2, dim3
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2, lbound3
      ! declare local variables and parameters
      integer :: status, lb1, lb2, lb3
      character(len=*), parameter :: subname = "str3_alloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(strArray) ) deallocate(strArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      lb3 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if( present(lbound3) ) lb3 = lbound3
      if(status == 0) allocate( strArray(lb1:dim1, lb2:dim2, &
           lb3:dim3), stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
810   call err_handle(err, 200, &
              called_from = subname//" in MOD "//modname)
      return
   end function str3_alloc
   !##################################################################
   integer(our_int) function str4_alloc(strArray, dim1, dim2, dim3, &
        dim4, err, lbound1, lbound2, lbound3, lbound4) result(answer)
      ! Allocates a character-string array of rank 4
      implicit none
      ! declare required arguments
      character(len=*), pointer :: strArray(:,:,:,:)
      integer, intent(in) :: dim1, dim2, dim3, dim4
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2, lbound3, &
         lbound4
      ! declare local variables and parameters
      integer :: status, lb1, lb2, lb3, lb4
      character(len=*), parameter :: subname = "str4_alloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(strArray) ) deallocate(strArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      lb3 = 1
      lb4 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if( present(lbound3) ) lb3 = lbound3
      if( present(lbound4) ) lb4 = lbound4
      if(status == 0) allocate( strArray(lb1:dim1, lb2:dim2, &
           lb3:dim3, lb4:dim4), stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
810   call err_handle(err, 200, &
              called_from = subname//" in MOD "//modname)
      return
   end function str4_alloc
   !##################################################################
   integer(our_int) function str5_alloc(strArray, dim1, dim2, dim3, &
        dim4, dim5, err, lbound1, lbound2, lbound3, lbound4, lbound5 &
        ) result(answer)
      ! Allocates a character-string array of rank 5
      implicit none
      ! declare required arguments
      character(len=*), pointer :: strArray(:,:,:,:,:)
      integer, intent(in) :: dim1, dim2, dim3, dim4, dim5
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2, lbound3, &
         lbound4, lbound5
      ! declare local variables and parameters
      integer :: status, lb1, lb2, lb3, lb4, lb5
      character(len=*), parameter :: subname = "str5_alloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(strArray) ) deallocate(strArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      lb3 = 1
      lb4 = 1
      lb5 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if( present(lbound3) ) lb3 = lbound3
      if( present(lbound4) ) lb4 = lbound4
      if( present(lbound5) ) lb5 = lbound5
      if(status == 0) allocate( strArray(lb1:dim1, lb2:dim2, &
           lb3:dim3, lb4:dim4, lb5:dim5), stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
810   call err_handle(err, 200, &
              called_from = subname//" in MOD "//modname)
      return
   end function str5_alloc
   !##################################################################
   integer(our_int) function str6_alloc(strArray, dim1, dim2, dim3, &
        dim4, dim5, dim6, err, lbound1, lbound2, lbound3, lbound4, &
        lbound5, lbound6) result(answer)
      ! Allocates a character-string array of rank 6
      implicit none
      ! declare required arguments
      character(len=*), pointer :: strArray(:,:,:,:,:,:)
      integer, intent(in) :: dim1, dim2, dim3, dim4, dim5, dim6
      type(error_type), intent(inout) :: err
      ! declare optional arguments
      integer, intent(in), optional :: lbound1, lbound2, lbound3, &
         lbound4, lbound5, lbound6
      ! declare local variables and parameters
      integer :: status, lb1, lb2, lb3, lb4, lb5, lb6
      character(len=*), parameter :: subname = "str6_alloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(strArray) ) deallocate(strArray, stat=status)
      if( status /= 0 ) goto 800
      lb1 = 1
      lb2 = 1
      lb3 = 1
      lb4 = 1
      lb5 = 1
      lb6 = 1
      if( present(lbound1) ) lb1 = lbound1
      if( present(lbound2) ) lb2 = lbound2
      if( present(lbound3) ) lb3 = lbound3
      if( present(lbound4) ) lb4 = lbound4
      if( present(lbound5) ) lb5 = lbound5
      if( present(lbound6) ) lb6 = lbound6
      if(status == 0) allocate( strArray(lb1:dim1, lb2:dim2, &
           lb3:dim3, lb4:dim4, lb5:dim5, lb6:dim6), stat=status ) 
      if(status /= 0 ) goto 810
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
810   call err_handle(err, 200, &
              called_from = subname//" in MOD "//modname)
      return
   end function str6_alloc
   !##################################################################
   integer(kind=our_int) function int1_dealloc(intArray, err) &
        result(answer)
      ! Deallocates an integer array of rank 1
      implicit none
      ! declare arguments
      integer(kind=our_int), pointer :: intArray(:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "int1_dealloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(intArray) ) deallocate(intArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
   end function int1_dealloc
   !##################################################################
   integer(kind=our_int) function int2_dealloc(intArray, err) &
        result(answer)
      ! Deallocates an integer array of rank 2
      implicit none
      ! declare arguments
      integer(kind=our_int), pointer :: intArray(:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "int2_dealloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(intArray) ) deallocate(intArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
   end function int2_dealloc
   !##################################################################
   integer(kind=our_int) function int3_dealloc(intArray, err) &
        result(answer)
      ! Deallocates an integer array of rank 3
      implicit none
      ! declare arguments
      integer(kind=our_int), pointer :: intArray(:,:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "int3_dealloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(intArray) ) deallocate(intArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
   end function int3_dealloc
   !##################################################################
   integer(kind=our_int) function int4_dealloc(intArray, err) &
        result(answer)
      ! Deallocates an integer array of rank 4
      implicit none
      ! declare arguments
      integer(kind=our_int), pointer :: intArray(:,:,:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "int4_dealloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(intArray) ) deallocate(intArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
   end function int4_dealloc
   !##################################################################
   integer(kind=our_int) function int5_dealloc(intArray, err) &
        result(answer)
      ! Deallocates an integer array of rank 5
      implicit none
      ! declare arguments
      integer(kind=our_int), pointer :: intArray(:,:,:,:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "int5_dealloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(intArray) ) deallocate(intArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
   end function int5_dealloc
   !##################################################################
   integer(kind=our_int) function int6_dealloc(intArray, err) &
        result(answer)
      ! Deallocates an integer array of rank 6
      implicit none
      ! declare arguments
      integer(kind=our_int), pointer :: intArray(:,:,:,:,:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "int6_dealloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(intArray) ) deallocate(intArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
   end function int6_dealloc
   !##################################################################
   integer(kind=our_int) function dbl1_dealloc(dblArray, err) &
        result(answer)
      ! Deallocates a double-precision real array of rank 1
      implicit none
      ! declare arguments
      real(kind=our_dble), pointer :: dblArray(:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "dbl1_dealloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(dblArray) ) deallocate(dblArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
   end function dbl1_dealloc
   !##################################################################
   integer(kind=our_int) function dbl2_dealloc(dblArray, err) &
        result(answer)
      ! Deallocates a double-precision real array of rank 2
      implicit none
      ! declare arguments
      real(kind=our_dble), pointer :: dblArray(:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "dbl2_dealloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(dblArray) ) deallocate(dblArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
   end function dbl2_dealloc
   !##################################################################
   integer(kind=our_int) function dbl3_dealloc(dblArray, err) &
        result(answer)
      ! Deallocates a double-precision real array of rank 3
      implicit none
      ! declare arguments
      real(kind=our_dble), pointer :: dblArray(:,:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "dbl3dealloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(dblArray) ) deallocate(dblArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
   end function dbl3_dealloc
   !##################################################################
   integer(kind=our_int) function dbl4_dealloc(dblArray, err) &
        result(answer)
      ! Deallocates a double-precision real array of rank 4
      implicit none
      ! declare arguments
      real(kind=our_dble), pointer :: dblArray(:,:,:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "dbl4dealloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(dblArray) ) deallocate(dblArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
   end function dbl4_dealloc
   !##################################################################
   integer(kind=our_int) function dbl5_dealloc(dblArray, err) &
        result(answer)
      ! Deallocates a double-precision real array of rank 5
      implicit none
      ! declare arguments
      real(kind=our_dble), pointer :: dblArray(:,:,:,:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "dbl5dealloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(dblArray) ) deallocate(dblArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
   end function dbl5_dealloc
   !##################################################################
   integer(kind=our_int) function dbl6_dealloc(dblArray, err) &
        result(answer)
      ! Deallocates a double-precision real array of rank 5
      implicit none
      ! declare arguments
      real(kind=our_dble), pointer :: dblArray(:,:,:,:,:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "dbl6dealloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(dblArray) ) deallocate(dblArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
   end function dbl6_dealloc
   !##################################################################
   integer(kind=our_int) function str1_dealloc(strArray, err) &
        result(answer)
      ! Deallocates a character-string array of rank 1
      implicit none
      ! declare arguments
      character(len=*), pointer :: strArray(:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "str1_dealloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(strArray) ) deallocate(strArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
   end function str1_dealloc
   !##################################################################
   integer(kind=our_int) function str2_dealloc(strArray, err) &
        result(answer)
      ! Deallocates a character-string array of rank 2
      implicit none
      ! declare arguments
      character(len=*), pointer :: strArray(:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "str2_dealloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(strArray) ) deallocate(strArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
   end function str2_dealloc
   !##################################################################
   integer(kind=our_int) function str3_dealloc(strArray, err) &
        result(answer)
      ! Deallocates a character-string array of rank 3
      implicit none
      ! declare arguments
      character(len=*), pointer :: strArray(:,:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "str3_dealloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(strArray) ) deallocate(strArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
   end function str3_dealloc
   !##################################################################
   integer(kind=our_int) function str4_dealloc(strArray, err) &
        result(answer)
      ! Deallocates a character-string array of rank 4
      implicit none
      ! declare arguments
      character(len=*), pointer :: strArray(:,:,:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "str4_dealloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(strArray) ) deallocate(strArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
   end function str4_dealloc
   !##################################################################
   integer(kind=our_int) function str5_dealloc(strArray, err) &
        result(answer)
      ! Deallocates a character-string array of rank 5
      implicit none
      ! declare arguments
      character(len=*), pointer :: strArray(:,:,:,:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "str5_dealloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(strArray) ) deallocate(strArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
   end function str5_dealloc
   !##################################################################
   integer(kind=our_int) function str6_dealloc(strArray, err) &
        result(answer)
      ! Deallocates a character-string array of rank 6
      implicit none
      ! declare arguments
      character(len=*), pointer :: strArray(:,:,:,:,:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer :: status
      character(len=*), parameter :: subname = "str6_dealloc"
      ! begin
      answer = RETURN_FAIL
	  status = 0
      if( associated(strArray) ) deallocate(strArray, stat=status)
      if( status /= 0 ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 201, &
              called_from = subname//" in MOD "//modname)
      return
   end function str6_dealloc
   !##################################################################
end module dynamic_allocation
!#####################################################################

!#####################################################################
module weighted_least_squares
   ! This module contains computational procedures for weighted
   ! least squares regression.
   use program_constants
   use error_handler
   use matrix_methods
   implicit none
   private ! by default
   ! declare public procedures
   public :: fit_wls
   ! parameters private to this module
   character(len=*), parameter :: modname = "weighted_least_squares"
   !##################################################################
contains
   !##################################################################
   subroutine fit_wls(x, y, w, beta, cov_unscaled, scale, err)
      ! Regresses y on x, using weights in w.
      !    beta = estimated coefficients
      !    cov_unscaled = inverse of (X^T W X) matrix
      !    scale = unbiased estimate of scale parameter sigma^2
      implicit none
      ! declare argumants
      real(kind=our_dble), intent(in) :: x(:,:), y(:), w(:)
      real(kind=our_dble), intent(out) :: beta(:), &
           cov_unscaled(:,:), scale
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      character(len=*), parameter :: subname = "fit_wls"
      real(kind=our_dble), allocatable :: xtwx(:,:), xtwy(:)
      integer(kind=our_int) :: n, p, i, j, k, status
      ! check arguments
      n = size(x,1)
      p = size(x,2)

      if( ( size(y) /= n ) .or. ( size(beta) /= p ) .or. &
           ( size(w) /= n) .or. ( size(cov_unscaled,1) /= p ) .or. &
           ( size(cov_unscaled,2) /= p ) ) goto 700
      ! form xtwy and lower triangle of xtwx
      allocate( xtwx(p,p), xtwy(p), stat=status)
      if( status /= 0 ) goto 710
      do j = 1, p
         xtwy(j) = sum( x(:,j) * w(:) * y(:) )

         do k = 1, j
            xtwx(j,k) = sum( x(:,j) * w(:) * x(:,k) )
         end do
      end do


      ! compute cov_unscaled and beta
      call cholesky_saxpy( xtwx, err )
      call invert_lower( xtwx, err )
      call premult_lower_by_transpose( xtwx, cov_unscaled, err )
      beta = matmul( cov_unscaled, xtwy )
      ! compute scale
      scale = 0.D0
      do i = 1, n
         scale = scale + w(i) * ( y(i) - sum( x(i,:) * beta ) )**2
      end do
      scale = scale / real( n - p, kind=our_dble )
      ! deallocate workspaces
      deallocate(xtwx, xtwy, stat=status)
      if( status /= 0 ) goto 720
      ! normal exit
      return
      ! error traps
700   call err_handle(err, 301, &
           called_from = subname//" in MOD "//modname)
      return
710   call err_handle(err, 200, &
           called_from = subname//" in MOD "//modname)
      return
720   call err_handle(err, 201, &
           called_from = subname//" in MOD "//modname)
      return
800   call err_handle(err, 1000, &
           called_from = subname//" in MOD "//modname)
      deallocate(xtwx, xtwy)
      return
   end subroutine fit_wls
   !##################################################################
end module weighted_least_squares
!#####################################################################
subroutine fls(ncase, npred, x, y, beta, se_beta)
      use error_handler
      use program_constants
      use weighted_least_squares
      implicit none
      ! declare arguments
	  integer(kind=our_int) :: ncase
	  integer(kind=our_int) :: npred
      real(kind=our_dble) :: x(ncase,npred)
      real(kind=our_dble) :: y(ncase)
      real(kind=our_dble) :: beta(npred)
	  real(kind=our_dble) :: se_beta(npred)
      ! declare local variables and parameters
      type(error_type) :: err
      character(len=*), parameter :: modname = "modelfit"
      character(len=*), parameter :: subname = "modelfit"
	  real(kind=our_dble)  :: n(ncase)
      real(kind=our_dble), allocatable :: w(:), z(:), cov_beta(:,:)
      integer(kind=our_int) :: ijunk, status, j, i, iter
      real(kind=our_dble) :: scale
      !  check arguments
      n = 1.D0
      allocate(w(ncase), z(ncase), cov_beta(npred,npred), stat=status)
      beta = 0.D0
      w = 1.0
      z = y
      call fit_wls( x, z, w, beta, cov_beta, scale, err)
      ! store the covariance matrix in results
      do i=1,npred
		se_beta(i) = sqrt(cov_beta(i,i))*sqrt(scale)
	  enddo
      if( allocated( w ) ) deallocate( w )
      if( allocated( z ) ) deallocate( z )
      if( allocated( cov_beta ) ) deallocate( cov_beta )
end subroutine fls
!#####################################################################
subroutine flogit(ncase, npred, x, y, beta, se_beta)
      use error_handler
      use program_constants
      use weighted_least_squares
      implicit none
      ! declare arguments
	  integer(kind=our_int) :: ncase
	  integer(kind=our_int) :: npred
      real(kind=our_dble) :: x(ncase,npred)
      real(kind=our_dble) :: y(ncase)
      real(kind=our_dble) :: beta(npred)
	  real(kind=our_dble) :: se_beta(npred)
      integer(kind=our_int) :: maxits
      real(kind=our_dble) :: eps
      ! declare local variables and parameters
      type(error_type) :: err, warn
      character(len=*), parameter :: modname = "modelfit"
      character(len=*), parameter :: subname = "modelfit"
	  real(kind=our_dble)  :: n(ncase)
      real(kind=our_dble), allocatable :: w(:), &
            pi(:), z(:), beta_old(:), log_odds(:), &
            odds(:), cov_beta(:,:)
      integer(kind=our_int) :: ijunk, status, j, i, iter, k
      real(kind=our_dble), parameter :: pmin = 1.0D-6
	  real(kind=our_dble) :: fmax
	  real(kind=our_dble) :: fmin
      real(kind=our_dble) :: scale
	  logical :: converged
      !  check arguments
      ! allocate local arrays
      fmax=log(1.0/pmin-1.0)
      fmin=-fmax
      maxits = 100
      eps = 1D-05
      n = 1.0
      allocate(w(ncase), pi(ncase), &
           z(ncase), beta_old(npred), log_odds(ncase), &
           odds(ncase), cov_beta(npred,npred), stat=status)
      beta = 0.0
      iter = 0
      converged = .false.
      do
         iter = iter + 1
         beta_old = beta
         log_odds = matmul(x, beta)
         do k = 1,ncase
            if(log_odds(k) > fmax) then
	           pi(k) = 1.0  ! to prevent overflow
	        elseif(log_odds(k) < fmin) then
		       pi(k) = 0.0  ! to prevent overflow
            else
               odds(k) = exp(log_odds(k))
               pi(k) = odds(k) / (1.0 + odds(k))
            endif
         enddo
         w = n * pi * ( 1.0 - pi )
         z = log_odds + ( y - n*pi ) / w
         call fit_wls( x, z, w, beta, cov_beta, scale, err)
         if( all( abs( beta - beta_old ) <= eps*abs(beta_old) )) &
              converged = .true.
         if( converged .or. (iter >= maxits) ) exit
      end do
      ! store the covariance matrix in results
      do i=1,npred
		se_beta(i) = sqrt(cov_beta(i,i))
	  enddo
      if( .not.converged ) then 
	    call intpr("Did not converge by iterations", -1, iter,1)
	  endif
      ! normal exit
      if( allocated( w ) ) deallocate( w )
      if( allocated( pi ) ) deallocate( pi )
      if( allocated( z ) ) deallocate( z )
      if( allocated( beta_old ) ) deallocate( beta_old )
      if( allocated( log_odds ) ) deallocate( log_odds )
      if( allocated( odds ) ) deallocate( odds )
      if( allocated( cov_beta ) ) deallocate( cov_beta )
end subroutine flogit
!#####################################################################
! --------------------------------------------------
SUBROUTINE lsr(no,x1,x2,y,thr,res)       
! --------------------------------------------------
use program_constants
IMPLICIT NONE
! - - - arg types - - -
INTEGER :: no
INTEGER :: i
real(kind=our_dble) :: y(no)
real(kind=our_dble) ::  x1(no)
real(kind=our_dble) ::  x2(no)
real(kind=our_dble) ::  res(4,3)
real(kind=our_dble) ::  thr
real(kind=our_dble) ::  z
real(kind=our_dble) ::  p
real ( kind = 8 ) ::  alnorm
! - - - local declarations - - -
real(kind=our_dble) , allocatable::  dmat(:,:)
real(kind=our_dble) , allocatable::  beta(:)
real(kind=our_dble) , allocatable::  se_beta(:)
real(kind=our_dble) , allocatable::  pvalue(:)
INTEGER:: status
LOGICAL :: flag
! - - - begin - - -
	allocate(dmat(no,4), beta(1:4), se_beta(1:4), pvalue(1:4), stat=status)
	dmat=0
	dmat(:,1) = 1.0
	dmat(:,2) = x1			  
	dmat(:,3) = x2
    dmat(:,4) = x1*x2
    CALL fls(no, 4, dmat, y, beta, se_beta)
    do i = 1, 4
    	z = beta(i)/se_beta(i)
	    p = 2.0 * alnorm ( abs(z), .true. )
		if( ((1.0+p) /= p) .and. (p <= thr)) then
			pvalue(i) = p
		else
			pvalue(i) = 1.0
		endif
	enddo
	res(:,1) = beta
	res(:,2) = se_beta
	res(:,3) = pvalue
    if( allocated( dmat ) ) deallocate( dmat )
    if( allocated( beta ) ) deallocate( beta )
    if( allocated( se_beta ) ) deallocate( se_beta )         
    if( allocated( pvalue ) ) deallocate( pvalue )        
END SUBROUTINE lsr

! --------------------------------------------------
SUBROUTINE logr(no,x1,x2,y,thr,res)       
! --------------------------------------------------
use program_constants
IMPLICIT NONE
! - - - arg types - - -
INTEGER :: no
INTEGER :: i
real(kind=our_dble) :: y(no)
real(kind=our_dble) ::  x1(no)
real(kind=our_dble) ::  x2(no)
real(kind=our_dble) ::  res(4,3)
real(kind=our_dble) ::  thr
real(kind=our_dble) ::  z
real(kind=our_dble) ::  p
real ( kind = 8 ) ::  alnorm
! - - - local declarations - - -
real(kind=our_dble) , allocatable::  dmat(:,:)
real(kind=our_dble) , allocatable::  beta(:)
real(kind=our_dble) , allocatable::  se_beta(:)
real(kind=our_dble) , allocatable::  pvalue(:)
INTEGER:: status
LOGICAL :: flag
! - - - begin - - -
	allocate(dmat(no,4), beta(1:4), se_beta(1:4), pvalue(1:4), stat=status)
	dmat=0
	dmat(:,1) = 1.0
	dmat(:,2) = x1			  
	dmat(:,3) = x2
    dmat(:,4) = x1*x2
    CALL flogit(no, 4, dmat, y, beta, se_beta)
    do i = 1, 4
    	z = beta(i)/se_beta(i)
	    p = 2.0 * alnorm ( abs(z), .true. )
		if( ((1.0+p) /= p) .and. (p <= thr)) then
			pvalue(i) = p
		else
			pvalue(i) = 1.0
		endif
    enddo
	res(:,1) = beta
	res(:,2) = se_beta
	res(:,3) = pvalue
    if( allocated( dmat ) ) deallocate( dmat )
    if( allocated( beta ) ) deallocate( beta )
    if( allocated( se_beta ) ) deallocate( se_beta )         
    if( allocated( pvalue ) ) deallocate( pvalue )      
END SUBROUTINE logr
	
	