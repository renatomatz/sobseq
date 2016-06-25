!> Module containing a state type for implementing parallel
!> Sobol series generators with strided and skip-ahead generation.
!> Note that this uses the gray code implementation, so
!> Generated numbers are shuffled compared to the original series.
module mod_sob_seq
  implicit none

  integer, parameter :: N_M = 32 ! Generate at most 2^32 points

  !> Type containing the state of a sobol sequence
  type sobol_state
    private
      integer :: v(N_M)   !< Direction numbers
      integer :: i = 1    !< Current number
      integer :: x = 0    !< Current value
      integer :: stride=0 !< Skip 2^this many values when generating
  contains
      procedure, public :: initialize !< Initialize direction numbers
      procedure, public :: skip_ahead !< Skip ahead to a specific position and return this value
      procedure, public :: next       !< Generate the next value in the sequence
  end type sobol_state
contains


!> Initialize the direction numbers using a primitive polynomial
subroutine initialize(state, s, a, m_in)
  implicit none
  class (sobol_state), intent(inout) :: state
  integer, intent(in) :: s !< Number of direction numbers
  integer, intent(in) :: a !< Coefficients of primitive polynomial
  integer, intent(in), dimension(s) :: m_in !< First direction numbers

  integer, dimension(N_M) :: m
  integer :: k, i, tmp

  m(1:s) = m_in

  do k=s+1, N_M
    tmp=ieor(2**s * m(k-s), m(k-s))
    do i = k-s+1, k-1
      tmp = ieor(m(i) * 2**(k-i) * ai(a, k-i), &
                 tmp)
    end do
    m(k) = tmp
  end do

  do k=1, N_M
    state%v(k) = m(k)/2**k
  end do
end subroutine initialize


!> Generate a value at a specific position i
function skip_ahead(state, i) result(output)
  implicit none
  class (sobol_state), intent(inout) :: state
  integer, intent(in) :: i
  real    :: output
  integer :: g ! Gray code representation of i
  integer :: j

  g = ieor(i,i/2)
  state%x = 0
  state%i = i
  do j=1,N_M
    if (btest(g,j-1)) state%x = ieor(state%x,state%v(j))
  enddo
  output = real(state%x) * 2.0**(-32)
end function skip_ahead


!> Generate the next (strided) value in a series
!> And update the state function
function next(state)
  use mod_bit_hilo
  implicit none
  class (sobol_state), intent(inout) :: state
  real :: next

  ! TODO: check if the strided case reduces to the nonstrided case and remove this,
  ! or check the speed difference, and provide two different functions
  if (state%stride .eq. 0) then
    state%x = ieor(state%x, state%v(i4_bit_lo0(state%i)))
  else
    state%x = ieor(state%x, ieor(state%v(state%stride), state%v(&
        i4_bit_lo0(ior(state%i, 2**state%stride - 1)))))
  endif

  state%i = state%i + 2**state%stride

  next = real(state%x) * 2.0**(-32)
end function next


!> Returns the value of the bit at position i (1-based index)
function ai(a,i)
implicit none
integer, intent(in) :: a, i
integer :: ai

if (btest(a,i-1)) then
  ai = 1
else
  ai = 0
end if
end function ai
end module mod_sob_seq