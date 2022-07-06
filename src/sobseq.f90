!> Module containing a state type for implementing parallel
!> Sobol series generators with strided and skip-ahead generation.
!> Note that this uses the gray code implementation, so
!> Generated numbers are shuffled compared to the original series.
module sobseq

    use iso_fortran_env

    implicit none

    private

    integer, dimension(50), parameter :: s_default = &
        [1, 2, 3, 3, 4, 4, 5, 5, 5, 5, &
         5, 5, 6, 6, 6, 6, 6, 6, 7, 7, &
         7, 7, 7, 7, 7, 7, 7, 7, 7, 7, &
         7, 7, 7, 7, 7, 7, 8, 8, 8, 8, &
         8, 8, 8, 8, 8, 8, 8, 8, 8, 8  ]

    integer, dimension(50), parameter :: a_default = &
        [0,  1,  1,  2,  1,  4,  2,  4,  7,  11,  &
         13, 14, 1,  13, 16, 19, 22, 25, 1,  4,   &
         7,  8,  14, 19, 21, 28, 31, 32, 37, 41,  &
         42, 50, 55, 56, 59, 62, 14, 21, 22, 38,  &
         47, 49, 50, 52, 56, 67, 70, 84, 97, 103  ]

    integer, dimension(50,8), parameter :: m_default = &
        transpose(reshape( &
            [1, 0, 0, 0,  0,  0,  0,   0,   &
             1, 3, 0, 0,  0,  0,  0,   0,   &
             1, 3, 1, 0,  0,  0,  0,   0,   &
             1, 1, 1, 0,  0,  0,  0,   0,   &
             1, 1, 3, 3,  0,  0,  0,   0,   &
             1, 3, 5, 13, 0,  0,  0,   0,   &
             1, 1, 5, 5,  17, 0,  0,   0,   &
             1, 1, 5, 5,  5,  0,  0,   0,   &
             1, 1, 7, 11, 19, 0,  0,   0,   &
             1, 1, 5, 1,  1,  0,  0,   0,   &
             1, 1, 1, 3,  11, 0,  0,   0,   &
             1, 3, 5, 5,  31, 0,  0,   0,   &
             1, 3, 3, 9,  7,  49, 0,   0,   &
             1, 1, 1, 15, 21, 21, 0,   0,   &
             1, 3, 1, 13, 27, 49, 0,   0,   &
             1, 1, 1, 15, 7,  5,  0,   0,   &
             1, 3, 1, 15, 13, 25, 0,   0,   &
             1, 1, 5, 5,  19, 61, 0,   0,   &
             1, 3, 7, 11, 23, 15, 103, 0,   &
             1, 3, 7, 13, 13, 15, 69,  0,   &
             1, 1, 3, 13, 7,  35, 63,  0,   &
             1, 3, 5, 9,  1,  25, 53,  0,   &
             1, 3, 1, 13, 9,  35, 107, 0,   &
             1, 3, 1, 5,  27, 61, 31,  0,   &
             1, 1, 5, 11, 19, 41, 61,  0,   &
             1, 3, 5, 3,  3,  13, 69,  0,   &
             1, 1, 7, 13, 1,  19, 1,   0,   &
             1, 3, 7, 5,  13, 19, 59,  0,   &
             1, 1, 3, 9,  25, 29, 41,  0,   &
             1, 3, 5, 13, 23, 1,  55,  0,   &
             1, 3, 7, 3,  13, 59, 17,  0,   &
             1, 3, 1, 3,  5,  53, 69,  0,   &
             1, 1, 5, 5,  23, 33, 13,  0,   &
             1, 1, 7, 7,  1,  61, 123, 0,   &
             1, 1, 7, 9,  13, 61, 49,  0,   &
             1, 3, 3, 5,  3,  55, 33,  0,   &
             1, 3, 1, 15, 31, 13, 49,  245, &
             1, 3, 5, 15, 31, 59, 63,  97,  &
             1, 3, 1, 11, 11, 11, 77,  249, &
             1, 3, 1, 11, 27, 43, 71,  9,   &
             1, 1, 7, 15, 21, 11, 81,  45,  &
             1, 3, 7, 3,  25, 31, 65,  79,  &
             1, 3, 1, 1,  19, 11, 3,   205, &
             1, 1, 5, 9,  19, 21, 29,  157, &
             1, 3, 7, 11, 1,  33, 89,  185, &
             1, 3, 3, 3,  15, 9,  79,  71,  &
             1, 3, 7, 11, 15, 39, 119, 27,  &
             1, 1, 3, 1,  11, 31, 97,  225, &
             1, 1, 1, 3,  23, 43, 57,  177, &
             1, 3, 7, 7,  17, 17, 37,  71   ], [8, 50]))

    integer, parameter :: wp = real64
    integer, parameter :: N_M = 31 ! Generate at most 2^31 points
    ! Problems with sign bit when generating N_M = 32

    !> Type containing the state of a sobol sequence
    type, public :: sobol_state
      private
        integer :: v(N_M)   !< Direction numbers
        integer :: i = 1    !< Current number
        integer :: x = 0    !< Current value
        integer :: stride=0 !< Skip 2^this many values when generating
    contains
        procedure, public :: skip_ahead => sd_skip_ahead  !< Skip ahead to a specific position and return this value
        procedure, public :: next => sd_next        !< Generate the next value in the sequence
        procedure, public :: next_strided => sd_next_strided!< Generate the next value in the sequence (strided version)
        procedure, public :: populate => sd_populate    !< Populate an array of any size with random sobol numbers.
    end type sobol_state

    interface sobol_state
        module procedure :: sd_initialize !< Initialize direction numbers
    end interface sobol_state

    type, public :: multi_dim_sobol_state
        private
        type(sobol_state), allocatable, dimension(:) :: states
        integer :: n_dim
    contains
        procedure, public :: skip_ahead => md_skip_ahead
        procedure, public :: next => md_next
        procedure, public :: next_strided => md_next_strided
        procedure, public :: populate => md_populate
    end type multi_dim_sobol_state

    interface multi_dim_sobol_state
        module procedure :: md_initialize
        module procedure :: md_initialize_default
        module procedure :: md_initialize_from_file
    end interface multi_dim_sobol_state

contains

! Single Dimension

!> Initialize the direction numbers using a primitive polynomial
function sd_initialize(s, a, m_in, stride) result(new_ss)
    integer, intent(in) :: s !< Number of direction numbers / Mathematical polynomial basis of degree s
    integer, intent(in) :: a !< Coefficients of primitive polynomial
    integer, intent(in), dimension(s) :: m_in !< First direction numbers
    integer, intent(in), optional :: stride
    type(sobol_state) :: new_ss

    integer, dimension(N_M) :: m
    integer :: k, i, tmp

    if (size(m_in) /= s) error stop "m_in must be of size s"

    m(1:s) = m_in

    do k=s+1, N_M
      tmp=ieor(2**s * m(k-s), m(k-s))
      do i=1,s-1
        tmp = ieor(m(k-i) * 2**i * ai(a, s-i), &
                   tmp)
      end do
      m(k) = tmp
    end do

    do k=1, N_M
      new_ss%v(k) = (2**(N_M-k))*m(k)
    end do

    new_ss%x = 0
    new_ss%i = 0
    new_ss%stride = 0
    if (present(stride)) new_ss%stride = stride

end function sd_initialize


!> Generate a value at a specific position i
function sd_skip_ahead(self, i) result(output)
  class (sobol_state), intent(inout) :: self
  integer, intent(in) :: i

  real(kind=wp) :: output
  integer :: g ! Gray code representation of i
  integer :: j, tmp

  g = ieor(i,i/2)
  self%x = 0
  self%i = i

  tmp = ai(g,1) * self%v(1)
  do j=2, N_M
    tmp = ieor(tmp,ai(g,j) * self%v(j))
  end do
  output = real(tmp, kind=wp) * 2.0_wp**(-N_M)
  self%x = tmp

end function sd_skip_ahead


!> Generate the next value in a series
!> And update the self function
function sd_next(self) result(next_elem)

  class (sobol_state), intent(inout) :: self

  real(kind=wp) :: next_elem

  self%x = ieor(self%x, self%v(i4_bit_lo0(self%i)))
  self%i = self%i + 1
  next_elem = real(self%x, kind=wp) * 2.0_wp**(-N_M)

end function sd_next

!> Generate the next value in a series
!> And update the self function
function sd_next_strided(self) result(next_elem)
  class (sobol_state), intent(inout) :: self

  real(kind=wp) :: next_elem

  self%x = ieor(self%x, ieor(self%v(self%stride), self%v(&
            i4_bit_lo0(ior(self%i, 2**self%stride - 1)))))
  self%i = self%i + 2**self%stride
  next_elem = real(self%x, kind=wp) * 2.0_wp**(-N_M)
end function sd_next_strided

subroutine sd_populate(self, arr)
    class (sobol_state), intent(inout) :: self
    real(kind=wp), dimension(:), intent(out) :: arr

    integer :: i

    do i=1, size(arr)
        arr(i) = self%next()
    end do
end subroutine sd_populate

! Multiple Dimension

!> Initialize the direction numbers using a primitive polynomial
function md_initialize(n_dim, s, a, m_in, stride) result(new_mdss)
    integer, intent(in) :: n_dim
    integer, dimension(:), intent(in) :: s !< Number of direction numbers / Mathematical polynomial basis of degree s
    integer, dimension(:), intent(in) :: a !< Coefficients of primitive polynomial
    integer, dimension(:,:), intent(in) :: m_in !< First direction numbers
    integer, intent(in), optional :: stride
    type(multi_dim_sobol_state) :: new_mdss

    integer, allocatable, dimension(:) :: m_row
    integer :: i

    if ((size(s) < n_dim).or.(size(a) < n_dim).or.(size(m_in, 1) < n_dim)) error stop

    new_mdss%n_dim = n_dim

    allocate(new_mdss%states(n_dim))
    allocate(m_row(s(size(s))))

    do i=1, n_dim
        ! this fixes runtime warning (temprary array)
        m_row = m_in(i,1:s(i))
        new_mdss%states(i) = sobol_state(s(i), a(i), m_row, stride)
    end do

end function md_initialize

function md_initialize_default(n_dim, stride) result(new_mdss)
    integer, intent(in) :: n_dim
    integer, optional, intent(in) :: stride
    type(multi_dim_sobol_state) :: new_mdss

    integer, dimension(n_dim) :: s
    integer, dimension(n_dim) :: a
    integer, allocatable, dimension(:,:) :: m

    integer :: max_cols

    max_cols = s_default(n_dim)
    allocate(m(n_dim, max_cols))
    s = s_default(1:n_dim)
    a = a_default(1:n_dim)
    m = m_default(1:n_dim,:max_cols)

    new_mdss = md_initialize(n_dim, s, a, m, stride)

end function md_initialize_default

function md_initialize_from_file(n_dim, file_name, stride) result(new_mdss)
    integer, intent(in) :: n_dim
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: stride
    type(multi_dim_sobol_state) :: new_mdss

    integer, dimension(n_dim) :: s
    integer, dimension(n_dim) :: a
    integer, allocatable, dimension(:,:) :: m

    integer :: err_stat, i, d
    integer, dimension(21) :: line
    character(len=80) :: err_msg

    open(unit=8, file=file_name, status="old", form="formatted", &
         iostat=err_stat, iomsg=err_msg)
    if (err_stat /= 0) then
        write (*,'(A,I6)') "Error opening file, stat: ", err_stat
        write (*,'(A)') trim(err_msg)
        error stop
    end if

    open(unit=9, file=file_name, status="old", form="formatted", &
         iostat=err_stat, iomsg=err_msg)
    if (err_stat /= 0) then
        write (*,'(A,I6)') "Error opening file, stat: ", err_stat
        write (*,'(A)') trim(err_msg)
        error stop
    end if

    do i=1,n_dim
        read (8,*) line
    end do
    allocate(m(n_dim, line(2)))
    do i=1,n_dim
        read (9,*) d, s(i), a(i), m(i,:)
    end do

    close(unit=8, iostat=err_stat, iomsg=err_msg)
    if (err_stat /= 0) then
        write (*,'(A,I6)') "Error closing file, stat: ", err_stat
        write (*,'(A)') trim(err_msg)
        error stop
    end if

    close(unit=9, iostat=err_stat, iomsg=err_msg)
    if (err_stat /= 0) then
        write (*,'(A,I6)') "Error closing file, stat: ", err_stat
        write (*,'(A)') trim(err_msg)
        error stop
    end if

    new_mdss = md_initialize(n_dim, s, a, m, stride)

end function md_initialize_from_file

!> Generate a value at a specific position i
function md_skip_ahead(self, i) result(output)
  class (multi_dim_sobol_state), intent(inout) :: self
  integer, intent(in) :: i

  real(kind=wp), dimension(self%n_dim) :: output

  integer :: j

  do j=1, self%n_dim
      output(j) = self%states(j)%skip_ahead(i)
  end do

end function md_skip_ahead


!> Generate the next value in a series
!> And update the state function
function md_next(self) result(next_elem)
  class (multi_dim_sobol_state), intent(inout) :: self

  real(kind=wp), dimension(self%n_dim) :: next_elem

  integer :: i

  do i=1, self%n_dim
      next_elem(i) = self%states(i)%next()
  end do
end function md_next

!> Generate the next value in a series
!> And update the state function
function md_next_strided(self) result(next_elem)
  class (multi_dim_sobol_state), intent(inout) :: self

  real(kind=wp), dimension(self%n_dim) :: next_elem

  integer :: i

  do i=1, self%n_dim
      next_elem(i) = self%states(i)%next_strided()
  end do
end function md_next_strided

subroutine md_populate(self, mat)
    class (multi_dim_sobol_state), intent(inout) :: self
    real(kind=wp), dimension(:,:), intent(out) :: mat

    integer :: i

    if (size(mat, 2) /= self%n_dim) error stop

    do i=1, size(mat, 1)
        mat(i, :) = self%next()
    end do
end subroutine md_populate

! Private Functions

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

!> Return the position of the lowest (rightmost) 0 bit in a int(4)
function i4_bit_lo0(num)
  implicit none

  integer(kind=4), intent(in) :: num
  integer :: i4_bit_lo0

  do i4_bit_lo0=1,bit_size(num)
    if (.not. btest(num,i4_bit_lo0-1)) return
  end do
end function i4_bit_lo0

end module sobseq
