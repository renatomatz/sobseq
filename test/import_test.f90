program test

    integer, parameter :: n_dim = 12

    integer, dimension(n_dim) :: s
    integer, dimension(n_dim) :: a
    integer, allocatable, dimension(:,:) :: m

    integer :: err_stat, i, d
    integer, dimension(21) :: line
    character(len=80) :: err_msg
    character(len=34) :: f

    f = "data/clean_new-joe-kuo-6.21201"

    open(unit=8, file=f, status="old", form="formatted", &
         iostat=err_stat, iomsg=err_msg)
    if (err_stat /= 0) then
        write (*,'(A,I6)') "Error opening file, stat: ", err_stat
        write (*,'(A)') trim(err_msg)
        error stop
    end if

    open(unit=9, file=f, status="old", form="formatted", &
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

    print *, "s:"
    print 100, s

    print *, "a:"
    print 100, a

    print *, "m:"
    print 200, transpose(m)

    100 format(12I3)
    200 format(5I3)

end program test
