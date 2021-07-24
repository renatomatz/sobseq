program timing

    use iso_fortran_env
    use sobseq

    type(sobol_state) :: rng

    real(real64), dimension(:), allocatable :: tmp

    integer, parameter :: s=1, a=0, m(1) = (/1/)
    integer, parameter :: N_samples = 1000000
    real(real64) :: t0, t1

    integer :: i

    allocate(tmp(N_samples))

    write(*,"(A,i10,A)", advance="no") "Testing time skip_ahead (", N_samples, "): "

    rng = sobol_state(s,a,m)

    call cpu_time(t0)
    do i=1,N_samples
        tmp(i) = rng%skip_ahead(i)
    end do
    call cpu_time(t1)

    write(*,*) (t1-t0)/N_samples, sum(tmp)/N_samples

    write(*,"(A,i10,A)",advance="no") "Testing time next (", N_samples, "): "

    rng = sobol_state(s,a,m)

    call cpu_time(t0)
    do i=1,N_samples
        tmp(i) = rng%next()
    end do
    call cpu_time(t1)

    write(*,*) (t1-t0)/N_samples, sum(tmp)/N_samples

end program timing
