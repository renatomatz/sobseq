program stride_tests

    use iso_fortran_env
    use sobseq
    use ogpf

    integer, parameter :: N_samples = 1024

    type(multi_dim_sobol_state) :: md_rng1
    type(multi_dim_sobol_state) :: md_rng2
    type(gpf) :: gp

    real(real64), dimension(:,:), allocatable :: seq
    real(real64), dimension(:), allocatable :: trash

    integer :: i

    allocate(seq(N_samples, 4))

    md_rng1 = multi_dim_sobol_state(4, 3)
    md_rng2 = multi_dim_sobol_state(4, 3)
    trash = md_rng2%skip_ahead(1)

    call gp%title("Sobol Sequence Visualization")
    call gp%xlabel("Dimension 1")
    call gp%ylabel("Dimension 2")

    do i=2,N_samples,2
        seq(i-1,:) = md_rng1%next_strided()
        seq(i,:)   = md_rng2%next_strided()
    end do

    call gp%plot(seq(:,1), seq(:,2:), "with points")

    call execute_command_line("rm ogpf_temp_script.gp")

end program stride_tests
