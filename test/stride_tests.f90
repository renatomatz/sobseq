program stride_tests

    use iso_fortran_env
    use sobseq
    use ogpf

    integer, parameter :: N_sobol = 3
    integer, parameter :: N_samples = 1024

    type(multi_dim_sobol_state), dimension(N_sobol) :: md_rng
    type(gpf) :: gp

    real(real64), dimension(:,:), allocatable :: seq
    real(real64), dimension(:), allocatable :: trash

    integer :: temp, stride
    integer :: i, j
    logical :: is_pow_2

    temp = N_sobol
    stride = 0
    is_pow_2 = .true.
    do
        is_pow_2 = is_pow_2.and.(mod(temp,2)==0)
        temp = ishft(temp,-1)
        stride = stride + 1
        if (temp == 1) exit
    end do
    if (.not.is_pow_2) stride = stride + 1

    allocate(seq(N_samples, 4))

    do i=1,N_sobol
        md_rng(i) = multi_dim_sobol_state(4, stride)
        trash = md_rng(i)%skip_ahead(i-1)
    end do

    call gp%title("Sobol Sequence Visualization")
    call gp%xlabel("Dimension 1")
    call gp%ylabel("Dimension 2")

    do i=N_sobol,N_samples,N_sobol
        do j=1,N_sobol
            seq(i-(j-1),:) = md_rng(j)%next_strided()
        end do
    end do

    call gp%plot(seq(:,1), seq(:,2:), "with points")

    call execute_command_line("rm ogpf_temp_script.gp")

end program stride_tests
