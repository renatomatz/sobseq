program stride_tests

    use iso_fortran_env
    use sobseq
    use ogpf

    integer, parameter :: N_sobol = 3
    integer, parameter :: N_samples = 1024

    type(multi_dim_sobol_state), dimension(N_sobol) :: md_rng
    type(multi_dim_sobol_state) :: check_md_rng
    type(gpf) :: gp

    real(real64), dimension(N_samples,4) :: seq
    real(real64), dimension(N_samples,4) :: check_seq
    real(real64), dimension(:), allocatable :: trash

    integer :: temp, stride, extra
    integer :: i, j, k

    temp = N_sobol
    stride = 0
    do
        temp = ishft(temp,-1)
        if (temp == 0) exit
        stride = stride + 1
    end do
    extra = N_sobol - 2**stride
    ! print *, N_sobol, stride, extra

    do i=1,N_sobol
        md_rng(i) = multi_dim_sobol_state(4, stride)
        trash = md_rng(i)%skip_ahead(i-1)
    end do
    check_md_rng = multi_dim_sobol_state(4)

    call gp%title("Sobol Sequence Visualization")
    call gp%xlabel("Dimension 1")
    call gp%ylabel("Dimension 2")

    gen: do i=1,N_samples,N_sobol
        do j=1,N_sobol
            k = i+j-1
            if (k > N_samples) exit gen
            seq(k,:) = md_rng(j)%next_strided()
            do k=1,extra
                trash = md_rng(j)%next()
            end do
        end do
    end do gen

    ! for some reason the draws are one step behind
    ! trash = check_md_rng%next()
    call check_md_rng%populate(check_seq)

    ! print 10, seq(:,1)
    ! 10 format(1F12.6)
    ! print 20, transpose(reshape([seq(:,1), check_seq(:,1)], [N_samples,2]))
    ! 20 format(2F12.6)
    ! if (.not.all(abs(seq-check_seq)<1e-5)) &
    !     error stop "missmatch in strided and non-strided sequences"

    call gp%plot(seq(:,1), seq(:,2:), "with points")

    call execute_command_line("rm ogpf_temp_script.gp")

end program stride_tests
