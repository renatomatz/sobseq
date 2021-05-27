program visual_tests

    use iso_fortran_env
    use sobseq
    use ogpf

    integer, parameter :: N_samples = 1024

    type(multi_dim_sobol_state) :: md_rng
    type(gpf) :: gp

    real(real64), dimension(:,:), allocatable :: seq

    allocate(seq(N_samples, 4))
        
    md_rng = multi_dim_sobol_state(4)

    call gp%title("Sobol Sequence Visualization") 
    call gp%xlabel("Dimension 1")
    call gp%ylabel("Dimension 2")

    call md_rng%md_populate(seq)
    call gp%plot(seq(:,1), seq(:,2:), "with points")

    call execute_command_line("rm ogpf_temp_script.gp")

end program visual_tests
