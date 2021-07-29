program visual_tests

    use iso_fortran_env
    use sobseq

    integer, parameter :: n_dim = 15
    integer, parameter :: n_samples = 10

    type(multi_dim_sobol_state) :: md_rng
    real(real64), dimension(:,:), allocatable :: seq

    allocate(seq(n_samples, n_dim))

    md_rng = multi_dim_sobol_state(n_dim, "data/clean_new-joe-kuo-6.21201")
    call md_rng%populate(seq)

    print '(20F8.4)', transpose(seq)

end program visual_tests
