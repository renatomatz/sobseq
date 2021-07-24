program basic_tests

    use iso_fortran_env
    use sobseq
    use cafut

    type(TestSuite) :: ts

    type(sobol_state) :: rng, rng2
    type(multi_dim_sobol_state) :: md_rng

    integer, parameter :: s=1, a=0, m(1) = (/1/)
    integer, parameter :: N_samples = 17

    real(real64), dimension(:), allocatable :: res, tgt
    real(real64), dimension(:,:), allocatable :: res2, tgt2

    integer :: i

    ts = TestSuite("Single Dimension")

    allocate(res(N_samples), tgt(N_samples))

    tgt = [0.500000000,     &
           0.250000000,     & 
           0.750000000,     & 
           0.375000000,     & 
           0.875000000,     & 
           0.125000000,     & 
           0.625000000,     & 
           0.312500000,     & 
           0.812500000,     & 
           6.25000000E-02,  &
           0.562500000,     &
           0.187500000,     &
           0.687500000,     &
           0.437500000,     &
           0.937500000,     &
           0.468750000,     &
           0.968750000]

    rng = sobol_state(s,a,m)
    res = [ (rng%skip_ahead(i), i=1, N_samples) ]
    call ts%add(TestRealArrVal("Skip-ahead"), res, tgt)

    rng = sobol_state(s,a,m)
    res = [ (rng%next(), i=1, N_samples) ]
    call ts%add(TestRealArrVal("Next"), res, tgt)

    deallocate(res, tgt)
    allocate(res(0:N_samples), tgt(0:N_samples))

    rng = sobol_state(s,a,m, stride=1)
    res(0) = rng%skip_ahead(0)

    rng2 = sobol_state(s,a,m, stride=1)
    res(1) = rng2%skip_ahead(1)

    do i=2,17
        if (mod(i,2) == 0) then
          res(i) = rng%next()
        else
          res(i) = rng2%next()
        endif
    end do

    tgt = [0.00000000 ,     &
           0.500000000,     & 
           0.500000000,     & 
           0.250000000,     & 
           0.250000000,     & 
           0.750000000,     & 
           0.750000000,     & 
           0.375000000,     & 
           0.375000000,     & 
           0.875000000,     &
           0.875000000,     &
           0.125000000,     &
           0.125000000,     &
           0.625000000,     &
           0.625000000,     &
           0.312500000,     &
           0.312500000,     &
           0.812500000]

    call ts%add(TestRealArrVal("Next Stride 1"), res, tgt)

    call ts%runTests()

    deallocate(res, tgt)

    ts = TestSuite("Multiple Dimensions")

    allocate(res(N_samples), tgt(N_samples))
    allocate(res2(N_samples, 2), tgt2(N_samples, 2))
        
    md_rng = multi_dim_sobol_state(2)
    
    call md_rng%populate(res2)

    tgt2(:,1) = [0.5000000000,     &
                 0.2500000000,     & 
                 0.7500000000,     & 
                 0.3750000000,     & 
                 0.8750000000,     & 
                 0.1250000000,     & 
                 0.6250000000,     & 
                 0.3125000000,     & 
                 0.8125000000,     & 
                 0.0625000000,     &
                 0.5625000000,     &
                 0.1875000000,     &
                 0.6875000000,     &
                 0.4375000000,     &
                 0.9375000000,     &
                 0.4687500000,     &
                 0.9687500000]

    tgt2(:,2) = [0.5000000000,     &
                 0.2500000000,     & 
                 0.7500000000,     & 
                 0.6250000000,     & 
                 0.1250000000,     & 
                 0.8750000000,     & 
                 0.3750000000,     & 
                 0.9375000000,     & 
                 0.4375000000,     & 
                 0.6875000000,     &
                 0.1875000000,     &
                 0.3125000000,     &
                 0.8125000000,     &
                 0.0625000000,     &
                 0.5625000000,     &
                 0.4687500000,     &
                 0.9687500000]

    res = res2(:,1)
    tgt = tgt2(:,1)

    call ts%add(TestRealArrVal("Populate Dim 1"), res, tgt) 

    res = res2(:,2)
    tgt = tgt2(:,2)

    call ts%add(TestRealArrVal("Populate Dim 2"), res, tgt) 

    deallocate(res, tgt, res2, tgt2)

    call ts%runTests()

end program basic_tests
