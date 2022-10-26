!
! National Tsing Hua University
!
! ASTR 660 Computational Astrophysics
!
! Created:  Kuo-Chuan Pan 2020
! Modified: Karen Yang 2022.10.15
!
! Problem:
!
!        Solving non-linear equations
!
program linear
    use linalg                  !module defined in linalg.f90 
    implicit none
    integer, parameter  :: N = 3
    real,dimension(N,N) :: lower, upper, A, P, Ainv
    real,dimension(N) :: b
    real,dimension(N) :: x
    real,dimension(4,4) :: aa,ll,uu,pp
    integer :: i,j

    ! lower triangle
    lower(1,1) = -1.0
    lower(1,2) =  0.0
    lower(1,3) =  0.0

    lower(2,1) = -6.0
    lower(2,2) = -4.0
    lower(2,3) =  0.0

    lower(3,1) =  1.0
    lower(3,2) =  2.0
    lower(3,3) =  2.0
    
    ! upper triangle
    upper(1,1) =  1.0
    upper(1,2) =  2.0
    upper(1,3) =  2.0

    upper(2,1) =  0.0
    upper(2,2) = -4.0
    upper(2,3) = -6.0

    upper(3,1) =  0.0
    upper(3,2) =  0.0
    upper(3,3) = -1.0

    ! A matrix
    A(1,1) =  1.0
    A(1,2) =  2.0
    A(1,3) =  2.0

    A(2,1) =  4.0
    A(2,2) =  6.0
    A(2,3) =  8.0

    A(3,1) =  4.0
    A(3,2) =  8.0
    A(3,3) = 10.0
    ! the vectore b
    b(1) =  3.0
    b(2) =  -6.0
    b(3) =  1.0

    !call solve_lower_triangular_matrix(N,lower,b,x)
    !call solve_upper_triangular_matrix(N,upper,b,x)
  
    call LU_decomposition(N,A,lower,upper)
    call solve_lu(N,A,b,x)

    !call mat_print("L",lower)
    !call mat_print("U",upper)
    call mat_print("A", A)

    !print *, "vector   b = ",b
    !print *, "solution x = ",x
    
    call mat_print("L", lower)
    call mat_print("U", upper)

end program linear 


