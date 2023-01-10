program barrier

    use OMP_LIB

    implicit none
    integer, parameter :: n = 10000
    integer :: i, tmp, result

    tmp = 0

!$omp parallel
    do i = 1, n
!$omp critical
        tmp = tmp + 1
!$omp end critical
    enddo
!$omp barrier     
!$omp master
    result = tmp  !wrong due to data race
    print*, 'result = ', result
!$omp end master

!$omp end parallel

end program barrier


