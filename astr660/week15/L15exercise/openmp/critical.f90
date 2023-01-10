program critical

    use OMP_LIB

    implicit none
    integer, parameter :: n = 100
    integer :: i, tmp

    tmp = 0

!$omp parallel do 
    do i = 1, n
!$omp critical
        tmp = tmp + 1  !wrong due to data race
!$omp end critical
    enddo
!$omp end parallel do

    print*, 'tmp = ', tmp

end program critical


