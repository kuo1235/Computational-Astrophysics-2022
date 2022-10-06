program pi

    print *, "Type the number of partition (N) you want to make"
    read *, N

    call compute_integral(N,area)

    print*, "PI=", 2. * area

end program


subroutine compute_integral(N,A)
    implicit none
    integer, intent(in)::N
    real, intent(out) :: A

    real :: x, dx, h, dA 
    real :: my_func
    integer :: i

    dx=2/real(N)
    A=0

    do i = 1, N

       x = -1 - 1/N + dx*real(i)
       h = my_func(x)
       dA= dx * h
       A = A + dA

    enddo

    return	

end subroutine compute_integral

real function my_func(x)
    real :: x
    my_func = sqrt(1.0-x**2)
    return
end function





