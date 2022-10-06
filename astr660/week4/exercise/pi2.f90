program pi


    implicit none 
    
    real :: x, dx, h, dA, area 
    real :: my_func
    integer :: i, N

    print *, "Type the number of partition (N) you want to make"
    read *, N

    dx=2/real(N)
    area=0

    do i = 1, N

        x = -1 - 1/N + dx*real(i)
        h = my_func(x)
        dA= dx * h
        area =area + dA

    enddo	

    print*, "PI=", 2. * area

end program

real function my_func(x)
    real :: x
    my_func = sqrt(1.0-x**2)
    return
end function





