program pi


    implicit none 
    
    real :: x, dx, h, dA, area 
    integer :: i, N

    print *, "Type the number of partition (N) you want to make"
    read *, N
    
    dx=2/real(N)
    area=0

    do i = 1, N

        x = -1 - 1/N + dx*real(i)
        h = sqrt(1.0 - x**2.0)
        dA= dx * h
        area =area + dA

    enddo	

    print*, "PI=", 2. * area

end program


