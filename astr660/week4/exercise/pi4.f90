program pi
    
    implicit none 
   
    integer :: i, N
    real, parameter :: test_pi = 4.0*atan(1.0)
    real :: error, area
    integer, parameter :: NMAX=7	
    integer, dimension(NMAX) :: n_iteration
      


    !---  set up N
    n_iteration = (/10,100,1000,10000,100000,1000000,10000000/)

    !--- open a file to store the error
    !--- open(u), where u is a valid unit number specifier
    open(unit=11, file="error.dat")
    
    !--- write the header
    !--- [unit=]: unit number; [FMT=] format-spec 
    write(unit=11, FMT=100) "# ", "N", "Error"

    !--- do the integral
    do i=1, NMAX
        N = n_iteration(i)
        call compute_integral(N,area)	
        print *, "N = ",n, " PI = ", 2.*area
        error = abs(2.*area-test_pi)/test_pi

        write(11,200) n_iteration(i), error

    enddo    
 
100 format(a2, a8, 2a24)
200 format(2x, i8, 2F24.14) 

    close(11)

end program pi


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





