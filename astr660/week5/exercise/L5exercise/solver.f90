module solver

    ! ---------------------------------------------------------
    ! ---------------------------------------------------------

    implicit none
    contains 

        !--------------------------------------------------
        !
        ! Implment the bisection method to solve the func
        !
        !
        ! Inputs:   func
        !
        ! Outputs:
        !           xs   : solution of func
        !           error: relative error
        !
        !--------------------------------------------------


        subroutine bisection(func, xs, err)
        implicit none
        real, external    :: func    ! the function to solve, tell the module that "func" will be declared in other files.
        real, intent(out) :: xs      ! solution
        real, intent(out) :: err     ! error
        real, save :: a = 1.0        ! bracking interval [a,b]
        real, save :: b = 3.0        ! bracking interval [a,b]
        real  :: fa, fx              ! f(a) and f(x)
           
        xs = (a+b)/2.        
                
        fa = func(a)
        fx = func(xs)

        if (sign(1.0, fa) .eq. sign(1.0, fx)) then
            a = xs
        else        
            b = xs
        end if  

        err = abs(fx) 

        end subroutine bisection


        !--------------------------------------------------
        !
        ! Implment the Newton's method to solve my_func
        !
        !
        ! Inputs: func, dfunc
        !
        ! Outputs:
        !           xs   : solution of my_func
        !           error: relative error
        !
        !--------------------------------------------------
        subroutine newton(func, dfunc, xs, err)
        implicit none
        real, external    :: func  ! the function to solve
        real, external    :: dfunc ! the first derivative of the function to solve
        real, intent(out) :: xs    ! solution
        real, intent(out) :: err   ! error
        real, save :: x = 2.0      ! trial value
        real :: fx
        real :: dfdx
       
        fx = func(x) 
        dfdx = dfunc(x) 


        if (sign(1.0,fx) .gt. 0.) then

        x = x - fx/dfdx
        xs = x

        end if
        
        err = abs(fx)

        end subroutine newton


        !--------------------------------------------------
        !
        ! Implment the Secant method to solve my_func
        !
        !
        ! Inputs: None
        !
        ! Outputs:
        !           xs   : solution of my_func
        !           error: relative error
        !
        !--------------------------------------------------
        subroutine secant(func, xs, err)
        implicit none
        real, external    :: func  ! the function to solve
        real, intent(out) :: xs    ! solution
        real, intent(out) :: err   ! error

        real,save :: x0 = 1.0  ! initial guess
        real,save :: x1 = 3.0  ! initial guess

        real :: fx0, fx1       ! f(x0) and f(x1)
        
        fx0 = func(x0)
        fx1 = func(x1)

        if (sign(1.0, fx1) .eq. sign(1.0, fx0)) then              
            x1 = x1 - fx1/((fx1-fx0)/(x1-x0))
            xs = x1
        else 
            x0 = x0 - fx0/((fx1-fx0)/(x1-x0))
            xs = x0
        end if

        err = abs(func(xs))
    


end subroutine secant


end module solver
