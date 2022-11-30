!---------------------------------------------------
!
! National Tsing Hua University
!
! ASTR 660 Computational Astrophysics
!
! Created:  Kuo-Chuan Pan 2020
! Modified: Karen Yang 2022.10.25
!
! Problem:
!
!        Solving boundary value problems

program shooting1

    use Solver 
    implicit none

    external :: my_func

    real :: x   
    real :: v, vx, vy
    real :: dt,t,tend
    integer :: i,n

    real :: a, b

    real :: try_y2
    real, dimension(2) :: y, ynext, yend
    

    dt    = 0.01  ! step size
    t    = 0.0   ! initial t

    v = 30.0

    a = 0.0 ! guessing the lower limit
    b = 30.0 ! guessing the upper limit
    
    n = 2

    yend(1) = 1.0 ! only use to dodgeing the if and do statement

    do while ( abs(yend(1)) >= 1.0E-6 )
    
        try_y2 = (a+b)/2.0
        
        t = 0.0
        y(1) = 0.0
        y(2) = try_y2

        vx = (v**2.0 - try_y2**2.0)**0.5

        tend = 50.0/vx
  
        do while ( ABS(t-tend) >= 1.0E-4 )
            
            if ((t+dt) .ge. tend) then
                dt = tend - t
            endif

            !call rk2(n, y, ynext, t, h, my_func)
            call rk4(n, y, ynext, t, dt, my_func)
            
            y = ynext
            t = t + dt
            
            yend = ynext

        enddo

        print *, "y = ", ynext(1)

        if (yend(1) < 0.0) then
            
            a = (a+b)/2.0
             
        else  

            b = (a+b)/2.0

        endif

    enddo

    print *, "Vy = ", -ynext(2) ! Since Vy is negative at tend
    print *, "The desired value is about 8.52"


end program shooting1

!-----------------------------------------------
!
!  Solving Boundary Value problem
!
!  u'' = 6t      0 < t < 1
!
!  with BC
!
!         u(t=0) = 1 and
!         u(t=1) = 1
!
!-----------------------------------------------

subroutine my_func(n, t, yin, k)
    implicit none
    integer, intent(in) :: n  ! number of ODEs
    real, intent(in)    :: t
    real, dimension(n), intent(in)  :: yin
    real, dimension(n), intent(out) :: k    ! dydt

    ! in this example  n = 2
    k(1) = yin(2)    
    k(2) = -9.8 ! g=9.8   
    return
end subroutine my_func







