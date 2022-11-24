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
    real, dimension(2) :: y, ynext

    dt    = 0.01  ! step size
    t    = 0.0   ! initial t

    v = 30.0

    !tend = 50.0/vx   ! final t

    !vy = ?
    !vx = (v**2 - vy**2)**0.5

    a = 25.0 ! guessing the lower limit
    b = 30.0 ! guessing the upper limit
    
    ! a trial value
    !try_y2 = (a+b)/2.0

    ! number of ODEs
    n = 2

    ! initial conditions
    !y(1) =  0.0      ! y(1) = u
    !y(2) =  try_y2   ! y(2) = u' 

    do while( t .lt. tend )

        if ((t+dt) .ge. tend) then
            dt = tend - t
        endif

        try_y2 = 8.52  !(a+b)/2.0

        y(2) = try_y2

        vx = (v**2 - try_y2**2)**0.5

        tend = 50.0/vx

        !call rk2(n, y, ynext, t, h, my_func)
        call rk4(n, y, ynext, t, dt, my_func)

        y = ynext
        t = t + dt

        !if 

        print *, t, ynext(1)

    enddo
    
    print *, "y = ", ynext(1)
    print *, "The desired value is 0.0"

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







