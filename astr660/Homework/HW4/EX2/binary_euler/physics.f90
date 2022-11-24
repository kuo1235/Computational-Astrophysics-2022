!---------------------------------------------------
! The physics module
!
module physics
    use Simulation_data
    implicit none
    contains
        
        subroutine initial()
            !
            ! setup initial conditions of each stars
            ! in this example we only have two stars
            !
            use constants, only : au, msun, pi, G
            implicit none
            integer :: i
            real :: m1, m2, force

            m1 = 1.0 * msun
            m2 = 2.0 * msun

            !
            ! Use Kepler's law to evaluate the orbital period
            ! and use Newton's law to evaluate the force between
            ! two stars
            !

            separation = 3.0*au
            period     = (((4.0 * pi**2.0) / (G*(m1+m2)) ) * (separation**3.0))**0.5
            force      = (G * m1 * m2) / (separation**2.0) 

            !
            ! setup initial conditions of star M1
            stars(1)%mass = m1
            stars(1)%x    = - 2.0*au
            stars(1)%y    = 0.0
            stars(1)%vx   = 0.0
            stars(1)%vy   = (2.0 * pi * 2.0*au / period) ! v1 = 2*pi*r1/T , r1 = m2/(m1+m2) * separation 
            stars(1)%ax   = ( force/m1 )
            stars(1)%ay   = 0.0 

            !
            ! setup initial conditions of star M2
            stars(2)%mass = m2
            stars(2)%x    = 1.0*au
            stars(2)%y    = 0.0
            stars(2)%vx   = 0.0 
            stars(2)%vy   = - ( 2.0 * pi * 1.0*au / period ) 
            stars(2)%ax   = ( - force/m2 )  
            stars(2)%ay   = 0.0 


        end subroutine initial

        subroutine update(dt)
            use constants
            implicit none
            real, intent(in)  :: dt
            integer :: i,j
            real    :: x, y, rsq, fx, fy
            real    :: radius, force, angle

            !
            ! In this example we use a first order scheme (Euler method)
            ! we approximate dx/dt = v  --> x^(n+1) - x^n = v * dt
            ! therefore, x at step n+1 = x^(n+1) = x^n + v * dt
            !
            ! the same approximation can be applied to dv/dt = a
            !
            do i= 1, 2
            ! update position to t = t + dt
                stars(i)%x = stars(i)%x + stars(i)%vx * dt           
                stars(i)%y = stars(i)%y + stars(i)%vy * dt

            ! update velocity to t = t + dt
                stars(i)%vx = stars(i)%vx + stars(i)%ax * dt
                stars(i)%vy = stars(i)%vy + stars(i)%ay * dt
            
            enddo
            ! update accelerations to t = t + dt
            
            x = stars(1)%x - stars(2)%x
            y = stars(1)%y - stars(2)%y

            rsq = x**2.0 + y**2.0

            angle = atan2(y,x) ! atan2(y,x)接受兩個參數，atan(slope)只接受一個數字 

            force = G * stars(1)%mass * stars(2)%mass / rsq

            fx = force * cos(angle)
            fy = force * sin(angle)

            stars(1)%ax = - fx/stars(1)%mass
            stars(1)%ay = - fy/stars(1)%mass

            stars(2)%ax =   fx/stars(2)%mass
            stars(2)%ay =   fy/stars(2)%mass

            return
        end subroutine update

end module physics

