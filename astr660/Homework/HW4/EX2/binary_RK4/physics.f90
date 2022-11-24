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

            real, dimension(2) :: k1, k2, k3, k4   


            !RK4 y_k+1 = y_k + h_k * (k1+2k2+2k3+k4) * 1/6
            
            k1(i) = 


            ! update position to t = t + dt
            do i=1, 2

                stars(i)%x = stars(i)%x + 

            ! update velocity to t = t + dt
            !TODO

            ! update accelerations to t = t + dt
            !TODO

            return
        end subroutine update

end module physics

