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

            stars2(1)%mass = m1
            stars3(1)%mass = m1
            stars4(1)%mass = m1
            !
            ! setup initial conditions of star M2
            stars(2)%mass = m2
            stars(2)%x    = 1.0*au
            stars(2)%y    = 0.0
            stars(2)%vx   = 0.0
            stars(2)%vy   = - ( 2.0 * pi * 1.0*au / period )
            stars(2)%ax   = ( - force/m2 )
            stars(2)%ay   = 0.0          

            stars2(2)%mass = m2
            stars3(2)%mass = m2
            stars4(2)%mass = m2

        end subroutine initial

        subroutine update(dt)
            use constants
            implicit none
            real, intent(in)  :: dt
            integer :: i,j
            real    :: x, y, rsq, fx, fy, h
            real    :: radius, force, angle


            !RK4 y_k+1 = y_k + h_k * (k1+2k2+2k3+k4) * 1/6
            
            ! Initial condition is given in previous subroutine
              
            ! update position to t = t + dt
           
            h = 0.5*dt
            
            !call update_acc(h, stars(1), stars(2), stars2(1), stars2(2))
            call update_euler(h, stars(1), stars(1), stars2(1))
            call update_euler(h, stars(2), stars(2), stars2(2))           ! update x2(i) y2(i) vx2(i) vy2(i)
            call update_acc(stars2(1), stars2(2), stars2(1), stars2(2))   ! update ax2(i) ay2(i)
            
            !print *, stars2(1)%vx, stars2(2)%vx

            h = 0.5*dt
            
            !call update_acc(h, stars2(1), stars2(2), stars3(1), stars3(2))
            call update_euler(h, stars(1), stars2(1), stars3(1))
            call update_euler(h, stars(2), stars2(2), stars3(2))          ! update x3(i) y3(i) vx3(i) vy3(i)
            call update_acc(stars3(1), stars3(2), stars3(1), stars3(2)) ! update ax3(i) ay3(i)

            h = dt
            
            !call update_acc(h, stars3(1), stars3(2), stars4(1), stars4(2))
            call update_euler(h, stars(1), stars3(1), stars4(1))
            call update_euler(h, stars(2), stars3(2), stars4(2))          ! update x4(i) y4(i) vx4(i) vy4(i)
            call update_acc(stars4(1), stars4(2), stars4(1), stars4(2)) ! update ax4(i) ay4(i)

            !call update_acc(h, stars4(1), stars4(2), stars2(1), stars2(2) !k3
 
            stars(1)%x = stars(1)%x + dt*(1.0*stars(1)%vx + 2.0*stars2(1)%vx + 2.0*stars3(1)%vx + 1.0*stars4(1)%vx)/6.0
            stars(2)%x = stars(2)%x + dt*(1.0*stars(2)%vx + 2.0*stars2(2)%vx + 2.0*stars3(2)%vx + 1.0*stars4(2)%vx)/6.0
              
            stars(1)%y = stars(1)%y + dt*(1.0*stars(1)%vy + 2.0*stars2(1)%vy + 2.0*stars3(1)%vy + 1.0*stars4(1)%vy)/6.0
            stars(2)%y = stars(2)%y + dt*(1.0*stars(2)%vy + 2.0*stars2(2)%vy + 2.0*stars3(2)%vy + 1.0*stars4(2)%vy)/6.0

            stars(1)%vx = stars(1)%vx + dt*(1.0*stars(1)%ax + 2.0*stars2(1)%ax + 2.0*stars3(1)%ax + 1.0*stars4(1)%ax)/6.0
            stars(2)%vx = stars(2)%vx + dt*(1.0*stars(2)%ax + 2.0*stars2(2)%ax + 2.0*stars3(2)%ax + 1.0*stars4(2)%ax)/6.0

            stars(1)%vy = stars(1)%vy + dt*(1.0*stars(1)%ay + 2.0*stars2(1)%ay + 2.0*stars3(1)%ay + 1.0*stars4(1)%ay)/6.0
            stars(2)%vy = stars(2)%vy + dt*(1.0*stars(2)%ay + 2.0*stars2(2)%ay + 2.0*stars3(2)%ay + 1.0*stars4(2)%ay)/6.0
            
            call update_acc(stars(1), stars(2), stars(1), stars(2))

            return
        end subroutine update

    subroutine update_euler(dt,s1,s2,s3)
            !!
            !! use s1 and s2 -> s3
            !!
            use constants
            implicit none
            real, intent(in) :: dt
            type(Star),intent(in)  :: s1
            type(Star),intent(in)  :: s2
            type(Star),intent(out) :: s3

            ! real    :: x, y, rsq, fx, fy
            ! real    :: radius, force, angle

            ! update position and velocity
           
            s3%x      = s1%x  +s2%vx *dt
            s3%y      = s1%y  +s2%vy *dt
            s3%vx     = s1%vx +s2%ax *dt
            s3%vy     = s1%vy +s2%ay *dt

        end subroutine update_euler 

        subroutine update_acc(s1,s2,s3,s4)
            !!
            !! use s1 -> s2
            !!
            use constants
            implicit none
            

            type(Star),intent(in)  :: s1, s2
            type(Star),intent(out) :: s3, s4

            real    :: angle,x,y,rsq,force,fx,fy

            x = s1%x - s2%x
            y = s1%y - s2%y

            !print *, x

            rsq = x**2.0 + y**2.0

            angle = atan2(y,x) ! atan2(y,x)接受兩個參數，atan(slope)只接受一個數字 

            force = G * s1%mass * s2%mass / rsq

            fx = force * cos(angle)
            fy = force * sin(angle)

            s3%ax = - fx/s1%mass
            s3%ay = - fy/s1%mass

            s4%ax = fx/s2%mass
            s4%ay = fy/s2%mass

        end subroutine update_acc

end module physics

