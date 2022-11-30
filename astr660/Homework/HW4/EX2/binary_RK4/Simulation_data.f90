module Simulation_data
    implicit none

    integer,parameter :: N=2    ! number of stars

    type Star
        integer :: id  
        real :: mass
        real :: x
        real :: y
        real :: vx
        real :: vy
        real :: ax
        real :: ay
    end type Star

    type(Star), dimension(N) :: stars, stars2, stars3, stars4
    real :: separation
    real :: period

end module Simulation_data
