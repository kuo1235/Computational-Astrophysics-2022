module Simulation_data
    implicit none
    integer, parameter :: imax   = 128   ! number of points in the x direction
    integer, parameter :: ibuf   = 1     ! number of ghost zones for BC.
    integer, parameter :: istart = 1     ! starting point
    integer, parameter :: iend   = imax  ! end point

    integer, parameter :: jmax   = 128   ! number of points in the y direction
    integer, parameter :: jbuf   = 1     ! number of ghost zones for BC.
    integer, parameter :: jstart = 1     ! starting point
    integer, parameter :: jend   = jmax  ! end point

    real, parameter  :: cx       = 1.0   ! x-velocity
    real, parameter  :: cy       = 1.0   ! y-velocity
    real, parameter  :: xmin     = 0.0   
    real, parameter  :: xmax     = 1.0
    real, parameter  :: ymin     = 0.0
    real, parameter  :: ymax     = 1.0   
    real, parameter  :: tend     = 2.5   ! final time

    real, parameter  :: cfl      = 0.4   ! cfl number

    real, save       :: dx, dy

    real, dimension(istart-ibuf:iend+ibuf,jstart-jbuf:jend+jbuf), save  :: u, uold
    real, dimension(istart-ibuf:iend+ibuf), save  :: x
    real, dimension(jstart-jbuf:jend+jbuf), save  :: y

    integer,parameter  :: io_interval = 10

end module Simulation_data
