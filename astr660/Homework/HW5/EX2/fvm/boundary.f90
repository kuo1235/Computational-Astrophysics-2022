subroutine boundary(v)
    !
    ! apply BC on array v
    !
    use Simulation_data
    implicit none
    real, dimension(istart-ibuf:iend+ibuf), intent(inout) :: v
    integer :: i

    ! apply boundary condition
    
    ! ghost zone ibuf=1

    v(istart-1) = v(iend)     !implement periodic BC for left boundary 
     
    v(iend+1) = v(istart)     !implement periodic BC for right boundary 
    
   

end subroutine boundary
