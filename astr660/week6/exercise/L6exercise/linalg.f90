module linalg

    ! ---------------------------------------------------------
    ! ref: https://rosettacode.org/wiki/LU_decomposition#Fortran
    ! ---------------------------------------------------------

    implicit none
    contains 


        subroutine mat_print(amsg,a)
            character(*), intent(in) :: amsg
            class    (*), intent(in) :: a(:,:)
            integer                  :: i
            print*,' '
            print*,amsg
            do i=1,size(a,1)
                select type (a)
                    type is (real(8)) ; print'(100f8.3)',a(i,:)
                    type is (integer) ; print'(100i8  )',a(i,:)
                end select
            end do
            print*,' '
        end subroutine

        subroutine solve_lower_triangular_matrix(N,L,b,x)
            implicit none

            integer, intent(in)  :: N
            real, dimension(N,N), intent(in)  :: L  ! lower triangle
            real, dimension(N)  , intent(in)  :: b  ! vector
            real, dimension(N)  , intent(out) :: x  ! solution
            real, dimension(N)  :: bs              

            integer :: i,j
             
            bs = b 
            do j=1,N
                if (L(j,j) == 0.) then
                    stop
                end if

                x(j) = bs(j) / L(j,j)    
                
                do i = j + 1, N
                    bs(i) = bs(i) - L(i,j) * x(j)    
                enddo    

            enddo

            return
        end subroutine solve_lower_triangular_matrix

        subroutine solve_upper_triangular_matrix(N,U,b,x)
            implicit none

            integer, intent(in)  :: N
            real, dimension(N,N), intent(in)  :: U  ! upper triangle
            real, dimension(N)  , intent(in)  :: b  ! vector
            real, dimension(N)  , intent(out) :: x  ! solution
            real, dimension(N)  :: bs

            integer :: i,j

            bs = b

            do j=N,1,-1
                if (U(j,j) == 0.) then
                    stop
                end if    
                x(j) = bs(j) / U(j,j)
                do i= 1, j - 1
                    bs(i) = bs(i) - U(i,j) * x(j)
                enddo        
    
            enddo

            return
        end subroutine solve_upper_triangular_matrix

        subroutine LU_decomposition(N,A,L,U)
            implicit none
            
            integer, intent(in)  :: N
            real, dimension(N,N), intent(in)  :: A    ! matrix
            real, dimension(N,N), intent(out)    :: L,U  ! matrix
            real, dimension(N,N) :: M, As
 
            integer :: i,j,k
             
            forall (j=1:N, i=1:N)
                L(i,j)=merge(1.0,0.0,i==j) ! if i=j, output:1, otherwise:0 -> creating identity matrix
                U(i,j)=0.0
            end forall

            As = A
            
            do k = 1,n-1
                if (As(k,k)==0.) then
                    stop
                end if
                
                do i = k+1,n
                    M(i,k)=As(i,k)/As(k,k)
                enddo

                do j = k+1,n
                    do i = k+1, n
                        As(i,j)=As(i,j)-M(i,k)*As(k,j)
                    enddo
                enddo                 
            enddo

            do i = 1, N
                L(i, :i-1) = M(i, :i-1)
                U(i,i: ) = As(i,i: ) 
            enddo

        end subroutine LU_decomposition

        subroutine solve_lu(N,A,b,x)
            ! solve: A x = b
            integer, intent(in)  :: N
            real, dimension(N,N), intent(in)  :: A  ! upper triangle
            real, dimension(N)  , intent(in)  :: b  ! vector
            real, dimension(N)  , intent(out) :: x  ! solution

            real, dimension(N,N) :: L, U, P
            real, dimension(N)   :: y, pb
            
            !A(1,1)=2.
            !A(1,2)=4.
            !A(1,3)=-2.
         




        end subroutine solve_lu

end module linalg
