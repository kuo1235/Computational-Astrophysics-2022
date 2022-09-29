program ex1_f

integer :: result, i
double precision, dimension(100) ::numbers

do i=1, 100
    numbers(i) = i
end do 

!print *, numbers

result = sum(numbers)
print *, result

end program ex1_f
