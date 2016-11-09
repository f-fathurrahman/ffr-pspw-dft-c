program fft_grid_size
  implicit none
  integer, parameter :: MAX_POW2=10,MAX_POW3=4,MAX_POW5=3,MAX_POW7=2
  integer :: i,j,k,l
  integer, allocatable :: grid_size(:)
  integer :: is

  allocate(grid_size((MAX_POW+1)**4))

  is = 1
  do i=0,MAX_POW
  do j=0,MAX_POW
  do k=0,MAX_POW
  do l=0,MAX_POW
    grid_size(is) = 2**i * 3**j * 5**k * 7**l
    write(*,'(4I5,I10)') i,j,k,l, grid_size(is)
    is = is + 1
  enddo
  enddo
  enddo
  enddo
  write(*,*) 'is = ', is
  write(*,*) 'size(grid_size) = ',size(grid_size)

  call Bubble(grid_size,size(grid_size))
  do i=1,size(grid_size)
    write(*,*) i, grid_size(i)
  enddo

  deallocate(grid_size)
  write(*,*) 'Program ended normally'
end program

!return p,q in ascending order
Subroutine Order(p,q)
integer p,q,temp
  if (p>q) then
    temp=p
    p=q
    q=temp
  end if
  return
end

!Buuble sorting of integer array A
Subroutine Bubble(A, n)
integer A(1:n)
  do i=1, n
    do j=n, i+1, -1
      call Order(A(j-1), A(j))
    end do
  end do
  return
end
