program t_ddot
  implicit none
  integer, parameter :: N=4
  real(8) :: x(N), y(N)
  real(8) :: res
  ! Function
  real(8) :: ddot

  x(:) = 1.d0
  y(:) = 2.d0

  res = ddot(N,x,1,y,1)

  write(*,*) 'res = ', res

end program
