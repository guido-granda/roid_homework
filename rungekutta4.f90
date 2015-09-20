program problem2
 implicit none
 real(kind=8):: a,b,h,k
 real(kind=8), dimension(2):: x1,func,t
 integer ::  n
 a=0.0
 b=100.0
 n=1000.0
 h=(a-b)/n
 x1(1)=0.0
 x1(2)=0.0
 k=1.0
 t(1)=a
 t(2)=a 
 call f(x1,t,k,func) 
 call runge4(func,t,x1,h,n)

end program problem2
subroutine runge4(f,t,x,h,nsteps)
          
	implicit none
	!n: number of equations of the sytem
	!nsteps: number of steps
   	real(kind=8), dimension(2),intent(inout) ::x,f
	real(kind=8), dimension(2) ::k1,k2,k3,k4
        real(kind=8), dimension(2) :: ta
        real(kind=8), dimension(2),intent(inout) :: t
        real(kind=8), intent(inout) ::  h
        integer ::nsteps,j
        ta=t
	do j=1,nsteps
        !  do i=1,2
		k1=h*f(t,x)
		k2=h*f(t+0.5*h,x+0.5*k1)
		k3=h*f(t+0.5*h,x+0.5*k2)
		k4=h*f(t+0.5*h,x+k3)
		x=x+1.0/6.0*(k1+2.0*k2+2.0*k3+k4)
		t=ta+j*h
		write(*,*) "n 		t		x"
	!  enddo
        enddo
end subroutine runge4

subroutine f(x,t,k,func)
   implicit none
   real(kind=8), intent(inout) :: k
   real(kind=8), dimension(2), intent(in) :: x
   real(kind=8), dimension(2) ,intent(out) :: func,t
   func(1)=x(2)
   func(2)=k*sin(x(1))
end subroutine f

