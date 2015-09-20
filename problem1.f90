program euler
 implicit none
 real(kind=8):: a,b,t,x1,x2,h,T1,T2
 real(kind=8), external ::f1,f2
 integer :: i,imax,n

 a=0.0     !  t_o
 b=100.0   ! tmax
 n=1000.0  ! number of steps
 h=(b-a)/n ! step size
 t=a
 x1=0.0    ! initial condition for theta
 x2=3.0   ! initial condition for dtheta/dt
  
 open(unit=1,file='euler_30.txt',form='formatted', status='replace',action='READWRITE' )

 print *, "	   time   		  theta			     dtheta/dt /n"
 
 do i=1,n+1
   write(unit=1,fmt="(2x,g14.6,2x,g14.6,2x,g14.6)") t,x1,x2
   write(*,'(2x,g14.6,2x,g14.6,2x,g14.6)' ) t, x1, x2
   t=t+h
   T1=f1(x2,t)
   T2=f2(x1,t)
   x1=x1+h*T1
   x2=x2+h*T2
 enddo
 close(unit=1)
end program euler

real(kind=8) function f1(x,t)
   implicit none
   real(kind=8), intent(in) :: x,t
   f1=x
end function f1

real(kind=8) function f2(x,t)
   implicit none
   real(kind=8), intent(in) :: x,t
   f2=-1.0*sin(x)
end function f2


