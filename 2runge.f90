program problem2
 implicit none
 real(kind=8):: t0,tmax,dt,t1
 real(kind=8), dimension(2):: u0,u1
 integer ::  n,n2
 external :: f
 t0=0.0
 tmax=100.0
 n=1000.0
 dt=(tmax-t0)/n
 u0(1)=0.0
 u0(2)=3.0
 n2=2
! k=1.0
! t(1)=a
! t(2)=a 
  
 open(unit=1,file='runge_30.txt',form='formatted', status='replace',action='READWRITE' )
 do
   write(unit=1,fmt="(2x,g14.6,2x,g14.6,2x,g14.6)") t0,u0(1),u0(2) 
   write(*,'(2x,g14.6,2x,g14.6,2x,g14.6)' ) t0, u0(1), u0(2)
   if (tmax <= t0) then
    exit
   end if
   t1=t0+dt 
   call rk4vec(t0,2,u0,dt,f,u1)
   t0=t1
   u0(1:n2)=u1(1:n2)
 end do
 close(unit=1)
 return

end program problem2
subroutine rk4vec ( t0, m, u0, dt, f, u )

!*****************************************************************************80
!
!! RK4VEC takes one Runge-Kutta step for a vector ODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T0, the current time.
!
!    Input, integer ( kind = 4 ) M, the dimension of the system.
!
!    Input, real ( kind = 8 ) U0(M), the solution estimate at the current time.
!
!    Input, real ( kind = 8 ) DT, the time step.
!
!    Input, external F, a subroutine of the form 
!      subroutine f ( t, m, u, uprime ) 
!    which evaluates the derivative UPRIME(1:M) given the time T and
!    solution vector U(1:M).
!
!    Output, real ( kind = 8 ) U(M), the fourth-order Runge-Kutta solution 
!    estimate at time T0+DT.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) dt
  external f
  real ( kind = 8 ) f0(m)
  real ( kind = 8 ) f1(m)
  real ( kind = 8 ) f2(m)
  real ( kind = 8 ) f3(m)
  real ( kind = 8 ) t0
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) u(m)
  real ( kind = 8 ) u0(m)
  real ( kind = 8 ) u1(m)
  real ( kind = 8 ) u2(m)
  real ( kind = 8 ) u3(m)
!
!  Get four sample values of the derivative.
!
  call f ( t0, m, u0, f0 )

  t1 = t0 + dt / 2.0D+00
  u1(1:m) = u0(1:m) + dt * f0(1:m) / 2.0D+00
  call f ( t1, m, u1, f1 )

  t2 = t0 + dt / 2.0D+00
  u2(1:m) = u0(1:m) + dt * f1(1:m) / 2.0D+00
  call f ( t2, m, u2, f2 )

  t3 = t0 + dt
  u3(1:m) = u0(1:m) + dt * f2(1:m)
  call f ( t1, m, u3, f3 )
!
!  Combine them to estimate the solution U at time T1.
!
  u(1:m) = u0(1:m) + dt * ( f0(1:m) + 2.0D+00 * f1(1:m) + 2.0D+00 * f2(1:m) &
    + f3(1:m) ) / 6.0D+00
  return
end
subroutine f(t,n,u,uprime)
   implicit none
   integer (kind=4) :: n
   real(kind=8) ::t,u(n),uprime(n) 
   real(kind=8) :: k
   k=1.0
   uprime(1)=u(2)
   uprime(2)=-k*sin(u(1))
   return
end subroutine f
