!*************************************************************************
! Openpipeflow.org.  If used in your work, please cite
! https://doi.org/10.48550/arXiv.1908.06730
! https://doi.org/10.1016/j.softx.2017.05.003 
!                                      Thanks in advance! Ashley 2019.
!*************************************************************************
!  > gfortran newton_Lorenz.f90 -llapack
!  > ./a.out
!*************************************************************************
! EXAMPLE: Periodic orbits of Lorenz system
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Newton vector, current guess new_x:
!    x(1)   = period T
!    x(2:4) = state (X,Y,Z)
!
!  Function to minimise F(x):
!    F(1)   = 0.
!    F(2:4) = X(x,T) - x,  difference between start and end point.
!
!  Newton update new_x := new_x + dx.  Constraint on update dx:
!    dx . \dot{x} = 0,  no update in direction of trajectory.
!
!  File guess.in :
!     T		initial guess for period
!     ndts	number of timesteps taken in one period
!     
!*************************************************************************
 module orbit
   implicit none
   save
   integer, parameter :: n = 4		! dimension of system, 
					!   including unknown period T
   integer            :: mgmres = 4	! max GMRES iterations
   integer            :: nits = 100	! max Newton iterations
   double precision   :: rel_err = 1d-8	! relative error |F|/|x|

   double precision :: del = -1d0	! These rarely need changing
   double precision :: mndl = 1d-20	!   for any problem
   double precision :: mxdl = 1d+20
   double precision :: gtol = 1d-3
   double precision :: epsJ = 1d-6

   double precision :: tol	
   integer :: info
   integer :: ndts	! number of timesteps to take in period T
   logical :: fixT	! Fix T for equilibrium, rather than PO solution

 end module orbit

!*************************************************************************
include 'NewtonHook.f90'
include 'GMRESm.f90'
 PROGRAM MAIN
!*************************************************************************
   use orbit
   use newton
   implicit none
   external :: getrhs, multJ, multJp, saveorbit 
   double precision, external :: dotprod
   double precision :: d

   allocate(new_x(n))
   allocate(new_fx(n))  

! 		 Initial guesses for PO, initial T,X,Y,Z, num timesteps
!   new_x = (/1.55, -13.7, -19.5, 27./)
   new_x = (/1.5, -13., -19., 27./)
   ndts = 300
   fixT = .false.

!		 Initial guess for equilibrium.  Take T short and fixed.
!   new_x = (/0.2, 13., 19., 28./)
!   ndts = 20
!   fixT = .true.

			! scale params by |x|
   d = dotprod(-1,new_x,new_x)
   tol  = rel_err * dsqrt(d)
   del  = del     * dsqrt(d)
   mndl = mndl    * dsqrt(d)
   mxdl = mxdl    * dsqrt(d)

   info = 1
   call newtonhook(getrhs, multJ, multJp, saveorbit, dotprod, &
                   mgmres, n, gtol, tol, del, mndl, mxdl, nits, info)
   print*, 'new_x = '
   print*, real(new_x)
			! check solution
   if(fixT) then
      d = sqrt((28d0-1d0)*(8d0/3d0)) 
      print*, 'Eqm = '
      print*, real(d), real(d), 28.-1.
   else    
      call steporbit(ndts,new_x, new_fx) ;
      print*, 'f(new_x) = '
      print*, real(new_fx(2:))
   end if

   stop

 contains

!*************************************************************************
 END PROGRAM MAIN
!*************************************************************************
!-------------------------------------------------------------------------
!  Lorenz system: dx/dt = f(x)
!            state:  x(1)==X, x(2)==Y, x(3)==Z.
!       parameters:  p(1)==s, p(2)==r, p(3)==b  (classic: s=10,r=28,b=8/3)
!-------------------------------------------------------------------------
 subroutine Lorenz_f(x, f)
    implicit none
    double precision, intent(in)  :: x(3)
    double precision, intent(out) :: f(3)
    double precision :: p(3) = (/10d0, 28d0, 8d0/3d0 /)
    f(1) = p(1)*(-x(1)+x(2)) ;
    f(2) = x(1)*(-x(3)+p(2)) - x(2) ;
    f(3) = x(1)*x(2) - p(3)*x(3) ;
 end
   
!-------------------------------------------------------------------------
!  function to be minimised   
!-------------------------------------------------------------------------
 subroutine getrhs(n_,x, y)
   use orbit
   implicit none
   integer,          intent(in)  :: n_
   double precision, intent(in)  :: x(n)
   double precision, intent(out) :: y(n)
   double precision :: x_(n), y_(n)

   x_ = x
   call steporbit(ndts,x_, y_)
   y = y_ - x					! diff
   y(1) = 0d0					! constraints, rhs=0

 end subroutine getrhs


!-------------------------------------------------------------------------
!  Action of Jacobian and constraints on update.
!  Use approximation  dF(x_n)/dx . dx = (F(x_n+eps.dx)-F(x_n))/eps
!-------------------------------------------------------------------------
 subroutine multJ(n_,x, y)
   use newton
   use orbit,    only : n, epsJ, ndts, fixT
   implicit none
   integer,          intent(in)  :: n_
   double precision, intent(in)  :: x(n)
   double precision, intent(out) :: y(n)   
   double precision, external :: dotprod
   double precision :: eps, s(n), dt
    				! (F(x0+eps.x)-F(x0))/eps
   eps = dsqrt(dotprod(1,x,x))
   if(eps==0d0)  stop 'multJ: eps=0 (1)'
   eps = epsJ * dsqrt(dotprod(1,new_x,new_x)) / eps
   if(eps==0d0)  stop 'multJ: eps=0 (2)'
   y = new_x + eps*x
   call getrhs(n_,y, s)
   y = (s - new_fx) / eps

   if(fixT) then		! no extra constraint if T fixed
      y(1) = 0d0
   else				! contstraint, 
				! no update in trajectory direction
      call steporbit(1,new_x, s)
      dt = new_x(1)/ndts
      s = (s - new_x) / dt
      y(1) = dotprod(-1,s,x)
   end if

 end subroutine multJ
 

!-------------------------------------------------------------------------
!  preconditioner for multJ.  Empty - no preconditioner required
!-------------------------------------------------------------------------
 subroutine multJp(n, x)
   implicit none
   integer,          intent(in)    :: n
   double precision, intent(inout) :: x(n)
 end subroutine multJp


!-------------------------------------------------------------------------
!  called at each newton iteration   
!-------------------------------------------------------------------------
 subroutine saveorbit()
   use newton
   use orbit
   implicit none
   double precision :: norm_x, p
   double precision, external :: dotprod

   print*, 'newton: iteration ', new_nits

   norm_x = dsqrt(dotprod(-1,new_x,new_x))
   print*, 'relative error = ', new_tol/norm_x

!  SAVE current solution, new_x
   
 end subroutine saveorbit
 
 
!-------------------------------------------------------------------------
! dot product.  can flag to exclude parameter T.  Could include weights
!-------------------------------------------------------------------------
 double precision function dotprod(n_,a,b)
   use orbit
   implicit none
   integer,          intent(in) :: n_
   double precision, intent(in) :: a(n), b(n)
   integer :: n1
   n1 = 1
   if(n_==-1) n1 = 2
   dotprod = dot_product(a(n1:n),b(n1:n))
 end function dotprod


!-------------------------------------------------------------------------
!  timestep
!-------------------------------------------------------------------------
 subroutine steporbit(ndts_,x, y)
   use orbit
   implicit none
   integer,          intent(in)  :: ndts_
   double precision, intent(in)  :: x(n)
   double precision, intent(out) :: y(n)
   double precision, save :: dt
   double precision :: a(3), a1(3), fa(3), fa1(3)
   integer :: i
  
   if(ndts_/=1) then		! Set timestep size dt=T/ndts_
      dt = x(1) / dble(ndts_)	! If only doing one step to calc \dot{x},
   end if			! then use previously set dt.

   a = x(2:)
			! second-order predictor-corrector method
   do i = 1, ndts_
      call Lorenz_f(a, fa)
      a1 = a + dt*fa
      call Lorenz_f(a1, fa1)
      a = a + 0.5d0*dt*(fa+fa1)
   end do
   
   y(1)  = 0d0
   y(2:) = a

 end subroutine steporbit
 
