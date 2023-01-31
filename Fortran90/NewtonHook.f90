!----------------------------------------------------------------------
! Openpipeflow.org.  If used in your work, please cite
! Willis, A. (2017) SoftwareX 6, 124-127.
! https://doi.org/10.1016/j.softx.2017.05.003 (open access)
!                                      Thanks in advance! Ashley 2019.
!----------------------------------------------------------------------
! Find root f(x) = 0 by Newton method with gmres-hookstep:
!    invert  df(x_n)/dx . s = f(x_n)  
! for s subject to constraint |s| < del  (size of 'trust region'), 
!    then     x_{n+1} = x_n - s
! To define subroutines f,df below, write as
!    df . s = f
! where f is a function to be minimised and operator df provides the
! action of the Jacobian, but may also contain any constraints on the 
! update s (e.g. no shifts in a homogenous direction), arranged such 
! that the rhs of the constraint (in f) is 0 when converged.
!----------------------------------------------------------------------
!
! Allocate(new_x(1:n),new_fx(1:n)), then put initial guess in new_x(1:n).  
! This variable is updated with the current best solution each iteration.
!
! f       evaluates y=f(x):  call f(n,x, y)
! df      evaluates y=Jacobian.x, see GMRESm.f90:  call df(n,x, y)
! dfp     preconditioner for df:  call df(n, x)
! sub     parameterless subroutine called at end of each iteration,
!         so that data could be saved etc.
! dp      dot product:  d = dp(n,a,b)
! m	  gmres dimension (also max num gmres its)
! n 	  dimension of x
! gtol    tolerence for gmres, typically 0.001
! tol     request |f(x)|<tol
! del     initial size of trust region for usefulness of df
!         set = 0d0 for no hookstep, <0d0 for del = |f(x0)|/10.
! mndl    min size of trust region
! mxdl    max size of trust region, if del<0 set no upper limit
! nits	  max num Newton its 
! info	  on input if =1 print* details
!         on exit if =0 sucessful
!                    =1 gmres failed
!                    =2 reached max iterations, nits
!                    =3 trust region got too small
!                    =4 unable to allocate memory
!							A.P.Willis 2008
!----------------------------------------------------------------------
 module newton
   implicit none
   save
   double precision, allocatable :: new_x(:), new_fx(:)
   double precision              :: new_tol,  new_del
   integer                       :: new_nits, new_gits
 end module newton


 subroutine newtonhook(f, df, dfp, sub, dp, m, n, gtol, &
                       tol, del, mndl, mxdl, nits, info)
   use newton
   implicit none
   external                        :: f, df, dfp, sub
   double precision, external      :: dp
   integer,          intent(in)    :: m, n
   double precision, intent(in)    :: gtol, tol, del, mndl, mxdl
   integer,          intent(in)    :: nits
   integer,          intent(inout) :: info
   double precision, allocatable :: v(:,:)
   double precision :: tol_,x_(n),fx_(n), tol__,x__(n),fx__(n),del__
   double precision :: mxdl_, ared, pred, snrm, s(n)
   double precision :: gres, gdel, h((m+1)*m)
   integer :: ginfo
   
   allocate(v(n,m+1), stat=ginfo)
   if(ginfo/=0) then 
      if(info==1) print*, 'newton: unable to allocate memory'
      info = 4
      return
   end if

   new_nits = 0
   new_gits = 0
   new_del  = del
   mxdl_    = mxdl
   ginfo    = info
   call f(n,new_x, new_fx)
   new_tol  = dsqrt(dp(n,new_fx,new_fx))
   if(del<0d0)  new_del = new_tol / 10d0
   if(del<0d0)  mxdl_   = 1d99
   if(info==1)  print*,'newton: nits=',new_nits,' res=',real(new_tol)
   call sub()
   x_   = new_x
   fx_  = new_fx
   tol_ = new_tol
   tol__ = 1d99

   if(new_tol<tol) then
      if(info==1) print*, 'newton: input already converged'
      info = 0
      deallocate(v)
      return
   end if
      				! - - - - Start main loop - - - - -  -
   do while(.true.)

      if(new_del<mndl) then
         if(info==1) print*, 'newton: trust region too small'
         info = 3
         deallocate(v)
         return      
      end if
         			! find hookstep s and update x
      s        = 0d0
      gres     = gtol * new_tol
      gdel     = new_del
      if(ginfo/=2) new_gits = m
      if(del==0d0) new_gits = 9999
      call gmresm(m,n,s,fx_,df,dfp,dp,h,v,gres,gdel,new_gits,ginfo)
      ginfo = info
      new_x = x_ - s
         			! calc new norm, compare with prediction
      call f(n,new_x, new_fx)
      new_tol = dsqrt(dp(n,new_fx,new_fx))
      snrm = dsqrt(dp(n,s,s))
      ared = tol_ - new_tol
      pred = tol_ - gdel
         
      if(info==1) then 
         print*,'newton: nits=',new_nits,' res=',real(new_tol)
         print*,'newton: gits=',new_gits,' del=',real(new_del)
         print*,'newton: |s|=',real(snrm),' pred=',real(pred)
         print*,'newton: ared/pred=',real(ared/pred)
      end if

      if(del==0d0) then
         if(info==1) print*, 'newton: took full newton step'
      else if(new_tol>tol__) then
         if(info==1) print*, 'newton: accepting previous step'
         new_x   = x__
         new_fx  = fx__
         new_tol = tol__
         new_del = del__
      else if(ared<0d0) then
         if(info==1) print*, 'newton: norm increased, try smaller step'
         new_del = snrm * 0.5d0
         ginfo   = 2
      else if(ared/pred<0.75d0) then
         if(info==1)  print*, 'newton: step ok, trying smaller step'
         x__     = new_x
         fx__    = new_fx
         tol__   = new_tol
         if(ared/pred> 0.1d0)  del__ = snrm
         if(ared/pred<=0.1d0)  del__ = snrm*0.5d0
         new_del = snrm * 0.7d0
         ginfo   = 2
      else if(snrm<new_del*0.9d0) then
         if(info==1) print*, 'newton: step good, took full newton step'
         new_del = min(mxdl_,snrm*2d0)
      else if(new_del<mxdl_*0.9d0) then
         if(info==1) print*, 'newton: step good, trying larger step'
         x__     = new_x
         fx__    = new_fx
         tol__   = new_tol
         del__   = new_del
         new_del = min(mxdl_,snrm*2d0)
         ginfo   = 2
      end if
         				! check if need to try another s
      if(ginfo==2) cycle
         				! end of iteration
      new_nits = new_nits + 1
      call sub()
      x_   = new_x
      fx_  = new_fx
      tol_ = new_tol
      tol__ = 1d99

      if(new_tol<tol) then
         if(info==1) print*, 'newton: converged'
         info = 0
         deallocate(v)
         return
      else if(new_nits==nits) then
         if(info==1) print*, 'newton: reached max its'
         info = 2
         deallocate(v)
         return
      end if
   
   end do

 end subroutine newtonhook

