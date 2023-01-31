%----------------------------------------------------------------------
! https://github.com/apwillis1/JFNK-Hookstep
!
! Please cite either of
! https://doi.org/10.48550/arXiv.1908.06730
! https://doi.org/10.1016/j.softx.2017.05.003 
!                                      Thanks in advance! Ashley 2023.
%----------------------------------------------------------------------
% Find root f(x) = 0 by Newton method with gmres-hookstep:
%    invert  df(x_n)/dx . s = f(x_n)  
% for s subject to constraint |s| < del  (size of 'trust region'), 
%    then     x_{n+1} = x_n - s
% To define subroutines f,df below, write as
%    df . s = f
% where f is a function to be minimised and operator df provides the
% action of the Jacobian, but may also contain any constraints on the 
% update s (e.g. no shifts in a homogenous direction), arranged such 
% that the rhs of the constraint (in f) is 0 when converged.
%----------------------------------------------------------------------
%
% Put initial guess in global new_x(1:n).  This variable is updated 
% with the current best solution each iteration.
%
% f       evaluates y=f(x):  y = f(n,x)
% df      evaluates y=Jacobian.x, see GMRESm.f90:  y = df(n,x)
% dfp     preconditioner for df, solve My=x:  y = df(n,x)
% sub     parameterless subroutine called at end of each iteration,
%         so that data could be saved etc.
% dp      dot product:  d = dp(n,a,b)
% m	  gmres dimension (also max num gmres its)
% n 	  dimension of x
% gtol    tolerence for gmres, typically 0.001
% tol     request |f(x)|<tol
% del     initial size of trust region for usefulness of df
%         set = 0d0 for no hookstep, <0d0 for del = |f(x0)|/10.
% mndl    min size of trust region
% mxdl    max size of trust region, if del<0 set no upper limit
% nits	  max num Newton its 
% info	  on input if =1 print* details
%         on exit if =0 sucessful
%                    =1 gmres failed
%                    =2 reached max iterations, nits
%                    =3 trust region got too small
%                    =4 unable to allocate memory
%----------------------------------------------------------------------
 function info = NewtonHook(f, df, dfp, sub, dp, m, n, gtol, ...
                            tol, del, mndl, mxdl, nits, info)
   global new_x		% current best guess, x
   global new_fx	%  f(x)
   global new_tol	% |f(x)|
   global new_del	% current size of trust region
   global new_nits	% number of newton iterations taken
   global new_gits	% number of GMRES iterations taken

   new_nits = 0 ;
   new_gits = 0 ;
   new_del  = del ;
   mxdl_    = mxdl ;
   ginfo    = info ;
   new_fx   = f(n,new_x) ;
   new_tol  = sqrt(dp(n,new_fx,new_fx)) ;
   if del<0d0  
      new_del = new_tol / 10d0 ;
      mxdl_   = 1d99
   end 
   if info==1 
      fprintf('newton: nits=%i  res=%e\n', new_nits, new_tol);
   end
   sub();
   x_   = new_x ;
   fx_  = new_fx ;
   tol_ = new_tol ;
   tol__ = 1d99 ;

   if new_tol<tol ;
      if info==1
         fprintf('newton: input already converged\n') ;
      end
      info = 0 ;
      return
   end
      				% - - - - Start main loop - - - - -  -
   while 1

      if new_del<mndl
         if info==1 
             fprintf('newton: trust region too small\n') ;
         end
         info = 3 ;
         return      
      end
         			% find hookstep s and update x
      s        = zeros(n,1) ;
      gres     = gtol * new_tol ;
      gdel     = new_del ;
      if ginfo ~= 2   
         new_gits = m ;
      end 
      if del==0d0
         new_gits = 9999 ;
      end
      [s,gres,gdel,new_gits,ginfo] = ...
          GMRESm(m,n,s,fx_,df,dfp,dp,gres,gdel,new_gits,ginfo) ;
      ginfo = info ;
      new_x = x_ - s ;
         			% calc new norm, compare with prediction
      new_fx  = f(n,new_x) ;
      new_tol = sqrt(dp(n,new_fx,new_fx)) ;
      snrm = sqrt(dp(n,s,s)) ;
      ared = tol_ - new_tol ;
      pred = tol_ - gdel ;
         
      if info==1 
         fprintf('newton: nits=%i  res=%e\n',new_nits,new_tol) ;
         fprintf('newton: gits=%i  del=%e\n',new_gits,new_del) ;
         fprintf('newton: |s|=%e  pred=%e\n',snrm,pred) ;
         fprintf('newton: ared/pred=%e\n',ared/pred) ;
      end

      if del==0d0 
         if info==1  
            fprintf('newton: took full newton step\n');
         end
      elseif new_tol>tol__
         if info==1 
            fprintf('newton: accepting previous step\n');
         end
         new_x   = x__ ;
         new_fx  = fx__ ;
         new_tol = tol__ ;
         new_del = del__ ;
      elseif ared<0d0
         if info==1 
            fprintf('newton: norm increased, try smaller step\n');
         end
         new_del = snrm * 0.5d0 ;
         ginfo   = 2 ;
      elseif ared/pred<0.75d0
         if info==1 
            fprintf('newton: step ok, trying smaller step\n');
         end
         x__     = new_x ;
         fx__    = new_fx ;
         tol__   = new_tol ;
         if ared/pred> 0.1d0
            del__ = snrm ;
         end
         if ared/pred<=0.1d0 
            del__ = snrm*0.5d0 ;
         end
         new_del = snrm * 0.7d0 ;
         ginfo   = 2 ;
      elseif snrm<new_del*0.9d0
         if info==1 
            fprintf('newton: step good, took full newton step\n');
         end
         new_del = min(mxdl_,snrm*2d0) ;
      elseif new_del<mxdl_*0.9d0
         if info==1 
            fprintf('newton: step good, trying larger step\n');
         end
         x__     = new_x ;
         fx__    = new_fx ;
         tol__   = new_tol ;
         del__   = new_del ;
         new_del = min(mxdl_,snrm*2d0) ;
         ginfo   = 2 ;
      end
         				% check if need to try another s
      if ginfo==2 
         continue
      end
         				% end of iteration
      new_nits = new_nits + 1 ;
      sub();
      x_   = new_x ;
      fx_  = new_fx ;
      tol_ = new_tol ;
      tol__ = 1d99 ;

      if new_tol<tol
         if info==1
            fprintf('newton: converged\n');
         end 
         info = 0 ;
         return
      elseif new_nits==nits
         if info==1
            fprintf('newton: reached max its\n');
         end
         info = 2 ;
         return
      end
   
   end

 end

