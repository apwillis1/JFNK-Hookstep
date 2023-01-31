%----------------------------------------------------------------------
% Openpipeflow.org.  If used in your work, please cite
% Willis, A. (2017) SoftwareX 6, 124-127.
% https://doi.org/10.1016/j.softx.2017.05.003 (open access)
%                                      Thanks in advance! Ashley 2019.
%----------------------------------------------------------------------
% solve A x = b for x ;  
% minimise |Ax-b| subject to constraint |x| < delta .
% requires lapack routines dgelsy, dgesvd.
%----------------------------------------------------------------------
% m	  gmres dimension
% n 	  dimension of x
% x	  on input:  guess for x, can be 0
%         on exit:  solution x, subject to constraint if del>0
% b	  input b
% matvec  performs y = A x:  y = matvec(N,x)
% psolve  preconditioner, solve M y = x:  y = psolve(N,x)
% dotprd  dot product, d = dotprd(n,a,b)
% h       Hessian matrix,  size (m+1)*m
% v       Krylov subspace, size n*(m+1)
% res	  on input: |Ax-b|/|b|<res; 
%         on exit:  residual reached
% del     on input: if(del>0) then the x returned is the hookstep
%         on exit:  norm of next b predicted by hook
% its	  on input: max num its; 
%         on exit:  number of its taken
% info	  on input: if(info==1) print* residuals
%                   if(info==2) recalc hookstep with new del
% 	  on exit:  0 sucessful, 1 method breakdown, 2 max its
%----------------------------------------------------------------------

 function [x,res,del,its,info] = ...
                  GMRESm(m,n,x,b,matvec,psolve,dotprd,res,del,its,info)
   persistent h
   persistent v
   persistent beta_
   persistent j_

   if info == 2  
      [y,del] = GMREShook(j_,h,m,beta_,del) ;
      z = v(:,1:j_)*y(1:j_) ;
      x = psolve(n,z) ;
      info = 0 ;
      return
   end	 

   tol = res ;
   imx = its ;
   its = 0 ;
   v   = zeros(n,m+1) ;

   while 1 %(restart)
   res_ = 1d99 ;
   stgn = 1d0 - 1d-14 ;
 
   beta_ = sqrt(dotprd(n,x,x)) ;
   if beta_ == 0d0 
      w = 0d0 ;
   else
      w = matvec(n,x) ;
   end 
   w = b - w ;
   beta_ = sqrt(dotprd(n,w,w)) ;
   v(:,1) = w / beta_ ;
     
   h = zeros(m+1,m) ;
   for j = 1:m
      j_ = j ;
      its = its + 1 ;
      z = v(:,j) ; 
      z = psolve(n,z) ;
      w = matvec(n,z) ;
      for i = 1:j
         h(i,j) = dotprd(n,w,v(:,i)) ;
         w = w - h(i,j)*v(:,i) ;
      end
      h(j+1,j) = sqrt(dotprd(n,w,w)) ;
      v(:,j+1) = w / h(j+1,j) ;
          
      p = zeros(j+1,1) ;
      p(1) = beta_ ;
      h_(1:j+1,1:j) = h(1:j+1,1:j) ;
      y = pinv(h_)*p;
      
      p = - h(1:j+1,1:j)*y ;
      p(1) = p(1) + beta_ ;
      res  = sqrt(sum(p.*p)) ;
      if info==1 
          fprintf('gmresm: it=%i  res=%e\n', its, res) ; 
      end
      
      done = ((res <= tol) || (its == imx) || (res > res_)) ;
      if done || j==m
         if del > 0d0  
            [y,del] = GMREShook(j,h,m,beta_,del) ;
         end
         z = v(:,1:j)*y(1:j) ;
         z = psolve(n,z) ;
         x = x + z ;
         if its==imx
            info = 2 ;
         end
         if res>res_
            info = 1 ;
         end
         if res<=tol
            info = 0 ;
         end 
         if done 
            return 
         end
         if del>0d0   
            fprintf('gmres: WARNING: m too small. restart affects hookstep.\n');
         end
      end
      res_ = res*stgn ;

   end %for
   end %while 
 
 end
 

