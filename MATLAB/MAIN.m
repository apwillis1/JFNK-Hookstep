%*************************************************************************
% https://github.com/apwillis1/JFNK-Hookstep
%
% Please cite either of
% https://doi.org/10.48550/arXiv.1908.06730
% https://doi.org/10.1016/j.softx.2017.05.003 
%                                      Thanks in advance! Ashley 2023.
%*************************************************************************
% EXAMPLE: Periodic orbits of Lorenz system.
%- - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - -
%  Newton vector, current guess x_n:
%    x(1)     = T
%    x(2:end) = (x,y,z)
%
%  Function to minimise F(x):
%    F(1)     = 0.
%    F(2:end) = X(x,T) - x, difference between start and end point.
%
%  Newton update: x_{n+1} = x_n + dx.  Constraint on update dx:
%    dx . \dot{x} = 0, no update in direction of trajectory.
%
%*************************************************************************
%*************************************************************************
% PROGRAM MAIN
%*************************************************************************
   global new_x		% Current best x
   global epsJ		% epsilon used in Jacobian approximation
   global ndts		% Number of timesteps taken in period T
   global fixT		% Fix T for equilibrium, rather than PO solution
   global p		% Parameters of dynamical system

   n       = 4 ;	% Dimension of system, including unknown params
   mgmres  = 4 ;	% max GMRES iterations
   nits    = 100 ;	% max Newton iterations
   rel_err = 1d-8 ;	% Relative error |F|/|x|

   del     = -1d0 ;	% These rarely need changing for any problem
   mndl    = 1d-20 ;
   mxdl    = 1d+20 ;
   gtol    = 1d-3 ;
   epsJ    = 1d-6 ;
			% Lorenz parameters
   p = [10; 28; 8/3]

% 		% Initial guesses for PO, initial T,X,Y,Z, num timesteps
%   new_x = [1.558; -13.76; -19.58; 27.] ;
%   new_x = [1.55; -13.7; -19.5; 27.] ;
   new_x = [1.5; -13.; -19.; 27.] ;
   ndts = 300 ;
   fixT = 0 ;

%		% Initial guess for equilibrium.  Take T short and fixed.
%   new_x = [0.2; 13.; 19.; 28.] ;
%   ndts = 20 ;
%   fixT = 1 ;

%		% plot initial guess
   start_x = new_x
   end_x = steporbit(ndts,start_x)
   x = zeros(4,ndts) ;
   x(:,1) = start_x ;
   for i = 1:ndts-1 ;
      x(:,i+1) = steporbit(1,x(:,i)) ;
   end
   hold on
   plot(x(2,:),x(4,:),'g','LineWidth',1) % X,Z
   plot(x(2,1),x(4,1),'go','LineWidth',1)

%		% scale parameters by |x| then call black box
   d = sqrt(dotprd(-1,new_x,new_x)) ;
   tol  = rel_err * d ;
   del  = del     * d ;
   mndl = mndl    * d ;
   mxdl = mxdl    * d ;

   info = 1 ;
   info = NewtonHook(@getrhs, @multJ, @multJp, @saveorbit, @dotprd, ...
                   mgmres, n, gtol, tol, del, mndl, mxdl, nits, info) ;

%		% plot final solution
   start_x = new_x
   end_x = steporbit(ndts,start_x)
   x = zeros(4,ndts) ;
   x(:,1) = start_x ;
   for i = 1:ndts-1 ;
      x(:,i+1) = steporbit(1,x(:,i)) ;
   end
   plot(x(2,:),x(4,:),'b','LineWidth',2) % X,Z
   plot(x(2,1),x(4,1),'bo','LineWidth',2)
   pause

   
%*************************************************************************
% END PROGRAM MAIN
%*************************************************************************


