%-------------------------------------------------------------------------
%  Action of Jacobian on update dx, and evaluation of constraint on update.
%  Use approximation    dF(x_n)/dx . dx = (F(x_n+eps.dx)-F(x_n))/eps
%-------------------------------------------------------------------------
 function y = multJ(n_,dx)
   global new_x
   global new_fx
   global epsJ
   global ndts
   global fixT
    				% (F(x0+eps.x)-F(x0))/eps
   eps = sqrt(dotprd(1,dx,dx)) ;
   eps = epsJ * sqrt(dotprd(1,new_x,new_x)) / eps ;
   y = new_x + eps*dx ;
   s = getrhs(n_,y) ;
   y = (s - new_fx) / eps ;

   if fixT			% no extra constraint if T fixed
      y(1) = 0d0 ;
   else				% constraint: dx . \dot{x} = 0
				% no update in trajectory direction
      s = steporbit(1,new_x) ;
      dt = new_x(1)/ndts ;
      s = (s - new_x) / dt ;
      y(1) = dotprd(-1,s,dx) ;
   end

 end 


