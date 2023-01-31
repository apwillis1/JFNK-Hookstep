%-------------------------------------------------------------------------
%  advance x by ndts_ timesteps
%-------------------------------------------------------------------------
 function y = steporbit(ndts_,x)
   persistent dt

   if ndts_ ~= 1		% Set timestep size dt=T/ndts_
      dt = x(1) / ndts_ ;	% If only doing one step to calc \dot{x},
   end				% then use previously set dt.

   a = x(2:end) ;
			% second-order predictor-corrector method
   for n = 1:ndts_
      fa  = Lorenz_f(a) ;
      a1  = a + dt*fa ;
      fa1 = Lorenz_f(a1) ;
      a  = a + 0.5d0*dt*(fa+fa1) ;
   end

   y = zeros(size(x)) ;
   y(2:end) = a ;
 end 
 

