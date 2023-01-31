%-------------------------------------------------------------------------
%  function to be minimised
%-------------------------------------------------------------------------
 function y = getrhs(n_,x)
   global ndts
   y_ = steporbit(ndts,x) ;
   y  = y_ - x ;				% difference
   y(1) = 0d0 ;					% constraints, rhs=0
 end


