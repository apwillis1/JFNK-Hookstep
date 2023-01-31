%-------------------------------------------------------------------------
%  called each newton iteration   
%-------------------------------------------------------------------------
 function saveorbit()
   global new_x
   global new_fx
   global new_tol
   global new_del
   global new_nits
   global new_gits

   fprintf('newton: iteration %i\n', new_nits) ;

   norm_x = sqrt(dotprd(-1,new_x,new_x)) ;   
   relative_err = new_tol/norm_x

   % SAVE current solution,  new_x
   
 end
 

