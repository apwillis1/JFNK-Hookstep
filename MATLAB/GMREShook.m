%-----------------------------------------------------------------
% Called by GMRESm.  Return hookstep vector y s.t. approx |y|=del.
% c.f. Viswanath (2008) arXiv:0809.1498
%-----------------------------------------------------------------
 function [y,del] = GMREShook(j,h,m,beta,del)

   a = h(1:j+1,1:j) ;

   [u,s,v] = svd(a) ; 
   s = diag(s) ;
   
   p = beta * u(1,1:j).' ;

   mu = max(s(j)*s(j)*1d-6,1d-99) ;
   qn = 1d99 ;
   while qn > del
      mu = mu * 1.1d0 ;
      q = p.*s ./ (mu+s.*s) ;
      qn = sqrt(sum(q.*q)) ;
   end

   y = v*q ;

   p = - h(1:j+1,1:j)*y ;
   p(1) = p(1) + beta ;
   del = sqrt(sum(p.*p)) ;
 
 end
