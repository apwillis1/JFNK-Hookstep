# Tutorial and Overview

Detailed notes can be found at https://doi.org/10.48550/arXiv.1908.06730 .  See Section 4.1 of the notes on adapting the code for your own use.  The tutorial of Section 4 seeks periodic solutions of the Lorenz equations.  Below is a brief overview of the tutorial and the code.

For an initial state ${\bf x}$, we say that ${\phi}^T({\bf x})$ is the result of time integration for a time $T$.  A periodic solution then satisfies ${\phi}^T({\bf x})={\bf x}$.  We set up the system to be solved as ${\bf F}({\bf x})={\phi}^T({\bf x})-{\bf x}$.  If $T$ is an unknown, we must append an extra condition to ${\bf F}$, but in the case of an equilibrium we can select any fixed value for $T$ (typically $O(1-10)$ for a non-dimensionalised system).  


The MATLAB and Fortran codes look surprisingly similar.  The subroutines do the following:
- `Lorenz_f.m`:  Defines the Lorenz evolution rule $(\dot{X},\dot{Y},\dot{Z})={\bf f}(X,Y,Z)$.
- `steporbit.m`: Evaluate ${\phi}^T({\bf x})$, i.e. do `ndts_` timesteps according to ${\bf f}$, where the input vector is constructed as  x(1)=T and x(2:4)=(X,Y,Z).  The timestep size is dt=T/ndts\_.
- `saveorbit.m`: Called at the end of each iteration Newton iteration, this saves best ${\bf x}$ so far.  `relative_err` $=||{\bf F}({\bf x})|| / ||{\bf x}||$.
- `MAIN.m`:  Set up initial guess ${\bf x}_0$ and call the 'black box' `NewtonHook.m`.

Other functions are called by `NewtonHook.m` and are unlikely to need changing for a problem of this type.
- `getrhs.m`: Evaluate right-hand side, i.e. ${\bf F}({\bf x})=\phi^T({\bf x})-{\bf x}$.  See (3.11) of the linked document.
- `multJ.m`: Evaluate multiplication by the Jacobian.  See (3.13) for the finite difference approximation.
- `multJp.m`: Preconditioner for multiplication, here an empty function.
- `dotprd.m`: Evaluate inner product $\langle\vec{a}|\vec{b}\rangle$.
- `GMRESm.m`: Generalized minimized residual method.  Section 3.2 of the notes.
- `GMREShook.m`: Calculate hookstep.  Section 3.3 of the notes.
