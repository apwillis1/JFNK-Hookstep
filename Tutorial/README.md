# Tutorial


This example seeks periodic solutions of the Lorenz equations, ${\phi}^T({\bf x})={\bf x}$.  We set up the system to be solved as ${\bf F}({\bf x})={\phi}^T({\bf x})-{\bf x}$.  Technically, since $T$ is also unknown, we append an extra condition to ${\bf F}$, but in the case of an equilibrium we could take $T$ fixed.  Please see https://doi.org/10.48550/arXiv.1908.06730 for details.


1. Download the MATLAB code (which looks surprisingly similar to the Fortran code).  The subroutines do the following:
- `Lorenz_f.m`:  Defines the Lorenz evolution rule $(\dot{X},\dot{Y},\dot{Z})={\bf f}(X,Y,Z)$.
- `steporbit.m`: Evaluate ${\phi}^T({\bf x})$, i.e. do `ndts_` timesteps according to ${\bf f}$, where the input vector is constructed as  x(1)=T and x(2:4)=(X,Y,Z).  The timestep size is dt=T/ndts\_.
- `saveorbit.m`: Called at the end of each iteration Newton iteration, this saves best ${\bf x}$ so far.  `relative_err` $=||{\bf F}({\bf x})|| / ||{\bf x}||$.
- `MAIN.m`:  Set up initial guess ${\bf x}_0$ and call the 'black box' `NewtonHook.m`.
