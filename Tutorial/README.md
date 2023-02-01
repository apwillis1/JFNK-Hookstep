# Tutorial and Adapting the code for your own use

**Detailed notes can be found at https://doi.org/10.48550/arXiv.1908.06730 .  The notes below are largely a repetition of Section 4.**

**The tutorial of Section 4 seeks periodic solutions of the Lorenz equations.**

For an initial state ${\bf x}$, we say that ${\phi}^T({\bf x})$ is the result of time integration for a time $T$.  A periodic solution then satisfies ${\phi}^T({\bf x})={\bf x}$.  We set up the system to be solved as ${\bf F}({\bf x})={\phi}^T({\bf x})-{\bf x}$.  If $T$ is an unknown, we must append an extra condition to ${\bf F}$, but in the case of an equilibrium we can select any fixed value for $T$ (typically $O(1-10)$ for a non-dimensionalised system).  

### Code overview

The MATLAB and Fortran codes look surprisingly similar.  The subroutines do the following:
- `Lorenz_f.m`:  Defines the Lorenz evolution rule $(\dot{X},\dot{Y},\dot{Z})={\bf f}(X,Y,Z)$.
- `steporbit.m`: Evaluate ${\phi}^T({\bf x})$, i.e. do `ndts_` timesteps according to ${\bf f}$, where the input vector is constructed as  x(1)=T and x(2:4)=(X,Y,Z).  The timestep size is dt=T/ndts\_.
- `saveorbit.m`: Called at the end of each iteration Newton iteration, this saves best ${\bf x}$ so far.  `relative_err` $=||{\bf F}({\bf x})|| / ||{\bf x}||$.
- `MAIN.m`:  Set up initial guess ${\bf x}_0$ and call the 'black box' `NewtonHook.m`.  [NewtonHook subroutine](https://openpipeflow.org/index.php?title=File:NewtonHook.f90)

Other functions are called by `NewtonHook.m` and are unlikely to need changing for a problem of this type.
- `getrhs.m`: Evaluate right-hand side, i.e. ${\bf F}({\bf x})=\phi^T({\bf x})-{\bf x}$.  See (3.11) of the linked document.
- `multJ.m`: Evaluate multiplication by the Jacobian.  See (3.13) for the finite difference approximation.
- `multJp.m`: Preconditioner for multiplication, here an empty function.
- `dotprd.m`: Evaluate inner product $\langle{\bf a}|{\bf b}\rangle$.
- `GMRESm.m`: Generalized minimized residual method.  Section 3.2 of the notes.  [GMRES subroutine](https://openpipeflow.org/index.php?title=File:GMRESm.f90).

- `GMREShook.m`: Calculate hookstep.  Section 3.3 of the notes.


## Tutorial

1. The following data are points on the periodic orbits of figure 8 of the notes, taken from Viswanath (2003).  $Z=27$ in call cases.
``` 
Orbit      X                  Y                T
   AB   −13.763610682134   −19.578751942452   1.5586522107162
  AAB   −12.595115397689   −16.970525307084   2.3059072639399
 AAAB   −11.998523280062   −15.684254096883   3.0235837034339
 AABB   −12.915137970311   −17.673100172646   3.0842767758221
```
2. In MATLAB, call `MAIN`.  It will plot the result of timestepping the initial guess for the 
AB orbit (green), call the `NewtonHook` subroutine, then plot the converged solution (blue).
Scroll back through the output, and 
compare `relative_err` for the initial guess at iteration 0 with the final relative error. 
3. Comment/uncomment other initial guesses `new_x` $={\bf x}_0$, or experiment with your own.  How do they affect the number of 
Newton iterations taken? Typically convergence takes $O(10)$ iterations, otherwise it will never converge.
4. Uncomment the initial guess for an equilibrium.  Here we assume a short fixed $T$, too short for a PO; $T$ is not permitted to change, otherwise 
$||{\bf \phi}^T({\bf x})-{\bf x}||$ could be reduced by simply taking $T\to 0$.
Check that `MAIN` can find the analytic equilibrium solution $(\pm\alpha,\pm\alpha,r-1)$, 
where $\alpha=\sqrt{(r-1)\,b}$. 

## Adapting the code for your own use

- For a very large system, for which you might consider parallelization 
(see final comment), you should probably use the Fortran90 version.
- Experiment with the Template/Example first, to get used to how the code is set up.  The initial guess is put in `new_x`.
- Note that at present, `new_x(1)` $=T$ (the period), and `new_x(2:end}` $={\bf x}$ (the state).
- The place to start editing code is `steporbit`.  If you already have 
an existing timestepping code, it could do something as simple 
as call it externally via system calls:
```matlab
 function y = steporbit(ndts_,x)
   persistent dt

   if ndts_ ~= 1                % Set timestep size dt=T/ndts_
      dt = x(1) / ndts_ ;       % If only doing one step to calc \dot{x},
   end                          % then use previously set dt.

   a = x(2:end) ;

   WRITE DATA TO FILES:
      dt     timestep size
      ndts_  number of steps to take
      a      initial condition
   
   LOAD STATE, TIMESTEP, SAVE STATE:
      system('run_my_code.exe')   
   
   LOAD TIMESTEPPED STATE: --> a

   y = zeros(size(x)) ;
   y(2:end) = a ;
 end 
```
- `saveorbit` is called at the end of each Newton iteration.  Add code here to save the current state `new_x`.
- If your inner product corresponds to $\langle{\bf a}|{\bf b}\rangle\ =\ {\bf a}^TW{\bf b}$ 
where $W$ is a diagonal matrix of positive weights, and here $T$ is the 
transpose, then pass ${\bf x}'=W^\frac{1}{2}{\bf x}$ to the code.
The existing functions that take inner products then need no modification.
- For parallel use with MPI+Fortran, the `NewtonHook` and `GMRES` codes do not need changing:
Split vectors over threads and let each thread pass its section
to `NewtonHook`.  The only place where an MPI call is required is an MPI\_Allreduce in the 
`dotprod` function.  To avoid all threads outputting information, set `info=1` on rank 0, and `info=0` on all other ranks.
- Further usage notes can be found at 
[NewtonHook subroutine](https://openpipeflow.org/index.php?title=File:NewtonHook.f90) and 
[GMRES subroutine](https://openpipeflow.org/index.php?title=File:GMRESm.f90).

