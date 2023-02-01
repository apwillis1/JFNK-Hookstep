# subroutine NewtonHook(...)

See also the newton-Lorenz template and comments within the Main file.

### Usage notes.

- In the Fortran version, prior to calling the subroutine `newtonhook(...)`, the arrays `new_x(:)` $\equiv{\bf x}_i$ 
and `new_fx(:)` $\equiv{\bf F}({\bf x}_i)$ need to be allocated with the dimension $n$:  `allocate(new_x(n),new_fx(n))`  
- Prior to calling `newtonhook(...)`, the array variable `new_x(:)` must be assigned the initial guess ${\bf x}_0$. 
As the calculation proceeds, `new_x(:)` is updated with the improved solution, having smaller corresponding `new_fx(:)` .
- `newtonhook(...)` also needs to be passed
   * a function that calculates the dot product of two vectors,
   * a function that calculates ${\bf F}({\bf x})$ for a given ${\bf x}$,
   * a function that calculates $\frac{{\bf dF}}{\bf {dx}}\ {\bf \delta x}$ at ${{\bf x}_i}$ for a given ${\bf \delta x}$.  See approximation in [README.md](./README.md), where ${\bf F}({\bf x}_i)$ is already stored in <tt>new_fx(:)</tt>.
   * a subroutine that replaces a vector ${\bf x}$ with the solution of $M{\bf x}'={\bf x}$ for ${\bf x}'$, where $M$ is a preconditioner matrix.  This may simply be an empty subroutine if no preconditioner is required, i.e. $M=I$.
   * a subroutine that is called at the end of each iteration.  This may be used to save the current state after each iteration, if desired.
- The functions above may require auxiliary data, in addition to the given ${\bf x}$ or ${\bf \delta x}$.  There are several possible solutions, e.g.:
   * place the required data in a module and access via `use mymodule`.
   * place the data in a <tt>common</tt> block.
   * place the data in `global` variables (MATLAB).
   * write the vector to disk, call a script that runs an external programme, then load the result.
- Extra constraints may be necessary, e.g. that determine an update to the period of an orbit.  
The function that evaluates $\frac{{\bf dF}}{{\bf dx}}\ {\bf \delta x}$ for a given ${\bf \delta x}$ should append to the 
result the evaluations of the constraints.  Correspondingly for each constraint, an extra <tt>0</tt> should be appended 
to ${\bf F}({\bf x}_i)$ (the evaluated constraint should equal zero when converged).  

### Parallel use

It is NOT necessary to edit this code for parallel (MPI) use:
* let each thread pass its subsection of the vector ${\bf x}$, 
* make the dot product  function <tt>mpi_allreduce</tt> the result of the dot product.
* to avoid multiple outputs to the terminal, set <tt>info=1</tt> on rank 0 and <tt>info=0</tt> for the other ranks.
