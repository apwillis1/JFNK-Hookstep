# subroutine GMRESm(...)

This implements the classic GMRES(m) method for solving the system 
$$A{\bf x}={\bf b}$$ for ${\bf x}$.  This implementation minimises the error $||A{\bf x}-{\bf b}||$ subject to the additional constraint $||{\bf x}||<\delta$.  (This constraint may be ignored by supplying $\delta<0$ to the call.)  Available for Fortran90, MATLAB/Octave.

### The GMRES(m) Method 

The main advantage of the GMRES method is that it only requires calculations of multiplies by $A$ for a given ${\bf x}$, it does not need to know $A$ itself.  This means that $A$ need not even be stored, and could correspond to a very complex linear 'action' on ${\bf x}$, e.g. a time integral with initial condition ${\bf x}$.  For a given starting vector ${\bf x}_0$, the method seeks solutions for ${\bf x}$ in $\mathrm{span}({\bf x}_0,\ A{\bf x}_0,\ A^2{\bf x}_0,...)$, but uses Gram-Schmidt orthogonalisation to improve the numerical suitability of this basis set.  The set of orthogonalised vectors is called the Krylov-subspace, and m is the maximum number of vectors stored.

Whereas m is traditionally a small number, e.g. 3 or 4, the additional constraint renders restarts difficult.  If the constraint is important, then m must be chosen sufficiently large to solve to the desired accuracy within m iterations.

### Preconditioning

The implementation can be supplied a preconditioner routine. (This can be avoided if combined with timestepping; see the remarks at [Arnoldi](https://openpipeflow.org/index.php?title=File:GMRESm.f90).

GMRES is likely to find it easier to solve $M^{-1}A\,x=M^{-1}b$ than the original system, if $M^{-1}$ is an approximate inverse for $A$.
For example, if $A$ is dominated by its diagonal elements, we might take $M$ to be the banded matrix consisting of the diagonal and the first sub- and super-diagonals of $A$.  Each GMRES iteration applied to the modified system now requires a muliplication by $A$ then by $M^{-1}$.  This is fine, as it is quick and easy to solve $Mx'=x$ for $x'$ for a banded matrix $M$.  Like $A$, we don't need to know the matrix $M^{-1}$ itself, only the result of multiplication by these matrices.

## Using the code

- The Fortran implementation requires the LAPACK package.
- The constraint $||{\bf x}||<\delta$ may be ignored by supplying negative <tt>del</tt>.
- In addition to scalar and array variables, the routine needs to be passed
   * an external function that calculates dot products,
   * an external subroutine that calculates the result of multiplication by $A$,
   * an external subroutine that replaces a given vector ${\bf x}$ with the solution ${\bf x}'$ of the system $M{\bf x}'={\bf x}$.  This may simply be an empty subroutine if no preconditioner is required, i.e. $M=I$.
- The functions above may require auxiliary data in addition to ${\bf x}$ or ${\bf \delta x}$.  Place this data in a module and access via <tt>use mymodule</tt> in the function/subroutine, or (MATLAB) place in `global` variables.

### Parallel use

It is NOT necessary to edit this code for parallel (MPI) use:
  * let each thread pass its subsection for the vector ${\bf x}$, 
  * make the dot product function <tt>mpi_allreduce</tt> the result of the dot product.
  * to avoid multiple outputs to the terminal, set <tt>info=1</tt> on rank 0 and <tt>info=0</tt> for the other ranks.
