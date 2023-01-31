# JFNK-Hookstep
Jacobian-Free Newton-Krylov solver with Hookstep-trust-region approach.

This is a solver for systems in the form ${\bf F}({\bf x})={\bf 0}$, where the dimension (length) $n$ of the vector ${\bf x}$ can be either large or small.  The trust-region approach enlarges the basin of attraction of a solution ${\bf x}_0$ by automatically adjusting the size of each step in the Newton iteration, while the Hookstep optimises the direction of the step to accellerate convergence. 

I hope that you will find that the code is simply written, and that it can be bolted on to any existing solver for ${\bf F}({\bf x})$, e.g. ${\bf F}$ is the result of timestepping $\{\bf x}$.  The template solves for periodic orbits of the Lorenz equations, $n=3$, while the same code has been used to compute nonlinear equilibria of pipe flow, $n=O(10^5)$.

The JFNK solver and GMRES code are currently available in MATLAB and Fortran90.  The MATLAB version will run under the free alternative Octave.

CITATION AND FURTHER INFO:  https://doi.org/10.48550/arXiv.1908.06730

AUTHOR:  Ashley P. Willis, School of Mathematics and Statistics, University of Sheffield. 
