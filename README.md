# JFNK-Hookstep
### Jacobian-Free Newton-Krylov solver with Hookstep-trust-region approach.

This is a solver for systems in the form ${\bf F}({\bf x})={\bf 0}$, where the dimension (length) $n$ of the vector ${\bf x}$ can be either large or small.  The trust-region approach enlarges the basin of attraction of a solution ${\bf x}_0$ by automatically adjusting the size of each step in the Newton iteration, while the Hookstep optimises the direction of the step to accelerate convergence.  

The method is 'Jacobian-Free' in the sense that the user only need supply a function that evaluates ${\bf F}({\bf x})$ for given ${\bf x}$, not the Jacobian matrix itself, despite that the Jacobian appears in the Newton iteration formula.  The Newton iteration formula is solved via the GMRES(m) algorithm, the implementation of which may be supplied a preconditioner.

I hope that you will find that the code is simply written, and that it can be bolted on to any existing solver for ${\bf F}({\bf x})$, e.g. ${\bf F}$ is the result of timestepping $\{\bf x}$.  The template solves for periodic orbits of the Lorenz equations, $n=3+1$ (dimension of the system + unknown period), while the same code has been used without modification to compute nonlinear equilibria of pipe flow, $n=O(10^6)$, via parallel (MPI) simulations.

The JFNK solver and GMRES codes could be run 'native' (MATLAB/Fortran90) or could be integrated with codes developed in other languages by asking ${\bf F}$ to save ${\bf x}$ to disk, execute the existing code via a terminal command, then loading the result.  The MATLAB version will run under the free alternative Octave.

Developed as part of www.openpipeflow.org

FURTHER DETAILS ON THE METHOD:  https://doi.org/10.48550/arXiv.1908.06730

CITATION:  https://doi.org/10.48550/arXiv.1908.06730 or https://doi.org/10.1016/j.softx.2017.05.003

AUTHOR:  Ashley P. Willis, School of Mathematics and Statistics, University of Sheffield.  https://www.sheffield.ac.uk/maths/research/fluid

THANKS:  Rich Kerswell, Predrag Cvitanovi\'c ([chaosbook.org](http://www.chaosbook.org)), John Gibson ([channelflow.org](http://www.channelflow.org)), 
Marc Avila and many others for their generous support in many forms.  Developed under EPSRC grants EP/K03636X/1, EP/P000959/1.


## Download

With git:
```
$ git clone https://github.com/apwillis1/JFNK-Hookstep.git
```
Manually: 
1. click on 'MATLAB' or 'Fortran90' above, 
2. click on a code file, 
3. click on 'Raw', 
4. Ctrl+S, or use browser option 'Save page as...' 

## Using the code

[**Tutorial** and **Adapting the code for your own use**](Tutorial/README.md).

[NewtonHook subroutine](./NewtonHook.md).

[GMRES subroutine](./GMRESm.md).


## Method overview

The codes implement the 'Jacobian-free Newton-Krylov (JFNK) method' for solving 
$${\bf F}({\bf x})={\bf 0},$$ 
where ${\bf x}$ and ${\bf F}$ are $n$-vectors, supplemented with a Hookstep--Trust-region approach.

This is a powerful method that can solve for ${\bf x}$ for a complicated nonlinear ${\bf F}({\bf x})$.  For example, to find an ''equilibrium'' solution or a ''periodic orbit'', let ${\bf F}({\bf x}) = {\bf X}({\bf x})-{\bf x}$, where ${\bf X}({\bf x})$ is the result of time integration of an initial condition ${\bf x}$.

### Newton-Raphson method

To find the roots $x$ of a function $f(x)$ in one dimension, given an initial guess $x_0$, the Newton-Raphson method generates improvements using the iteration $x_{i+1}=x_i-f(x_i)/f'(x_i)$, where the dash denotes the derivative.  The iteration can be re-expressed as
$$x_{i+1}=x_i+\delta x_i \quad\mbox{where}\quad f'(x_i)\ \delta x_i = -f(x_i).$$
The extension of Newton's method to an $n$-dimensional system is then
$$(a)\quad {\bf x}_{i+1} = {\bf x}_i + {\bf \delta x}_i  \quad\mbox{where}\quad (b)\quad \frac{{\bf \partial F}}{{\bf \partial x}}\ {\bf \delta x}_i = -{\bf F}({\bf x}_i).$$
The $n\times n$ Jacobian matrix $\frac{{\bf \partial F}}{{\bf \partial x}}$ is evaluated at ${\bf x}_i$.  In order to apply the update $(a)$, the linear system $(b)$ needs to be solved for the unknown ${\bf \delta x}_i$.

### Jacobian-Free Newton method

The Jacobian matrix is usually difficult to evaluate.  However, the problem $(b)$ is in the form $A\ {\bf \delta x}={\bf b}$, which can be solved using the Krylov-subspace method GMRES(m).  The GMRES algorithm does not need to know the matrix $A$ itself, only the result of multiplying a vector by $A$.  For a given starting vector ${\bf \delta x}_0$, e.g. ${\bf \delta x}_0={\bf b}$, the method seeks solutions for ${\bf \delta x}$ in $\mathrm{span}({\bf \delta x}_0,\ A\ {\bf\delta x}_0,\ A^2\ {\bf\delta x}_0,...)$, but uses Gram-Schmidt orthogonalisation to improve the numerical suitability of this basis set.  The set of orthogonalised vectors is called the Krylov-subspace, and m is the maximum number of vectors stored.

Iterations of the GMRES algorithm for the problem $(b)$ involve calculating products of the Jacobian with given ${\bf \delta x}$, which may be approximated, e.g.
$$(c)\qquad \frac{{\bf \partial F}}{{\bf \partial x}} {\bf \delta x} \approx \frac{1}{\epsilon}({\bf F}({\bf x}_i+\epsilon\ {\bf \delta x})-{\bf F}({\bf x}_i))$$ 
for some small value $\epsilon$.  (Try $\epsilon$ such that $(\epsilon\ ||{\bf \delta x}||)\ /\ ||{\bf x}_i||=10^{-6}$.)  The important point is that we do not need to know the Jacobian; **only a routine for evaluating** ${\bf F}({\bf x})$ **is required**.

Note that provided that each step of the Newton method, ${\bf \delta x}$, is essentially in the correct direction, the method is expected to converge.  Therefore the tolerance specified in the accuracy of the solution for ${\bf \delta x}$ in each Newton step (calculated via the GMRES method) typically need not be so stringent as the tolerance placed on the Newton method itself for the solution ${\bf x}$.  For example, we might seek a relative error for the Newton solution $||{\bf F}({\bf x})||/||{\bf x}||=O(10^{-8})$, but a relative error for the GMRES steps $||A{\bf \delta x}-{\bf b}||/||{\bf \delta x}||=O(10^{-3})$ is likely to be sufficient for calculation of the steps ${\bf \delta x}$.

### Hookstep approach

To improve the domain of convergence of the Newton method, it is commonplace to limit the size of the step taken.  One approach is simply to take a 'damped' step in the direction of the solution to $(b)$, i.e. step by $\alpha\ {\bf \delta x}_i$ , where $\alpha \in (0,1]$.  

In the hookstep approach, we minimise subject to the condition that the magnitude of the Newton step is limited, $||{\bf \delta x}_i|| < \delta$ 
where $\delta$ is the size of the 'trust region'

$$(d)\quad \min_{{\bf \delta x}_i:\ ||{\bf \delta x}_i||<\delta} \  || {\bf F}({\bf x}_i) + \frac{{\bf \partial F}}{{\bf \partial x}} {\bf \delta x}_i || \ .$$
Given the minimisation, the hookstep ${\bf \delta x}_i$ is expected to produce a better result than a simple damped step of the same size.  It is also expected to perform much better in 'valleys', where it produces a bent/hooked step to a point along the valley, 
rather than jumping from one side of the valley to the other.

The hookstep can be calculated with little extra work to the GMRES method, provided that the size of Krylov-subspace, m, is chosen sufficiently large to solve to the desired accuracy within m GMRES iterations.

For a given ${\bf \delta x}_i$, the reduction in error predicted by the minimisation $(d)$ can be compared with the actual reduction in $||{\bf F}({\bf x})||$.  According to the accuracy of the prediction, the size of the trust region $\delta$ can be adjusted automatically.

### Preconditioning 

The GMRES implementations can be supplied with a preconditioner routine (see [GMRES subroutine](./GMRESm.md)), but is unlikely to be necessary when the method is combined with time integration.

### Further details

Extended overview https://doi.org/10.48550/arXiv.1908.06730 .  See section 4 for a [Tutorial and Adapting the code for your own use](Tutorial/README.md)

Further details on the method at [channelflow.org](http://channelflow.org/dokuwiki/doku.php?id=docs:math:newton_krylov_hookstep).
