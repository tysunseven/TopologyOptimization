<script type="text/javascript" async
  src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js">
</script>
# Problem 4:  Method of Moving Asymptotes (MMA)
- [Problem 4:  Method of Moving Asymptotes (MMA)](#problem-4--method-of-moving-asymptotes-mma)
  - [MMA](#mma)
    - [Approximated optimization problem](#approximated-optimization-problem)
    - [The update strategy for (l\_j) and (u\_j)](#the-update-strategy-for-l_j-and-u_j)
    - [Relation to Common Forms of Optimization Problems](#relation-to-common-forms-of-optimization-problems)
    - [Solving the Approximated Optimization Problem](#solving-the-approximated-optimization-problem)
      - [Newton's method for finding the extremum of a function ( f )](#newtons-method-for-finding-the-extremum-of-a-function--f-)
      - [Slack variables method for optimization problem with inequality constraints](#slack-variables-method-for-optimization-problem-with-inequality-constraints)
      - [KKT conditions for optimization problem with inequality constraints](#kkt-conditions-for-optimization-problem-with-inequality-constraints)
      - [KKT conditions for Slack variables method](#kkt-conditions-for-slack-variables-method)
      - [Perturbed KKT/Primal-dual interior-point method](#perturbed-kktprimal-dual-interior-point-method)
    - [4. MMA的优点](#4-mma的优点)
    - [6. Summary](#6-summary)

## MMA
The MMA considers the following optimization problem:

$$
\text{Minimize } \quad f_0(x) + a_0 \cdot z + \sum_{i=1}^{m} \left( c_i \cdot y_i + \frac{d_i}{2} \cdot y_i^2 \right)
$$

$$
\begin{aligned}
\text{Subject to} \quad & f_i(x) - a_i \cdot z - y_i \leq 0,            & i = 1, 2, \ldots, m \\
                        & x_j^{\text{min}} \leq x_j \leq x_j^{\text{max}}, & j = 1, 2, \ldots, n \\
                        & z \geq 0,                                        &  \\
                        & y_i \geq 0,                                       & i = 1, 2, \ldots, m
\end{aligned}
$$

With appropriate choices of the coefficients \\( a_i \\), \\( c_i \\), and \\( d_i \\), it can be shown that this formulation is equivalent to other common forms of optimization problems. The details of this equivalence will be presented in subsequent sections. The overall framework of MMA is similar to the OC method, where an initial guess is made, and based on this guess, an approximate optimization problem is formulated. This approximate problem is then solved, and the solution is fed back into the main program. If the obtained solution satisfies the convergence criteria or meets the predefined requirements, the optimization loop is terminated. However, if the solution does not satisfy these conditions, it is used as the new starting point for the next iteration of the MMA process, which continues until a satisfactory solution is found or the maximum number of iterations is reached. Next, we will first describe what the approximated optimization problem looks like for a given guess, and then we will demonstrate how to solve the approximated problem using a primal-dual Newton method.
```matlab
function [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2, ...
f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);
%   m    = The number of general constraints.
%   n    = The number of variables x_j.
%  xmin  = Column vector with the lower bounds for the variables x_j.
%  xmax  = Column vector with the upper bounds for the variables x_j.
%  a0    = The constants a_0 in the term a_0*z.
%  a     = Column vector with the constants a_i in the terms a_i*z.
%  c     = Column vector with the constants c_i in the terms c_i*y_i.
%  d     = Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
```
After the value of \\( x \\) is updated, \\( x_{\min} \\) and \\( x_{\max} \\) are updated to the minimum value of \\( x \\) minus a certain value and the maximum value of \\( x \\) plus a certain value. Of course, these values cannot exceed the range of 0 to 1.
### Approximated optimization problem
Assuming that the guess at the \\(k\\)th iteration is \\(x^{(k)}\\) (with the initial guess manually specified and the \\(k\\)-th guess being the solution obtained from the previous iteration), MMA employs a method called *fractional approximation* to approximate the original nonlinear objective function and constraints into the following form:

$$
f_i^{(k)}(x) = \sum_{j=1}^{n} \left( \frac{p_{ij}^{(k)}}{u_j^{(k)} - x_j} + \frac{q_{ij}^{(k)}}{x_j - l_j^{(k)}} \right)+r_i^{(k)},\quad 0\leqslant i\leqslant m
$$

The coefficients \\( u_j^{(k)} \\) and \\( l_j^{(k)} \\), referred to as asymptotes (which is where the name "Method of Moving Asymptotes" originates, as these values are adjusted in each iteration), are determined in a somewhat ad hoc manner. The specific method for selecting these asymptotes will be discussed in detail later. After \\( u_j^{(k)} \\) and \\( l_j^{(k)} \\) are selected, the values of \\( p_{ij}^{(k)} \\), \\( q_{ij}^{(k)} \\), and \\( r_i^{(k)} \\) are determined such that \\( f_i^{(k)}(x) \\) and \\( f_i(x) \\) have the same value and first derivative at \\( x^{(k)} \\).
$$
p_{ij}^{(k)} = \left(u_j^{(k)} - x_j^{(k)}\right)^2  \frac{\partial f_i}{\partial x_j}^+(x^{(k)}),\quad q_{ij}^{(k)} = \left(x_j^{(k)} - l_j^{(k)}\right)^2  \frac{\partial f_i}{\partial x_j}^-(x^{(k)}).
$$

It should be noted that at most one of \\( p_{ij}^{(k)} \\) and \\( q_{ij}^{(k)} \\) can be non-zero. To invoke MMA, we need to calculate the values of \\(\frac{\partial f_i}{\partial x_j}(x^{(k)})\\) in the main program and pass them as parameters. In addition, the current iteration number \\(k\\) and the current guess \\(x^{(k)}\\) also need to be passed as parameters.
```matlab
%  iter  = Current iteration number ( =1 the first time mmasub is called).
%  xval  = Column vector with the current values of the variables x_j.
%  df0dx = Column vector with the derivatives of the objective function
%          f_0 with respect to the variables x_j, calculated at xval.
%  dfdx  = (m x n)-matrix with the derivatives of the constraint functions
%          f_i with respect to the variables x_j, calculated at xval.
%          dfdx(i,j) = the derivative of f_i with respect to x_j.
```

Based on the provided code, we can step-by-step derive how the expressions for \\( p_{ij}^{(k)} \\) and \\( q_{ij}^{(k)} \\) are constructed. 
1. **Defining Variables**:
   ```matlab
   ux1 = upp-xval;  ux2 = ux1.*ux1;  xl1 = xval-low;  xl2 = xl1.*xl1;
   ```
   $$ux1 = u_j^{(k)} - x_j^{(k)},\quad ux2 = (u_j^{(k)} - x_j^{(k)})^2,\quad xl1 = x_j^{(k)} - l_j^{(k)},\quad xl2 = (x_j^{(k)} - l_j^{(k)})^2  $$

2. **Computing the Positive and Negative Parts of the Derivatives**:
   ```matlab
   p0 = max(df0dx,0);  q0 = max(-df0dx,0);
   ```
$$p0_j = \frac{\partial f_0}{\partial x_j}^+(x^{(k)}),\quad q0_j = \frac{\partial f_0}{\partial x_j}^-(x^{(k)}) $$

3. **Adjustment Term**:
   ```matlab
   pq0 = 0.001*(p0 + q0) + raa0*xmamiinv;
   p0 = p0 + pq0;  q0 = q0 + pq0;
   ```

   Here, an adjustment term \\( pq0 \\) is added to ensure numerical stability. Although this part is not explicitly shown in the final expressions, it plays an important role in numerical computations. Where
   ```matlab
   raa0 = 0.00001;
   xmamiinv = eeen./xmami;

   eeen = ones(n,1);
   xmami = xmax-xmin;
   xmamieps = 0.00001*eeen;
   xmami = max(xmami,xmamieps);
   ```

    It can be seen that if \\(x_{\text{max}} - x_{\text{min}}\\) is less than 0.00001, the final result of \\(raa0 \times xmamiinv\\) is 1. If \\(x_{\text{max}} - x_{\text{min}}\\) is greater than 0.00001, the final result of \\(raa0 \times xmamiinv\\) is a value less than 1, and how much smaller it is depends on the scale of \\(x_{\text{min}}\\) relative to 0.00001. So in my opinion, when \\(x_{\text{max}} - x_{\text{min}}\\) is close to 0.00001, the value of \\(pq0\\) cannot be ignored. So the role of the \\( pq0 \\) term is still somewhat confusing.

4. **Combining the Upper/Lower Bounds with the Derivatives**:
   ```matlab
   p0 = p0.*ux2;  q0 = q0.*xl2;
   ```
   $$
p_{ij}^{(k)} = \left(u_j^{(k)} - x_j^{(k)}\right)^2  \frac{\partial f_i}{\partial x_j}^+(x^{(k)}),\quad q_{ij}^{(k)} = \left(x_j^{(k)} - l_j^{(k)}\right)^2  \frac{\partial f_i}{\partial x_j}^-(x^{(k)}).
$$

Earlier, we discussed how the objective function is approximated. Now, we will discuss how the range of the decision variable \\(x\\) is determined in the approximated optimization problem. It is clear that we need to satisfy $$ l_j^{(k)} < \alpha_j^{(k)} < x_j^{(k)} < \beta_j^{(k)} < u_j^{(k)} $$In practice, we choose
$$
\alpha_j^{(k)} = \max\left\{x_j^{\text{min}}, (1-\text{albefa})\times l_j^{(k)} + \text{albefa}\times x_j^{(k)},x_j^{(k)}-\text{move}\times(x^{\text{max}}_j-x^{\text{min}}_j)\right\}
$$

$$
\beta_j^{(k)} = \min\left\{x_j^{\text{max}}, (1-\text{albefa})\times u_j^{(k)} + \text{albefa}\times x_j^{(k)},x_j^{(k)}+\text{move}\times(x^{\text{max}}_j-x^{\text{min}}_j)\right\}
$$Where in mmasub.m,
```matlab
move = 1.0;  albefa = 0.1;
```


### The update strategy for \\(l_j\\) and \\(u_j\\)
In optimization algorithms, the magnitude of the second derivatives typically reflects the curvature of the objective or constraint functions. In the Method of Moving Asymptotes (MMA), the selection of the asymptotes \\( l_j^{(k)} \\) and \\( u_j^{(k)} \\) directly influences the magnitude of the second derivatives of the approximated functions, thus affecting the curvature and convergence properties of the optimization problem.


For any point \\( x \\) such that \\( l_j^{(k)} < x_j < u_j^{(k)} \\), the second derivative of the approximated function \\( f_i^{(k)} \\) is given by:

$$
\frac{\partial^2 f_i^{(k)}}{\partial x_j^2} = \frac{2p_{ij}^{(k)}}{(u_j^{(k)} - x_j)^3} + \frac{2q_{ij}^{(k)}}{(x_j - l_j^{(k)})^3},\quad \frac{\partial^2 f_i^{(k)}}{\partial x_j \partial x_l} = 0,\quad j\neq l
$$

Since \\( p_{ij}^{(k)} \geq 0 \\) and \\( q_{ij}^{(k)} \geq 0 \\), \\( f_i^{(k)} \\) is a convex function. Specifically, at \\( x = x^{(k)} \\), the second derivative simplifies to:

$$
\frac{\partial^2 f_i^{(k)}}{\partial x_j^2} = 2\frac{\partial f_i}{\partial x_j}\frac{1}{u_j^{(k)} - x_j^{(k)}},\quad  \text{if } \frac{\partial f_i}{\partial x_j} > 0, \quad \frac{\partial^2 f_i^{(k)}}{\partial x_j^2}=-2\frac{\partial f_i}{\partial x_j}\frac{1}{x_j^{(k)} - l_j^{(k)}},  \quad \text{if } \frac{\partial f_i}{\partial x_j} < 0
$$

Next, we will test a simple univariate function \\( f(x) \\) and \\( x_0 \\) to observe the effects of \\( u \\) and \\( l \\). It should first be noted that at most one of \\( p_{ij}^{(k)} \\) and \\( q_{ij}^{(k)} \\) can be non-zero at the same time. When the derivative of \\( f \\) at \\( x_0 \\) is greater than zero, we observe the term involving \\( u \\), and when the derivative of \\( f \\) at \\( x_0 \\) is less than zero, we observe the term involving \\( l \\). 
$$f''(x_0) \Delta x = -f'(x_0).$$

Since the approximation function shares the same value and first derivative as the original function, the Newton's method update formula above indicates that the closer the active term in \\( u \\) or \\( l \\) is to \\( x \\), the larger the value of \\( f''(x_0) \\) becomes. Consequently, the step size of the update becomes shorter, making our strategy more conservative. If the function behavior near \\( x_0 \\) is complex, an overly large step size might cause us to move directly from a monotonically decreasing interval to a monotonically increasing interval, potentially missing a local minimum point. In such cases, we should adjust our update strategy by reducing the step size, which corresponds to increasing the distance of \\( u \\) and \\( l \\) relative to \\( x_0 \\). In the code, one way to identify a complex function behavior near \\( x_0 \\) is by comparing the trend of the estimated solution \\( x \\) across two consecutive iterations. If the trend is reversed, it indicates that we have moved from an increasing/decreasing interval to a decreasing/increasing interval. Conversely, if the trend of the estimated solution \\( x \\) remains the same across consecutive iterations, it is reasonable to assume that we are in a relatively safe monotonic interval, allowing us to take a larger step size, which corresponds to decreasing the distance of \\( u \\) and \\( l \\) relative to \\( x_0 \\).

With the above analysis, it becomes easy to understand the update rules for \\( u \\) and \\( l \\) in the code. For the initial iterations (e.g., \\( k = 0 \\) and \\( k = 1 \\)), a simple choice for the asymptotes is given by:
$$
l_j^{(k)} = x_j^{(k)} - \text{asyinit} \times (x_{\max, j} - x_{\min, j}) \quad \text{and} \quad u_j^{(k)} = x_j^{(k)} + \text{asyinit} \times (x_{\max, j} - x_{\min, j})
$$

where asyinit is a fixed real number such as 0.01.
1. **If the process tends to oscillate**: When the variable \\( x_j \\) shows signs of oscillation between iterations \\( k-1 \\) and \\( k-2 \\), the asymptotes need to be stabilized by moving the asymptotes closer to the current iteration point:
$$ l_j^{(k)}=x_j^{(k)} - \text{asydecr} \times (x_j^{(k-1)} - l_j^{(k-1)}),\quad u_j^{(k)}=x_j^{(k)} + \text{asydecr} \times (x_j^{(k-1)} - u_j^{(k-1)}) $$

where asydecr is a fixed real number such as 0.7. We also set upper and lower bounds for the changes in \\( l \\) and \\( u \\). This is clearly evident in the code and will not be elaborated further here.

3. **If the process is monotone and slow**: In cases where the process is converging slowly and monotonically, the asymptotes should be relaxed by moving them away from the current iteration point:
$$ l_j^{(k)}=x_j^{(k)} - \text{asyincr} \times (x_j^{(k-1)} - l_j^{(k-1)}),\quad u_j^{(k)}=x_j^{(k)} + \text{asyincr} \times (x_j^{(k-1)} - u_j^{(k-1)}) $$

where asyincr is a fixed real number such as 1.2.
### Relation to Common Forms of Optimization Problems
To match an MMA problem to the standard NLP form, consider how the parameters \\( a_i \\), \\( c_i \\), and \\( d_i \\) are chosen:

1. **Bounds on Decision Variables**: In the MMA formulation, the bounds on \\( x_j \\) are directly analogous to the bounds \\( x_j^{\text{min}} \leq x_j \leq x_j^{\text{max}} \\) in the NLP form. These bounds ensure that the solution remains within a feasible region.

2. **Inequality Constraints**: The constraints \\( f_i(x) \leq 0 \\) in the NLP form correspond to the inequality constraints in MMA. To ensure that these constraints are properly integrated into the NLP framework, the parameters \\( a_i \\), \\( c_i \\), and \\( d_i \\) are chosen in a way that reflects the specific characteristics of the MMA model.

   - **Setting \\( a_i = 0 \\)**: This ensures that the variable \\( z \\) is set to zero in the optimal solution of the MMA problem, simplifying the problem and aligning it with the NLP form.
   - **Making \\( y_i \\) "expensive"**: By setting \\( d_i = 0 \\) and choosing \\( c_i \\) to be a large number, the variable \\( y_i \\) becomes costly to include in the solution, effectively penalizing its use. This aligns the MMA problem with the NLP form by minimizing unnecessary variables.

3. **Correspondence to Optimal Solutions**: When these adjustments are made, the resulting MMA problem matches the NLP form, allowing the optimal solution of the MMA problem to correspond directly to the optimal solution of the NLP problem. This ensures that the problem is not only tractable but also optimally aligned with the structure and requirements of standard NLP formulations.



### Solving the Approximated Optimization Problem
The benefit of using fractional approximation is that it transforms the problem into a convex form, making it easier to solve using standard optimization techniques.

#### Newton's method for finding the extremum of a function \\( f \\)
Given an initial guess \\( x_0 \\), perform a second-order Taylor expansion of \\( f \\) at \\( x_0 \\) and ignore higher-order terms.
$$
f(x_0 + \Delta x) \approx f(x_0) + \nabla f(x_0)^T \Delta x + \frac{1}{2} \Delta x^T H(x_k) \Delta x.
$$

We want to find the value of \(\Delta x\) such that \(\nabla f(x + \Delta x)\) is zero. Therefore, we take the gradient of the above expression with respect to \(\Delta x\) and set the left-hand side to zero, resulting in
$$
\nabla f(x_0) + H(x_0) \Delta x = 0\Longleftrightarrow H(x_0) \Delta x = -\nabla f(x_0).
$$

Then we set \\(x_1 = x_0 + \Delta x\\) and proceed to the next iteration.
#### Slack variables method for optimization problem with inequality constraints
Consider an optimization problem with inequality constraints
$$
\begin{aligned}
\text{Minimize } \quad & f(x) \\
\text{Subject to} \quad & g_i(x) \leq 0, & i = 1, \ldots, m \\
                        & h_j(x) = 0, & j = 1, \ldots, p
\end{aligned}
$$

By introducing slack variables \\( s \\) as decision variables, where \\( s_i \geq 0 \\), the original optimization problem can be transformed into an equivalent optimization problem as follows
$$
\begin{aligned}
\text{Minimize } \quad & f(x) \\
\text{Subject to} \quad & s_i \geq 0, & i = 1, \ldots, m \\
                        & g_i(x) + s_i = 0, & i = 1, \ldots, m \\
                        & h_j(x) = 0, & j = 1, \ldots, p
\end{aligned}
$$

It can be observed that an inequality has been transformed into an equality plus an inequality, but the new inequality is a direct constraint on the decision variable \\( s \\), which is somewhat different in role from the original inequality.

#### KKT conditions for optimization problem with inequality constraints
Consider an optimization problem with inequality constraints
$$
\begin{aligned}
\text{Minimize } \quad & f(x) \\
\text{Subject to} \quad & g_i(x) \leq 0, & i = 1, \ldots, m \\
                        & h_j(x) = 0, & j = 1, \ldots, p
\end{aligned}
$$

The corresponding Lagrangian function is
$$
L(x, u, v) = f(x) + \sum_{i=1}^m u_i g_i(x) + \sum_{j=1}^p v_j h_j(x).
$$

The **Karush-Kuhn-Tucker conditions** or **KKT conditions** are:

- **Stationarity**:
  $$
  \frac{\partial L}{\partial x}=\frac{\partial f}{\partial x}(x)+\sum_{i=1}^mu_i\frac{\partial g_i}{\partial x}(x)+\sum_{j=1}^pv_j\frac{\partial h_j}{\partial x}(x)=0.
  $$

- **Complementary Slackness**:
  $$
  u_i \cdot g_i(x) = 0 \quad \text{for all } i
  $$

- **Primal Feasibility**:
  $$
  g_i(x) \leq 0, \quad h_j(x) = 0 \quad \text{for all } i, j
  $$

- **Dual Feasibility**:
  $$
  u_i \geq 0 \quad \text{for all } i
  $$

#### KKT conditions for Slack variables method
As an example, we calculate the KKT conditions for 
$$
\begin{aligned}
\text{Minimize } \quad & f(x) \\
\text{Subject to} \quad & s_i \geq 0, & i = 1, \ldots, m \\
                        & g_i(x) + s_i = 0, & i = 1, \ldots, m \\
                        & h_j(x) = 0, & j = 1, \ldots, p
\end{aligned}
$$
The corresponding Lagrangian function is
$$
L(x, s, u, v, \lambda) = f(x) + \sum_{i=1}^m u_i (g_i(x)+s_i) + \sum_{j=1}^p v_j h_j(x)-\sum_{i=1}^m\lambda_i s_i.
$$
- **Stationarity**:
  $$
  \frac{\partial L}{\partial x}=\frac{\partial f}{\partial x}(x)+\sum_{i=1}^mu_i\frac{\partial g_i}{\partial x}(x)+\sum_{j=1}^pv_j\frac{\partial h_j}{\partial x}(x)=0,\quad \frac{\partial L}{\partial s_i}=u_i-\lambda_i=0.
  $$

- **Complementary Slackness**:
  $$
  \lambda_is_i = 0 \quad \text{for all } i
  $$

- **Primal Feasibility**:
  $$
  -s_i \leq 0, \quad g_i(x)+s_i=0,\quad h_j(x) = 0 \quad \text{for all } i, j
  $$

- **Dual Feasibility**:
  $$
  \lambda_i \geq 0 \quad \text{for all } i
  $$

The most specific equation is \\(\lambda_i = \mu_i\\). Typically, in the KKT system of equations, we directly use this equation to reduce the number of variables.
#### Perturbed KKT/Primal-dual interior-point method
- Suppose we have an optimization problem with \\(n\\) decision variables, \\(m\\) inequality constraints, and \\(p\\) equality constraints. The Lagrangian function for this problem has \\(n+m+p\\) parameters. The KKT conditions include \\(p\\) equality constraints, \\(n\\) stationarity conditions, and \\(m\\) complementary slackness conditions, totaling \\(n+m+p\\) conditions. Additionally, we have \\(m\\) inequality constraints.

- If you decide to introduce \\(l\\) slack variables, the number of decision variables becomes \\(n+l\\), with \\(m\\) inequality constraints and \\(p+l\\) equality constraints. The Lagrangian function should now have \\(n+l+(m-l+l)+p+l\\) variables: \\(n\\) for \\(x\\), \\(l\\) for \\(s\\), \\(m-l+l\\) for the Lagrange multipliers of the inequality constraints \\((m-l\\) for those without slack variables and \\(l\\) for those with), \\(p\\) for the original equality constraints, and \\(l\\) for the equality constraints introduced by the slack variables. You will have \\(n+l\\) stationarity conditions, \\(m-l+l\\) complementary slackness conditions, and \\(p+l\\) equality constraints. By using \\(l\\) of the stationarity conditions, you equate the Lagrange multipliers of the inequality constraints with slack variables to those of the equality constraints introduced by the slack variables. In total, there are \\(n+l+m+p\\) variables and \\(n+l+m+p\\) equality constraints. Additionally, we still have \\(m\\) inequality constraints.

- Let’s denote all decision variables as \\(x\\) and all constraints as a vector-valued function \\(F\\). We aim to solve \\(F(x) = 0\\). Here, following the spirit of the Newton method, we first find a starting point \\(x_0\\) that satisfies the inequality constraints:
$$F(x_0 + \Delta x) = F(x_0) + \nabla F(x_0)\Delta x = 0$$This determines the update direction \\(\Delta x\\), after which an appropriate step length is taken to satisfy both the inequality constraints and ensure convergence.

- Since the starting point \\(x_0\\) is within the feasible region, and the iterated values remain within the feasible region, this method falls under the category of interior-point methods. Since the decision variables of the primal problem and the dual variables (i.e., the Lagrange multipliers) are updated simultaneously, it is called a primal-dual method.

- To avoid numerical instability or slow convergence due to the right-hand side of the complementary slackness conditions being zero, we introduce a perturbation \\(\epsilon\\) on the right-hand side of these equations. During the iteration process, \\(\epsilon\\) is gradually reduced. When it is reduced to a certain threshold, the solution \\(F(x) = 0\\) is considered to be found. This completes the solution of the subproblem in MMA, and the updated \\(x\\) values are returned to the main program for the next iteration. Note that there are two layers of iteration here.

With the above analysis, we can clearly understand the transition from 7.4(a)-7.4(n) to 7.5(a)-7.5(n).

- In 7.4, \\( d \\), \\( j \\), and \\( k \\) are inequality constraints, but \\( j \\) and \\( k \\) are constraints directly on the values of the decision variables, while \\( d \\) can be considered a more fundamental constraint. Therefore, we introduce slack variables \\( s \\) for the inequality constraints \\( d \\).
- It can be seen that in 7.4, the Lagrange multipliers corresponding to the constraints \\( d \\) are denoted by \\( \lambda \\). In 7.5, the corresponding equality constraints also use \\( \lambda \\), and due to the reasons explained earlier, the Lagrange multipliers for the inequality constraint \\( s > 0 \\) for the slack variables are also equal to \\( \lambda \\). This explains why the coefficient in front of \\( s \\) in 7.5(i) is \\( \lambda \\).
- In 7.5, \\( a \\)-\\( c \\) are the stationarity conditions for the original decision variables, \\( d \\) represents the equality constraints after introducing the slack variables, \\( j \\) to \\( m \\) represent the original inequality constraints and the inequality constraints on the dual variables, and \\( n \\) represents the new inequality constraints for the slack variables and their corresponding dual variables. \\( e \\) to \\( i \\) are the complementary slackness conditions for these inequality constraints, perturbed by \\( \epsilon \\).
- On page 15 of the paper, the Newton method is used to solve the equation \\( F(x) = 0 \\).

### 6. Summary

The MMA (Method of Moving Asymptotes) effectively addresses the complexity and numerical instability in large-scale nonlinear optimization problems by introducing the concept of moving asymptotes. In each iteration, it constructs a simplified, easy-to-solve convex subproblem and controls the convergence process by dynamically adjusting the asymptotes, ensuring the efficiency and stability of the algorithm. Due to these advantages, the MMA method has been widely applied in fields such as structural optimization and mechanical design.

