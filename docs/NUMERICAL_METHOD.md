# Numerical method and implementation conventions

## Model

For a power-law exponent $p>1$, the implementation uses

$$
S(A)=\nu_0(\delta+|A^{\mathrm{sym}}|)^{p-2}A^{\mathrm{sym}},
\qquad
F(A)=(\delta+|A^{\mathrm{sym}}|)^{(p-2)/2}A^{\mathrm{sym}}.
$$

The experiment is posed on $\Omega=(0,1)^2$. The manufactured velocity is
divergence free but does not have a homogeneous normal trace on the square.
Consequently, the implemented boundary conditions are the inhomogeneous
extension

$$
u\cdot n=g_n, \qquad (S(Du)n)_\tau=g_\tau.
$$

The data $g_n$ and $g_\tau$ are obtained from the exact manufactured
velocity and stress.

## Spatial discretization

Both variants use the lowest-order Taylor-Hood pair

$$
\widehat V_h=(P_2^c(\mathcal T_h))^2,
\qquad \widehat Q_h=P_1^c(\mathcal T_h),
$$

and one global Real element to enforce zero mean pressure.

### Strong normal condition

On the axis-aligned square, $u_h\cdot n=g_n$ is imposed by prescribing the
first velocity component on the vertical sides and the second component on
the horizontal sides. Tangential components remain unconstrained.

This implementation is intentionally geometry specific. Extending it to a
general polygon requires a normal-trace-aware method; the multiplier variant
already provides one.

### Lifted multiplier normal condition

The code unknown is a lifted multiplier represented by the normal trace of

```python
FiniteElement("RT", triangle, 2)["facet"]
```

whose facet normal trace has the required piecewise-linear order. The residual
contains

$$
(u_h\cdot n,m_h\cdot n)_\Gamma
-(\ell_h\cdot n,v_h\cdot n)_\Gamma
-(g_n,m_h\cdot n)_\Gamma.
$$

The complete known stress vector $S(Du_{\mathrm{ex}})n$, rather than only its
tangential component, is placed on the right-hand side. Decomposing this load
into tangential and normal parts gives

$$
\lambda_{\mathrm{article}}
=-(\ell_h\cdot n)-n\cdot S(Du_{\mathrm{ex}})n.
$$

Thus $\ell_h\cdot n$ is not itself the physical normal-stress multiplier.
For the exact manufactured solution it equals $-q$. This lifting is useful
because the additional unknown then has the same explicitly controlled
regularity as the pressure. The helper
`pstokes_fem.boundary.article_multiplier` performs the reconstruction whenever
the physical multiplier is needed.

```python
lifted_multiplier = solution.split(deepcopy=True)[3]
stress_load = manufactured.averaged_stress(t_left, t_right)
lambda_article = article_multiplier(lifted_multiplier, stress_load, mesh)
```

The restricted RT element also has facet degrees of freedom on interior
facets. They are fixed by

$$
(\ell_h^+\cdot n^+,m_h^+\cdot n^+)_{\mathcal S_h^{\mathrm{int}}}.
$$

This term only extends the boundary multiplier by zero into the interior; it
does not alter the velocity-pressure problem. The rank test in
`tests/test_multiplier_nullspace.py` verifies the corresponding trace mass
matrix.

## Time discretization

For $t_m=m\tau$, the target residual uses backward Euler,

$$
\left(\frac{u_h^m-u_h^{m-1}}{\tau},v_h\right)
+(S(Du_h^m),Dv_h)-(p_h^m,\operatorname{div}v_h).
$$

The default load rule `right` evaluates the volume force and natural traction
at $t_m$, reproducing the supplied code. The rules `midpoint`, `gauss2`, and
`gauss3` approximate their interval averages. Essential normal data are always
interpolated at $t_m$.

## Nonlinear continuation

Continuation solves are used only as Newton initializers:

- strong variant: $p=2$, then $p=1.75$, then the target exponent;
- multiplier variant: $p=2$, then the target exponent.

Only the final stage is recorded. Therefore, continuation changes solver
robustness but not the target discrete equations.

## Manufactured pressure

The pressure has zero spatial mean and uses the article's scaling

$$
c_q=10^{-3}\quad(p=1.5),
\qquad
c_q=10^3\quad(p=2.5).
$$

Other exponents require an explicit `--pressure-amplitude` argument.

## Discrete Leray projector experiment

Section 6.3 uses the meshes obtained by dividing the unit square into
`i * i` squares and splitting each square diagonally, for `i = 1, ..., 39`.
Thus `h_i = sqrt(2) / i` for DOLFIN's `UnitSquareMesh(i, i)`.

For every nodal basis function `phi_j` of the unconstrained vector P2 space,
the implementation first computes its L2 projection `w_h = P_Vh phi_j` onto
the velocity space with homogeneous normal trace. In the multiplier variant,
this is the saddle-point projection with the same facet-restricted RT trace
space and interior zero lift used above.

To compute the complementary projector, the code solves for a potential
`z_h`, its discrete gradient `g_h`, and a real mean multiplier `r_h`:

$$
 (g_h,\nabla q_h)+(z_h,s_h)+(r_h,q_h)
 +(g_h,v_h)+(z_h,\operatorname{div}v_h)
 =-(\operatorname{div}w_h,q_h).
$$

The velocity test and trial spaces carry the homogeneous normal condition,
strongly or through the RT trace saddle point. Consequently,

$$
 \mathcal P_h^\perp w_h=g_h,
 \qquad
 \mathcal P_h w_h=w_h-g_h,
$$

which is the realization of
`grad_h Delta_h^{-1} div_h` used in the article. The reported value is

$$
 \max_j
 \frac{\|\mathcal J_h P_{V_h}\varphi_j\|_{L^r}}
      {\|P_{V_h}\varphi_j\|_{L^r}},
 \qquad r\in\{2,p,p'\}.
$$

This maximum over the projected nodal basis is a stability indicator, not the
induced Lr operator norm. Quadrature degree 6 reproduces the prototype; higher
orders can be evaluated in the same basis traversal.

## Reported errors

The velocity quantity is intentionally

$$
\max_m\|u_h^m-u(t_m)\|_{L^2(\Omega)}
+\left(\tau\sum_m
\|F(Du_h^m)-F(Du(t_m))\|_{L^2(\Omega)}^2\right)^{1/2}.
$$

The $L^\infty(I;L^2(\Omega))$ term is the quantity used to produce the
published tables. An $L^2(I;L^2(\Omega))$ label would be a typographical
error for these data.
