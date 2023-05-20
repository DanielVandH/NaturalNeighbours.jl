```@meta
CurrentModule = NaturalNeighbours
```

In this section, we give some of the mathematical background behind natural neighbour interpolation, and other interpolation methods provided from this package. The discussion here will be limited, but you can see this [thesis](https://kluedo.ub.rptu.de/frontdoor/deliver/index/docId/2104/file/diss.bobach.natural.neighbor.20090615.pdf) or this [Wikipedia article](https://en.wikipedia.org/wiki/Natural_neighbor_interpolation) for more information. The discussion that follows is primarily sourced from the linked thesis. These ideas are implemented by the `interpolate` function.

# Voronoi Tessellation 

We need to first define the _Voronoi tessellation_. We consider some set of points $\boldsymbol X = \{\boldsymbol x_1, \ldots, \boldsymbol x_m\} \subseteq \mathbb R^2$. The Voronoi tessellation of $\boldsymbol X$, denoted $\mathcal V(\boldsymbol X)$, is a set of convex polygons $\{V_{\boldsymbol x_1}, \ldots, V_{\boldsymbol x_m}\}$, also called _tiles_, such that 

```math
\begin{equation*}
V_{\boldsymbol x_i} = \{\boldsymbol x \in \mathbb R^2 : \|\boldsymbol x - \boldsymbol x_i\| \leq \|\boldsymbol x - \boldsymbol x_j\|,~ \boldsymbol x_j \in \boldsymbol X\}.
\end{equation*}
```

In particular, any point in $V_{\boldsymbol x_i}$ is closer to $\boldsymbol x_i$ than it is to any other point in $\boldsymbol X$. We will also denote $V_{\boldsymbol x_i}$ by $V_i$. [DelaunayTriangulation.jl](https://github.com/DanielVandH/DelaunayTriangulation.jl) is used to build $\mathcal V(\boldsymbol X)$. An example of a Voronoi tessellation is shown below.

```@raw html
<figure>
    <img src='../figs/example_tessellation.png', alt='Voronoi Tessellation'><br>
</figure>
```

# Natural Neighbours 

See that the tiles of the tessellation in the figure above intersect along a line, called the _Voronoi facet_, that we denote by $\mathcal F_{ij} = \mathcal V_i \cap \mathcal V_j$. Whenever $\mathcal F_{ij} \neq \emptyset$, we say that $\boldsymbol x_i$ and $\boldsymbol x_j$ are _natural neighbours_ in $\boldsymbol X$. We denote the set of natural neighbours to a point $\boldsymbol x \in \boldsymbol X$ by $N(\boldsymbol x) \subseteq \boldsymbol X$, and we denote the corresponding indices by $N_i = \{j : \boldsymbol x_j \in N(\boldsymbol x_j)\}$.

# Natural Neighbour Coordinates 

We represent points locally using what are known as natural neighbour coordinates, which are based on the nearby Voronoi tiles. In particular, we make the following definition: Any set of convex coordinates $\boldsymbol \lambda$ (convex means that $\lambda_i \geq 0$ for each $i$) of $\boldsymbol x$ with respect to the natural neighbours $N(\boldsymbol x)$ of $\boldsymbol x$ that satisfies: 

1.  $\lambda_i > 0 \iff \boldsymbol x_i \in N(\boldsymbol x)$,
2.  $\boldsymbol \lambda$ is continuous with respect to $\boldsymbol x$, 

is called a set of _natural neighbour coordinates_ of $\boldsymbol x$ in $\boldsymbol X$, or just the natural neighbour coordinates of $\boldsymbol x$, or the _local coordinates_ of $\boldsymbol x$.

# Natural Neighbour Interpolation 

Now that we have some definitions, we can actually define the steps involved in natural neighbour interpolation. We are supposing that we have some data $(\boldsymbol x_i, z_i)$ for $i=1,\ldots,m$, and we want to interpolate this data at some point $\boldsymbol x_0 \in \mathcal C(\boldsymbol X)$, where $\boldsymbol X$ is the point set $(\boldsymbol x_1,\ldots,\boldsymbol x_m)$ and $\mathcal C(\boldsymbol X)$ is the convex hull of $\boldsymbol x$. We let $Z$ denote the function values $(z_1,\ldots,z_m)$.

The steps are relatively straight forward.

1. First, compute local coordinates $\boldsymbol \lambda(\boldsymbol x_0)$ with respect to the natural neighbours $N(\boldsymbol x_0)$.
2. Combine the values $z_i$ associated with $\boldsymbol x_i \in N(\boldsymbol x)$ using some blending function $\varphi(\boldsymbol \lambda, Z)$.

To consider the second step, note that a major motivation for working with local coordinates is the following fact: The coordinates $\boldsymbol \lambda$ that we compute can be used to represent our point $\boldsymbol x_0$ as $\boldsymbol x_0 = \sum_{i \in N_0} \lambda_i(\boldsymbol x_0)\boldsymbol x_i$, a property known as the local coordinates property ([Sibson, 1980](https://doi.org/10.1017/S0305004100056589)), or the natural neighbour coordinates property if $\boldsymbol \lambda$ is convex (as we assume them to be). 

In particular, the coordinates $\boldsymbol \lambda$ determine how much each point $\boldsymbol x_i \in N(\boldsymbol x_0)$ contributes to the representation of our query point $\boldsymbol x_0$, hence the term "natural". So, a natural interpolant is to simply take this linear combination and replace $\boldsymbol x_i$ by $z_i$, giving the scattered data interpolant

```math
f(\boldsymbol x_0) = \sum_{i \in N_0} \lambda_i z_i.
```

Note that the natural neighbour coordinates property only works for points in the convex hull of $\boldsymbol X$ (otherwise $\boldsymbol \lambda$ is not convex), hence the restriction $\boldsymbol x_0 \in \mathcal C(\boldsymbol X)$.

# Some Local Coordinates

Let us now define all the coordinates we provide in this package.

## Nearest Neighbours 

To represent a point $\boldsymbol x_0$, we can use what are known as _nearest neighbour coordinates_, which simply assigns a weight of $1$ to the generator of the tile that $\boldsymbol x_0$ is in:

```math
\lambda_i^{\text{NEAR}} = \begin{cases} 1 & \boldsymbol x_0 \in \mathcal V_i, \\ 0 & \text{otherwise}. \end{cases}
```

The resulting scatterd data interpolant $f^{\text{NEAR}}$ is then just 

```math
f^{\text{NEAR}}(\boldsymbol x) = z_i, 
```

where $\boldsymbol x \in \mathcal V_i$. An example of what this interpolant looks like is given below.

```@raw html
<figure>
    <img src='../figs/fnear_example.png', alt='Nearest Neighbour Interpolation'><br>
</figure>
```

## Laplace Coordinates 

Here we introduce _Laplace coordinates_, also known as _non-Sibsonian coordinates_. To define these coordinates, let us take some tessellation $\mathcal V(\boldsymbol X)$ and see what happens when we add a point into it.

```@raw html
<figure>
    <img src='../figs/new_tile.png', alt='Tessellation with a new point'><br>
</figure>
```

In the figure above, the tiles with the black boundaries and no fill are the tiles of the original tessellation, and we show the tile that would be created by some query point $\boldsymbol x_0$ (the magenta point) with a blue tile. We see that the insertion of $\boldsymbol x_0$ into the tessellation has intersected some of the other tiles, in particular it has modified only the tiles corresponding to its natural neighbours in $N(\boldsymbol x_0).

For a given generator $\boldsymbol x_i \in N(\boldsymbol x_0)$, we see that there is a corresponding blue line passing through its tile. Denote this blue line by $\mathcal F_{0i}$, and let $r_i = \|\boldsymbol x_0 - \boldsymbol x_i\|$. With this definition, we define 

```math
\lambda_i^{\text{LAP}} = \frac{1}{\sum_{j \in N_0} \hat\lambda_j^{\text{LAP}}}, \quad \hat\lambda_i^{\text{LAP}} = \frac{\ell(\mathcal F_{0i})}{r_i},
```

where $\ell(\mathcal F_{0i})$ is the length of the facet $\mathcal F_{0i}$. In particular, $\hat\lambda_i^{\text{LAP}}$ is the ratio of the blue line inside the tile and the distance between the generator $\boldsymbol x_i$ and the query point $\boldsymbol x_0$. These coordinates are continuous in $\mathcal C(\boldsymbol X)$ with derivative discontinuities at the data sites $\boldsymbol X$. The resulting interpolant $f^{\text{LAP}}$ inherits these properties, where 

```math
f^{\text{LAP}}(\boldsymbol x_0) = \sum_{i \in N_0} \lambda_i^{\text{LAP}}z_i.
```

An example of what this interpolant looks like is given below.

```@raw html
<figure>
    <img src='../figs/flap_example.png', alt='Laplace Interpolation'><br>
</figure>
```

## Sibson Coordinates 

Now we consider Sibson's coordinates. These coordinates are similar to Laplace's coordinates, except we consider the areas rather than lengths for the facets. In particular, let us return to the figure above, reprinted below for convenience:

```@raw html
<figure>
    <img src='../figs/new_tile.png', alt='Tessellation with a new point'><br>
</figure>
```

The idea is to consider how much area this new blue tile "steals" from the tiles of its natural neighbours. Based on the following identity ([Sibson, 1980](https://doi.org/10.1017/S0305004100056589)),

```math
\text{Area}(\mathcal V_{\boldsymbol x_0}) \boldsymbol x = \sum_{i \in N_0} \text{Area}(\mathcal V_{\boldsymbol x} \cap \mathcal V_{\boldsymbol x_i}^{\boldsymbol x_0}),
```

where $\mathcal V_{\boldsymbol x_0}$ is the new tile shown in blue, and $\mathcal V_{\boldsymbol x_i}^{\boldsymbol x_0}$ is the tile associated with $\boldsymbol x_i$ in the original tessellation, i.e. prior to the insertion of $\boldsymbol x_0$, we define _Sibson's coordinates_ as 

```math
\lambda_i^{\text{SIB}} = \frac{1}{\sum_{j \in N_0} \hat\lambda_j^{\text{SIB}}}, \quad \hat\lambda_i^{\text{SIB}} = \text{Area}(\mathcal V_{\boldsymbol x_0} \cap \mathcal V_{\boldsymbol x_i}^{\boldsymbol x_0}).
```

A clearer way to write this is as 

```math
\lambda_i^{\text{SIB}} = \frac{A(\boldsymbol x_i)}{A(\boldsymbol x)},
```

where $A(\boldsymbol x_i)$ is the area of the intersection between the original tile for $\boldsymbol x_i$ and the new tile at $\boldsymbol x_0$, and $A(\boldsymbol x)$ is the total area of the new blue tile. These coordinates are $C^1$ continuous in $\mathcal C(\boldsymbol X) \setminus \boldsymbol X$, with derivative discontinuities at the data sites, and so too is the interpolant 

```math
f^{\text{SIB}}(\boldsymbol x_0) = \sum_{i \in N_0} \lambda_i^{\text{SIB}}z_i.
```

We may also use $f^{\text{SIB}0}$ and $\lambda_i^{\text{SIB}0}$ rather than $f^{\text{SIB}}$ and $\lambda_i^{\text{SIB}}$, respectively. 

An example of what this interpolant looks like is given below.

```@raw html
<figure>
    <img src='../figs/fsib0_example.png', alt='Sibson Interpolation'><br>
</figure>
```

Our implementation of these coordinates follows [this article](https://gwlucastrig.github.io/TinfourDocs/NaturalNeighborTinfourAlgorithm/index.html) with some simple modifications.

## Sibson-1 Coordinates 

Here we describe an extension of Sibson's coordinates, which we may also call Sibson-0 coordinates, which is $C^1$ at the data sites (but still $C^1$ in $\mathcal C(\boldsymbol X) \setminus \boldsymbol X$). A limitation of it is that it requires an estimate of the gradient $\boldsymbol \nabla_i$ at the data sites $\boldsymbol x_i$, which may be estimated using the derivative generation techniques describd in the sidebar. 

To define the interpolant for these new coordinates, denoted $f^{\text{SIB}1}$, we first define:

```math
\begin{align*}
r_i &= \|\boldsymbol x_0-\boldsymbol x_i\|, \\
\gamma_i &= \frac{\lambda_i^{\text{SIB}0}}{r_i}, \\
\zeta_i &= z_i + (\boldsymbol x_0 - \boldsymbol x_i)^T\boldsymbol\nabla_i, \\
\zeta &= \frac{\sum_{i\in N_0} \gamma_i\zeta_i}{\sum_{i\in N_0} \gamma_i}, \\
\alpha &= \frac{\sum_{i \in N_0} \lambda_i^{\text{SIB}0}r_i}{\sum_{i \in N_0} \gamma_i}, \\
\beta &= \sum_{i \in N_0} \lambda_i^{\text{SIB}0}r_i^2.
\end{align*}
```

Our interpolant is then defined by 

```math
f^{\text{SIB}1}(\boldsymbol x_0) = \frac{\alpha f^{\text{SIB}0} + \beta\zeta}{\alpha + \beta}.
```

This interpolant exactly reproduces spherical quadratics $\boldsymbol x \mapsto \mu(\boldsymbol x - \boldsymbol a)^T(\boldsymbol x - \boldsymbol a)$.

An example of what this interpolant looks like is given below.

```@raw html
<figure>
    <img src='../figs/fsib1_example.png', alt='Sibson-1 Interpolation'><br>
</figure>
```

Notice that the peak of the function is much smoother than it was in the other interpolant examples.

## Triangle Coordinates

Now we define triangle coordinates. These are not actually natural coordinates (they are not continuous in $\boldsymbol x$), but they just give a nice comparison with other methods. The idea is to interpolate based on the barycentric coordinates of the triangle that the query point is in, giving rise to a piecewise linear interpolant over $\mathcal C(\boldsymbol X)$.

Let us take our query point $\boldsymbol x=(x,y)$ and suppose it is in some triangle $V$ with coordinates $\boldsymbol x_1 = (x_1,y_1)$, $\boldsymbol x_2 = (x_2,y_2)$, and $\boldsymbol x_3=(x_3,y_3)$. We can then define:

```math
\begin{align*}
\Delta &= (y_2-y_3)(x_1-x_3)+(x_3-x_2)(y_1-y_3), \\
\lambda_1 &= \frac{(y_2-y_3)(x-x_3)+(x_3-x_2)(y-y_3)}{\Delta}, \\
\lambda_2 &= \frac{(y_3-y_1)(x-x_3) + (x_1-x_3)(y-y_3)}{\Delta}, \\
\lambda_3 &= 1-\lambda_1-\lambda_2.
\end{align*}
```

With these definitions, our interpolant is 

```math
f^{\text{TRI}}(\boldsymbol x) = \lambda_1z_1 + \lambda_2z_2 + \lambda_3z_3.
```

(Of course, the subscripts would have to be modified to match the actual indices of the points rather than assuming them to be $\boldsymbol x_1$, $\boldsymbol x_2$, and $\boldsymbol x_3$.) This is the same interpolant used in [FiniteVolumeMethod.jl](https://github.com/DanielVandH/FiniteVolumeMethod.jl).

An example of what this interpolant looks like is given below.

```@raw html
<figure>
    <img src='../figs/tri_example.png', alt='Triangle Interpolation'><br>
</figure>
```

# Extrapolation

An important consideration is extrapolation. Currently, all the methods above assume that the query point $\boldsymbol x_0$ is in $\mathcal C(\boldsymbol X)$, and the interpolation near the boundary of $\mathcal C(\boldsymbol X)$ often has some weird effects. There are many approaches available for extrapolation, such as with [ghost points](https://doi.org/10.1016/j.cad.2008.08.007), although these are not implemented in this package (yet!).

The approach we take for any point outside of $\mathcal C(\boldsymbol X)$, or on $\partial\mathcal C(\boldsymbol X)$, is to find the ghost triangle that $\boldsymbol x_0$ is in (ghost triangles are defined [here](https://danielvandh.github.io/DelaunayTriangulation.jl/dev/boundary_handling/#Ghost-Triangles) in the DelaunayTriangulation.jl documentation), which will have some boundary edge $\boldsymbol e_{ij}$. (Similarly, if $\boldsymbol x_0 \in \partial \mathcal C(\boldsymbol X)$, $\boldsymbol e_{ij}$ is the boundary edge that it is on.) We then simply use two-point interpolation, letting 

```math
f(\boldsymbol x_0) \approx \lambda_iz_i + \lambda_jz_j,
```

where $\lambda_i = 1-t$, $\lambda_j = t$, $\ell = \|x_i - \boldsymbol x_j\|$, and $t = [(x_0 - x_i)(x_j - x_i) + (y_0 - y_i)(y_j - y_i)]/\ell^2$. Note also that in this definition of $t$ we have projected $\boldsymbol x_0$ onto the line through $\boldsymbol x_i$ and $\boldsymbol x_j$ -- this projection is not necessarily on $\boldsymbol e_{ij}$, though, so $t$ will not always be in $[0, 1]$, meaning the coordinates are not guaranteed to be (and probably won't be) convex.

This extrapolation will not always be perfect, but it is good enough until we implement more sophisticated methods. If you want to disable this approach, just use the `project = false` keyword argument when evaluating your interpolant.

Similarly, if you have points defining a boundary of some domain that isn't necessarily convex, the function `identify_exterior_points` may be useful to you, provided you have represented your boundary as defined [here in DelaunayTriangulation.jl](https://danielvandh.github.io/DelaunayTriangulation.jl/dev/boundary_handling/#Boundary-Specification). See the Switzerland example in the sidebar for more information.