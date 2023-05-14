# NaturalNeighbourInterp

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://DanielVandH.github.io/NaturalNeighbourInterp.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://DanielVandH.github.io/NaturalNeighbourInterp.jl/dev/)
[![Build Status](https://github.com/DanielVandH/NaturalNeighbourInterp.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DanielVandH/NaturalNeighbourInterp.jl/actions/workflows/CI.yml?query=branch%3Amain)

This is a package for performing [natural neighbour interpolation](https://en.wikipedia.org/wiki/Natural_neighbor_interpolation) over planar data sets. This method of (scattered data) interpolation takes in some data $X = \{(x_i,y_i)\}_{i=1}^m \subset \mathbb R^2$ with corresponding data values $Z = \{z_i\}_{i=1}^m$ and constructs a function $f \colon \mathbb R^2 \to \mathbb R$ such that $f(x_i, y_i) = z_i$, $i=1,\ldots,m$, based on the _Voronoi tessellation_ of $X$. We use [DelaunayTriangulation.jl](https://github.com/DanielVandH/DelaunayTriangulation.jl) to construct the Voronoi tessellations.

To explain this in more detail, let us make some definitions:

1. _Voronoi tessellation_: Let $X = \{(x_i, y_i)\}_{i=1}^m \subseteq \mathbb R^2$ be a set of points in the plane. We define the _Voronoi tessellation_ $\mathcal V(X) := \{\mathcal V_i\}_{i=1}^m$ to be the partition of $\mathbb R^2$ into convex polygons $\mathcal V_i$ such that all points $(x, y) \in \mathcal V_i$ are closer to the _generator_ $(x_i, y_i)$ than to any other generator $(x_j, y_j)$, $i \neq j$, meaning
$$
\mathcal V_i = \{(x, y) \in \mathbb R^2 \mid (x - x_i)^2 + (y - y_i)^2 \leq (x - x_j)^2 + (y - y_j)^2, i \neq j\},
$$

2. _Natural neighbours_: Two points $(x_i, y_i)$ and $(x_j, y_j)$ are _natural neighbours_ if the intersection $\mathcal V_i \cap \mathcal V_j$ is not empty. The set of natural neighbours to $(x_i, y_i)$ will be denoted $\mathcal N_i$.

Using natural neighbours, we can define a set of generalised barycentric coordinates $\lambda_i(x_0, y_0)$ (also called _weights_ or _local coordinates_) about the natural neighbours $\mathcal N_0$ of some query point $(x_0, y_0)$, allowing us to write:

1. $(x_0, y_0) = \sum_{i \in \mathcal N_0} \lambda_i(x_0, y_0)x_i$. (This is the _local coordinates_ property.)
2. $1 = \sum_{i \in \mathcal N_0} \lambda_i(x_0, y_0)$. (This is the _partition of unity_ property.)
3. For all $i \in \mathcal N_0$, $\lambda_i(x_0, y_0) \geq 0$, meaning $(x_0, y_0)$ is a _convex combination_ of the $x_i \in \mathcal N_i$ or the coordinates $\lambda_i$ are _convex_.

Using these three properties, our interpolation $f$ is defined by 

$$
f(x_0, y_0) = \sum_{i \in \mathcal N_0} \lambda_i(x_0, y_0)z_i.
$$

A basic definition for $\lambda_i(x_0, y_0)$ is through _Sibson's coordinates_, which defines 

$$
\hat\lambda_i = \operatorname{Area}(\mathcal V_0' \cap \mathcal V_i), \quad \lambda_i = \left(\sum_{j \in \mathcal N_0} \hat\lambda_j\right)^{-1}\hat\lambda_i.
$$

Here, $\mathcal V_0'$ is defined as the _virtual tile_ of $(x_0, y_0)$, and is the polygon that would appear in the Voronoi tessellation $\mathcal V(X \cup \{(x_0, y_0)\})$ with $(x_0, y_0)$ inserted. Thus, $\hat\lambda_i$ is the amount of area _stolen_ from $\mathcal V_i$ from the insertion of $(x_0, y_0)$, and $\lambda_i$ is the normalised version.