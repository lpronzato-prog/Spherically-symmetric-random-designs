# Space-Filling Designs: MATLAB Implementation

## Overview
This repository provides MATLAB implementations for the construction of random designs with small $L_s$-Mean Quantisation Error of a
spherically symmetric distribution $\mu$ in a space of large dimension $d$. Three cases are considered: $\mu$ is uniform in the unit ball $B_d(0,1)$, $\mu$ is uniform on the Euclidean sphere $S_{d-1}(0,1)$, or $\mu$ is spherically normal $N(0,I_d/d)$. The designs consists of $n$ points independently distributed, uniformly in a ball $B_d(0,R)$ or on a sphere $S_{d-1}(0,R)$, or with a normal distribution $N(0,\sigma I_d/d)$. 

The script Example_random_designs_sphere_and_ball.m provides an example of constructions. All `.m` files required to determine the **optimal values** of $R$ and $\sigma$ as functions of $d$, $n$, and $s$ are provided.
 
---

## Disclaimer
- These MATLAB functions are provided as-is, in the hope that they will be useful, but **without any warranty or guarantee of fitness for a particular purpose**.
- The implementations of the algorithms from the paper arXiv:2605.12568 have not been optimised for efficiency or numerical accuracy. For example, using a Cholesky decomposition for the kernel matrices involved could improve computational efficiency and numerical stability.
- All functions have been written without the assistance of AI tools, and there is likely room for improvement.

---
## Acknowledgments
If you use these functions in your research, please cite the following paper:
"L. Pronzato, A. Zhigljavsky: Non-asymptotic quantisation of spherically symmetric distributions, https://arxiv.org/abs/2605.12568"
