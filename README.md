# Quantisation of spherically symmetric distributions: MATLAB Implementation

## Overview
This repository provides MATLAB implementations for the construction of random designs with small $L_s$-Mean Quantisation Error of a
spherically symmetric distribution $\mu$ in a space of large dimension $d$. Three cases are considered: $\mu$ is uniform in the unit ball $B_d(0,1)$, $\mu$ is uniform on the Euclidean sphere $S_{d-1}(0,1)$, or $\mu$ is spherically normal $N(0,I_d/d)$. 


The script Example_random_designs_sphere_and_ball.m provides examples of constructions for the case where the designs consist of $n$ points independently distributed, uniformly in a ball $B_d(0,R)$ or on a sphere $S_{d-1}(0,R)$, or with a normal distribution $N(0,\sigma I_d/d)$. 
All `.m` files required to determine the **optimal values** of $R$ and $\sigma$ as functions of $d$, $n$, and $s$ are provided.

The script Example_construction_optimal_mixure.m provides examples of constructions of optimal distributions for the quantisation of the uniform measure in the unit ball $B_d(0,1)$ and of the spherically normal distribution $N(0,I_d/d)$. The optimal distribution is uniform on $m$ concentric spheres, where  $m$ depends on $n$ for given $s$ and $d$. The optimality of the distribution is assessed by plotting the sensitivity function as a function of the radius.
 
---

## Disclaimer
- These MATLAB functions are provided as-is, in the hope that they will be useful, but **without any warranty or guarantee of fitness for a particular purpose**.
- The implementations of the algorithms from the paper mentionned below have not been optimised for efficiency or numerical accuracy. 
- All functions have been written without the assistance of AI tools, and there is likely room for improvement.

---
## Acknowledgments
If you use these functions in your research, please cite the following paper:
"L. Pronzato, A. Zhigljavsky: Optimal random quantisers for spherically symmetric distributions, 
https://arxiv.org/pdf/2610.11772 
https://hal.science/hal-05786349