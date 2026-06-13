# Bi-directional Dynamic Programming for the k-color Shortest Path Problem

This repository contains the implementation and experimental evaluation of a **Bounded Bi-directional Dynamic Programming (BDP)** approach to solve the NP-hard $k$-color Shortest Path Problem ($k$-CSPP).

## Problem Description
The $k$-color Shortest Path Problem ($k$-CSPP) arises in various real-world applications such as reliable telecommunication network design, multimodal transportation, obstacle minimization in robotics, and wavelength-routed optical networks. 

Given a directed graph with positively weighted and colored edges, a source node $s$, a target node $t$, and a positive integer $k$, the objective is to find the minimum-cost path from $s$ to $t$ that uses **at most $k$ distinct colors**.

## Proposed Approach: Bounded Bi-directional DP
This project explores a bi-directional alternative to existing mono-directional dynamic programming methods. The BDP algorithm:
* **Explores simultaneously:** Extends labels both forward from the source $s$ and backward from the destination $t$.
* **Efficient Fathoming:** Uses an Upper Bound (UB) computed via the Color-Constrained Dijkstra (CCDA) heuristic to fathom infeasible states early.
* **Half-way Point Strategy:** Halts exploration dynamically when the path cost reaches half of the computed UB. This reduces redundant label generation while ensuring no optimal solutions are lost.
* **Dominance Rules:** Prunes the search space by eliminating paths that are strictly dominated in terms of both distance and color set inclusivity.

## Implementation Details
* **Language:** Julia 1.12.6
* **Modeling Framework:** JuMP
* **Solver (for baselines):** IBM ILOG CPLEX 22.1 (Configured with 1 thread)

## Instances
The algorithm is evaluated against existing methodologies (Mono-directional DP, Valid Inequalities Branch-and-Cut, CPLEX) on several benchmark sets:
1.  **Random and Grid Digraphs (Ferone et al., 2019)** 
2.  **Layered-based Digraphs (Castro et al., 2024)** 
3.  **Labyrinth Instances (Newly Generated)** 

### Report of the Project
[Link to Report](https://github.com/echo85/k-color-shortest-path/blob/main/report.pdf)
