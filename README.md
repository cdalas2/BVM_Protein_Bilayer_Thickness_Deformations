# CSCI596FINAL
<img src="CloverSurfaces.png " width="250" height="750">

> A boundary value method for protein-induced bilayer deformation calculations

---

### Table of Contents

- [Description](#description)
- [Benchmarks](#benchmarks)
- [How To Use](#computational-problem)
- [References](#references)
- [Author Info](#author-info)

---

## Description

#### Physical problem

Using a classical elastic continuum model for lipid bilayers, what is the lipid bilayer deformation induced by a protein of arbitrary shape?

#### Framing of the problem in mathematics 

<img src="equations.png " width="1500" height="250">
We truncation the infinite series that represents general form of the bilayer deformation field to N terms and use the truncated form in applying the boundary conditions along the protein-bilayer interface to form a linear system of boundary equations which can be solved for the coefficients A_n,B_n.


#### Computational obstacles and solution

Since Basset functions grow exponentially with increasing order n, floating point overflow issues and matrix conditioning issues manifest. We use Arblib library for arbitrary precision floating point calculations, choosing an appropriate precision to avoid overflow and to offset numerical instability when solving the linear system of boundary equations. 

[Back To The Top](#CSCI596FINAL)

---

#### Software

- C++
- Arblib library (version 2.17 or newer)
- Carlos_Membrane_Project
- OpenMP
- OpenMPI
- Paraview

[Back To The Top](#CSCI596FINAL)

---

## Benchmarks

#### Boundary value method benchmarked against the finite element method
The finite element method converges with decreasing average mesh length.
Using paraview we can create a pipeline to extract the length data of the mesh elements used in the finite elements method. We average the mesh lengths import the length data into MATLAB and average them. 

<img src="FEM_mesh.png " width="1000" height="625">

<img src="FEM_edges.png " width="1000" height="1000">

<img src="FEM_celldata.png " width="1000" height="1000">

We found that an accuracy of about 0.02% can be expected with finite elements of average mesh lengths of about 0.1 nm.
Testing the boundary value method for several protein shapes, we found that as the boundary value method converges to an agreement with the finite element method of  about 0.02% with increasing truncation length N, which was noted earlier as the expected accuracy of the finite element method with the actual solution.

#### Thread count speed up and efficiency benchmarks

[Back To The Top](#CSCI596FINAL)

## How To Use

1. $ make clover
2. $ bash job

#### Installation

1. Download github repo. This will contain the necessary make file, bash script, and source code files
2. Install Arblib library by Frederick Johansson (version 2.17 or newer). You can build from source [1] or download as a package through anaconda3.
3. Install OpenMP.
4. Install OpenMPI (if you want to do multiple protein runs)

[Back To The Top](#CSCI596FINAL)

---

## References
[Back To The Top](#CSCI596FINAL)

---

## Author Info

- LinkedIn - [Carlos Alas LI](https://www.linkedin.com/in/carlos-alas-6a4643160/)
- ResearchGate - [Carlos Alas RG](https://www.researchgate.net/profile/Carlos_Alas3)

[Back To The Top](#CSCI596FINAL)
