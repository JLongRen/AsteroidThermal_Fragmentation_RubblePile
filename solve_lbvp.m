function [u] = solve_lbvp(L,f,B,g,N) 
% author: Jialong Ren
% date: 9.28 2019
% $$\mathcal{L}(u)=f \quad x\in \Omega $$ 
% % with boundary conditions 
% % $$\mathcal{B}(u)=g \quad x\in\partial\Omega$$. 
% Input: 
% L = matrix representing the discretized linear operator of size N by N, where N is the number of degrees of fredom 
% f = column vector representing the discretized r.h.s. and contributions  due non-homogeneous Neumann BC��s of size N by 1 
% B = matrix representing the constraints arising from Dirichlet BC��s of  size Nc by N 
% g = column vector representing the non-homogeneous Dirichlet BC��s of size % Nc by 1. 
% N = matrix representing a orthonormal basis for the null-space of B and % of size N by (N-Nc). 
% Output: 
% u = column vector of the solution of size N by 1 

up=B'*(B*B'\g);%particular solu
u0=N*((N'*L*N)\N'*(f-L*up)); %associated hom. solu
u=u0+up;

