function [B,N,fn] = build_bnd(Param,Grid,I) 
% author: Jialong Ren
% date: 9.28 2019
% Description: 
% This function computes the operators and r.h.s vectors for both Dirichlet 
% and Neumann boundary conditions. 
% 
% Input: 
% Grid = structure containing all pertinent information about the grid. 
% Param = structure containing all information about the physical problem 
% in particular this function needs the fields 
% Param.dof_dir = Nc by 1 column vector containing 
% the dof¡¯s of the Dirichlet boundary. 
% Param.dof_neu = column vector containing 
% the dof¡¯s of the Neumann boundary. 
% Param.qb = column vector of prescribed fluxes on Neuman bnd. 
% I = identity matrix in the full space 
% 
% Output: 
% B = Nc by N matrix of the Dirichlet constraints 
% N = (N-Nc) by (N-Nc) matrix of the nullspace of B 
% fn = N by 1 r.h.s. vector of Neuman contributions 


B=I(Param.dof_dir,:);
N=I; N(:,[Param.dof_dir])=[];

fn =zeros(Grid.N,1);
if isfield(Param,'dof_neu')
fn(Param.dof_neu)=Param.qb*Grid.A(Param.dof_neu)./ Grid.V(Param.dof_neu); %turn flux to source form
end
if isfield(Param,'dof_rad')
fn(Param.dof_rad)=Param.qr*Grid.A(Param.dof_f_neu)./ Grid.V(Param.dof_neu); %turn flux to source form
end