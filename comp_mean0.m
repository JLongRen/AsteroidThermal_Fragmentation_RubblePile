function [Kd] = comp_mean0(K,p,Grid) 

        K1=K(1:Grid.Nx-1);
        K2=K(2:Grid.Nx);
        mean = zeros(Grid.Nx+1,1); 
        mean(2:Grid.Nx) = (0.5*(K1.^p+K2.^p)).^(1/p);
        mean(1)=K(1); 
        mean(Grid.Nx+1)=K(Grid.Nx);
        Kd = sparse(diag (mean));

end