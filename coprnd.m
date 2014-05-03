function [ u ] = coprnd(rho, m, n )
    R(1:n,1:n) = rho; 
    R(1:n+1:n*n) = 1;
    u = copularnd('Gaussian',R,m); 
end



        
