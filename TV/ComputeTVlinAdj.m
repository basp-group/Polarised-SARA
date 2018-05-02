function Dtz = ComputeTVlinAdj(z,nx,ny)
 
Z = [reshape(z(1:nx*ny),ny,nx);reshape(z(nx*ny+1:end),ny,nx)];
Zv = Z(1:ny,:);        
Zh = Z(ny+1:2*ny,:);   
U = Zv;                                        
U(:,1:(nx-1)) = Zv(:,1:(nx-1))-Zv(:,2:nx);    
V = Zh;                                      
V(1:(ny-1),:) = Zh(1:(ny-1),:)-Zh(2:ny,:);    
Dtz = U + V;                                  
Dtz = Dtz(:);
 
end
 