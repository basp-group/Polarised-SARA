function Dx = ComputeTVlin(x)
%  function Dx = ComputeTVlin(X)
 

nx = 100;
ny = 100;

X = reshape(x,ny,nx);
nx = size(X,1);
ny = size(X,2);
Du = X;
Du(:,2:nx) = X(:,2:nx)-X(:,1:(nx-1));   
Dv = X;
Dv(2:ny,:) = X(2:ny,:)-X(1:(ny-1),:);   
Dx = [Du(:);Dv(:)];                     
 
end