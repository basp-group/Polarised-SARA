function z = prox_tv(coef,x,lambda,nx,ny)
 
u = x(1:nx*ny);
v = x(nx*ny+1:end);

zu = zeros(nx*ny,1);     
zv = zu;
coeff = coef*lambda;


 
sqrtuv = sqrt(u.^2 + v.^2);


zu(sqrtuv>coeff) = (1-coeff./sqrtuv(sqrtuv>coeff)).*u(sqrtuv>coeff);    
zv(sqrtuv>coeff) = (1-coeff./sqrtuv(sqrtuv>coeff)).*v(sqrtuv>coeff);

% zu(sqrtuv>coeff) = (1-coeff(sqrtuv>coeff)./sqrtuv(sqrtuv>coeff)).*u(sqrtuv>coeff);    
% zv(sqrtuv>coeff) = (1-coeff(sqrtuv>coeff)./sqrtuv(sqrtuv>coeff)).*v(sqrtuv>coeff);

z = [zu;zv];
end