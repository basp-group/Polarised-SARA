function [u,v,w,R]=util_uv_read_old(stringname, gen_data, nbr)
% from Arwa
% uvread reads realistic uv coverages
% Input:
% stringname   - u,v,w, cooordinates file
% nber            - number of uv points to keep
% Output:
% u,v,w            - u,v,w coordinates

uvw = load(stringname);
xyz = uvw.uvw;
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);

if gen_data == 0
     flag = uvw.flag_I;
%     flag = uvw.flag;
    flag = flag(:);
    x = x(~flag);
    y = y(~flag);
    z = z(~flag);
end


% get rid of 0 components
xdummy=x(abs(x)+abs(y)>0);
ydummy=y(abs(x)+abs(y)>0);
zdummy=z(abs(x)+abs(y)>0); 
xyz=[xdummy ydummy zdummy];
sz = size(zdummy, 1);

if nargin==3
   index=randi([1 sz],1,nbr);
   x=xyz(index,1);
   y=xyz(index,2);
   z=xyz(index,3);    
end

xyz=[x'; y';z'];

uvw = xyz;
bmax=max(sqrt(uvw(1,:).^2+uvw(2,:).^2));

uvw=pi*uvw./(bmax);
u=uvw(1,:)';
v=uvw(2,:)';
w=uvw(3,:)';

R=length(u);

end