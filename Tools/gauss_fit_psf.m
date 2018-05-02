function psf_fit = gauss_fit_psf(psf,varargin)
if isempty(varargin)
    scale = 0.8;
else
    scale =cell2mat(varargin);
end
%%
scale = 1;
N  =  size(psf);
h = N(1);
w = N(2);
[X, Y] = meshgrid(1:h,1:w);
X = X(:); Y=Y(:); Z = psf(:);
x0 = floor(w/2)+mod(w/2+1,2); % guess position (center seems a good place to start)
y0 = floor(h/2)+mod(h/2+1,2); 

a = @(theta, sigmax,sigmay) cos(theta)^2/(2*sigmax^2) + sin(theta)^2/(2*sigmay^2);
b = @(theta, sigmax, sigmay) -sin(2*theta)/(4*sigmax^2) + sin(2*theta)/(4*sigmay^2);
c =@(theta, sigmax, sigmay)  sin(theta)^2/(2*sigmax^2) + cos(theta)^2/(2*sigmay^2);
% Z = A*exp( - (a*(X-x0).^2 - 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) ;

% 2D gaussian fit object
gauss2 = fittype( @(sigmax, sigmay,theta,x,y) ...
    exp(-(a(theta,sigmax,sigmay).*(x-x0).^2-2*b(theta,sigmax,sigmay).*(x-x0).*(y-y0)+c(theta, sigmax,sigmay).*(y-y0).^2)),...
    'independent', {'x', 'y'},'dependent', 'z' );
% sigmax = 2; % guess width
% sigmay = 2; % guess width

sigmax =2.1233;
sigmay =2.1233;

% theta = pi/10;
theta = 0;

% compute fitf
% sf = fit([X,Y],double(Z),gauss2,'StartPoint',[sigmax, sigmay,theta]);
% figure(6); clf; plot(sf,[X,Y],Z);
% sigmax = scale * sf.sigmax;
% sigmay = scale * sf.sigmay;
% theta = sf.theta;
sigmax = scale * sigmax;
sigmay = scale * sigmay;
gfit = @(x,y) exp(-(a(theta,sigmax,sigmay).*(x-x0).^2-2*b(theta,sigmax,sigmay).*(x-x0).*(y-y0)+c(theta, sigmax,sigmay).*(y-y0).^2));

xx = 1:h;
yy = 1:w;
[X,Y] = meshgrid(xx,yy);
psf_fit = gfit(X,Y);


psf_fit =psf_fit./max(psf_fit(:));
%%
% end
 
xx = conv2(im{1},psf_fit,'same');

xx2 = xx;

xx = xx.*(max(im{1}(:))/max(xx(:)));

num = abs(xx - im{1}).^2;
num = sum(num(:));
den = abs(im{1}).^2;
den = sum(den(:));
nrmse = sqrt(num/den);

