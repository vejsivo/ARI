function sphspectrum(p0,r,varargin)
%SPHSPECTRUM plot the spectral set for a spherical polynomial family.
%
% SPHSPECTRUM(p0,r) plot a spectral set for a spherical polynomial 
% family p(s,q) = p0(s) + sum_{k=1}^{n}q_{k}s^{k}, where norm(q,2)<=r
% for a range of values of the parameter r.
%
% SPHSPECTRUM(p0,r,x,y) plot a spectral set for a spherical polynomial 
% family p(s,q) = p0(s) + sum_{k=1}^{n}q_{k}s^{k}, where norm(q,2)<=r
% for a range of values of the parameter r.The parameters x and y are 
% vectors of values at real and imaginary axes for which the gridding is 
% performed.
%
% SPHSPECTRUM(...,'fill') fills the contours describing the spectral set.

% Author(s):    Z. Hurak   08-23-2002
% Copyright (c) 2002 by PolyX, Ltd.
%-----------------------------------------------------------------------------

global PGLOBAL;
try,
    PGLOBAL.FORMAT;
catch,
    painit
end

tol = PGLOBAL.ZEROING;

%=============================================================================
% Tests for proper inputs: 
%=============================================================================

if nargin < 2
    error('Not enought input arguments.');
end

try,
    p0 = pol(p0);
catch,
    error(lasterr);
end

n = p0.deg;
type = 'cont';

switch length(varargin)
  case 0
      flag_boundsGiven = logical(0);
  case 1
      type = varargin{:};
      if ~strcmp(type,'fill')
          error('The third argument not valid.');
      end
      flag_boundsGiven = logical(0);
  case 2
      [x,y] = deal(varargin{:});
      if ~(isnumeric(x) & isnumeric(y))
          error('Incorrect type of the third and fourth arguments.');
      end
      flag_boundsGiven = logical(1);
  case 3
      [x,y,type] = deal(varargin{:});
      if ~(isnumeric(x) & isnumeric(y))
          error('Incorrect type of the third and fourth arguments.');
      end
      if ~strcmp(type,'fill')
          error('The fifth argument not valid.');
      end  
      flag_boundsGiven = logical(1);
  otherwise
      error('Too many input arguments.');
end

p0_re = real(roots(p0));
p0_im = imag(roots(p0));
rootRemin = min(p0_re);
rootRemax = max(p0_re);
rootImmin = min(p0_im);
rootImmax = max(p0_im);

if flag_boundsGiven
    if (x(1)>rootRemin) | (x(end)<rootRemax) | (y(1)>rootImmin) | (y(end)<rootImmax)
        warning('The assigned bounds do not cover all the roots of the nominal polynomial.');
    end
end
	
%=============================================================================
% Guessing at the bounds on the real and imaginary parts: 
%==========================================================================

if ~flag_boundsGiven
    M = 2; N = 40;
    xmin = rootRemin-M*abs(rootRemin);
    xmax = rootRemax+M*abs(rootRemax);
    ymin = rootImmin-M*abs(rootImmin);
    ymax = rootImmax+M*abs(rootImmax);
    x = linspace(xmin,xmax,N);
    y = linspace(ymin,ymax,N);
end

rA = length(y); cA = length(x);

%=============================================================================
% Computing the values of the spherical boundary function: 
%=============================================================================

[X,Y] = meshgrid(x,y); A = X+j*Y;

phi = zeros(rA,cA);

for ii = 1:rA
    for jj = 1:cA
        if abs(imag(A(ii,jj))) <= tol
            phi(ii,jj) = polyval(p0,A(ii,jj))/sqrt(polyval(sum(mono(0:n).^2),A(ii,jj)));
        else
            p0_vec = [real(polyval(p0,A(ii,jj)));imag(polyval(p0,A(ii,jj)))];
            z_i = polyval(mono(0:n).',A(ii,jj));
            re_z_i = real(z_i);
            im_z_i = imag(z_i);
            detAAT = sum(re_z_i.^2)*sum(im_z_i.^2)-sum(re_z_i.*im_z_i)^2; 
            adjAAT(1,1) = sum(im_z_i.^2);
            adjAAT(1,2) = -sum(re_z_i.*im_z_i);
            adjAAT(2,1) = -sum(re_z_i.*im_z_i);
            adjAAT(2,2) = sum(re_z_i.^2);
            phi(ii,jj) = sqrt(p0_vec.'*adjAAT*p0_vec/detAAT);
        end
    end
end

if length(r) == 1
    r = [r r];
end

if strcmp(type,'cont')
    contour(X,Y,phi,r);xlabel('Real');ylabel('Imaginary');
else
    contourf(X,Y,phi,r);xlabel('Real');ylabel('Imaginary');
end


%=============================================================================
% End of ..\sphspecsetplot.m 
%=============================================================================