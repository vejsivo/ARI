function [R,W] = tsyp(p0,w,epsilon)
%TSYP   Uses Tsypkin-Polyak function to find the robustness margin for a continuous 
%       interval polynomial.
%
%  TSYP(P0,W,EPSILON) uses Tsypkin-Polyak function to find and graphically display 
%  the robustness margin R. P0 is a nominal polynomial, EPSILON is a set of scale factors,
%  and W is a vector of frequencies. Then the robustness margin can be defined as maximal R 
%  such that the polynomial family pr(s,q) = p0(s) + R*sum([-epsilon(i), epsilon(i)]*s^i), 
%  where the sum is taken from 0 to n, is robustly stable.
%
%  TSYP(P0,W) uses the coeficients of the P0 polynomial to set the scale factors.
%
%  TSYP(P0) estimates the range of frequencies based on the fastest root of the nominal polynomial.
%
%  TSYP(P0,[],EPSILON)...
%
%  R = TSYP(...) Supresses the graphical output and assigns the robustness margin to R
%  
%  [R,W] = TSYP(P0) As the second output returns the vector of frequencies, which can
%  be further used by KHPLOT, for instance.
%
%  See also: KHPLOT, KHARIT.

%  See Barmish (1994): New Tools for Robustness of Linear Systems, pp.97.

%-----------------------------------------------------------------------------------------
%  Author(s):    Z. Hurak   04-10-2000
%  Copyright (c) 2000 by Polyx, Ltd.
%  $Revision: 2.5.0 $    $Date: 11-22-2000 $
%-----------------------------------------------------------------------------------------



global PGLOBAL;
eval('PGLOBAL.FORMAT;', 'painit;');

if strcmp(PGLOBAL.VERBOSE, 'yes'); %conversion from 'yes/no' to 'on/off'
   verbose = 'on';
else
   verbose = 'off';
end

%====================================================================
%Tests for proper inputs:
%====================================================================
if nargin < 1
   error('Not enough input arguments.');

elseif nargin == 2
   epsilon = p0{:};
   if ~isa(p0,'pol') | ~isnumeric(w)
      error('Invalid arguments.')
   end

elseif nargin == 3
   if ~isnumeric(epsilon) | ~isa(p0,'pol') | ~isnumeric(w)
      error('Invalid arguments.')
   end
   if isempty(w)
      wmax = ceil(max(abs(roots(p0))));
      w = linspace(0,wmax,10000);
   end
   if length(epsilon) ~= (p0.deg+1)
      error('The size of the scale factors vector does not agree with the nominal polynomial.');
   end
   if any(epsilon<0)
      error('The scales can only be nonnegative.');
   end

elseif nargin == 1
   if ~isa(p0,'pol')
      error('Invalid arguments.')
   end
   epsilon = p0{:};
   wmax = ceil(max(abs(roots(p0))));
   w = linspace(0,wmax,10000);
end

sy = symbol(p0);
if (sy ~= 's') & (sy ~= 'p')
   error('Continuous polynomial only.');
end

if ~isstable(p0)
   error('Nominal polynomial must be stable.');
end


%=========================================================================
%Finds the minimum infinity norm of the Tsypkin-Polyak function
%=========================================================================

w = w + eps;

Gtp = gtp(w,p0,epsilon);
Gtp_inf = max(abs(real(Gtp)),abs(imag(Gtp)));

[Gtp_inf_ordered,pp] = sort(Gtp_inf); %bounds on the frequency
w1 = min(w(pp(1)),w(pp(3)));
w2 = max(w(pp(1)),w(pp(3)));

Gtp_inf = inline('max(abs(real(gtp(om,po,ep))),abs(imag(gtp(om,po,ep))))',...
   'om','po','ep');   %infinity norm of the Tsypkin-Polyak function

options = optimset('Display',verbose,'Diagnostics',verbose);
[wmin,Gtp_inf_min] = fminbnd(Gtp_inf,w1,w2,options,p0,epsilon);
                      %finds the minimum of the infinity norm on the interval <w1,w2>

Rw = 0.997*Gtp_inf_min;                         %tolerance (3e-3);
R0 = abs(p0{0})/epsilon(1);                     %zero frequency condition
RN = abs(p0{deg(p0)})/(epsilon(deg(p0)+1)+eps); %degree invariance condition
R = min(RN,min(Rw,R0));         

%=========================================================================
%Stability test
%=========================================================================

pmin = p0 - R*pol(epsilon,p0.deg);     %'safety' test if the resulting interval polynomial 
pmax = p0 + R*pol(epsilon,p0.deg);     %is stable
[stab,k1,k2,k3,k4] = kharit(pmin,pmax);
if (stab ~= 1)
   warning(sprintf('Resulting margin does not guarantee robust stability of the interval polynomial.\nRun again with extended frequency range and/or denser gridding.'));
end

%=========================================================================
%Graphical output
%=========================================================================

if nargout == 0
   plot(Gtp);
   axis equal;
   ax_lim = axis;
   title('Tsypkin-Polyak plot of the \fontname{ScriptS}l_\propto\fontname{Arial}robustness margin for an interval polynomial');
   xlabel('Real');
   ylabel('Imaginary');
   axis manual;
   hold on,
   plot([-R-j*R R-j*R R+j*R -R+j*R -R-j*R],'r');
   Rmax = max(abs(ax_lim(1:2)));
   Imax = max(abs(ax_lim(3:4)));
   plot([-Rmax,Rmax,NaN,0,0],[0,0,NaN,-Imax,Imax],':k');
   hold off
   
end

if nargout == 2
   W = w;
end

%end .. tsyp
