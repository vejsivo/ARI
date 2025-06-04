function [y,x,gopt] = mixeds(n,m,d,a1,b1,a2,b2,gmin,gmax,accuracy,p11,p12)
%MIXEDS  Solution of the SISO mixed sensitivity problem
%
% The command
%    [y,x,gopt] = MIXEDS(n,m,d,a1,b1,a2,b2,gmin,gmax,accuracy,tol,'show')
% solves the mixed sensitivity problem for the SISO plant P = n/d
% with weighting functions V = m/d, W1 = a1/b1 and W2 = a2/b2.
%
% The input parameters gmin and gmax are lower and upper bounds for
% the minimal H-infinity norm, respectively.
%
% The parameter accuracy specifies how closely the minimal norm is
% to be approached.
%
% The optional parameter tol = [tolcncl tolstable tolspf tollr] 
% defines four tolerances. 
% -The tolerance tolcncl is used in canceling identical pole-zero pairs 
%  in the transfer function C = y/x of the optimal compensator. Its 
%  default value is 1e-4.
% -The tolerance tolstable is used in the various stability tests. It
%  has the default value 1e-8.
% -The tolerance tolspf is used in the spectral factorization, and has 
%  default value 1e-8. 
% -The tolerance tollr is used in the left-to-right and right-to-left
% conversions. Its default value is 1e-8.
%
% If the optional input argument 'show' is present then the successive 
% values of gamma during the binary search are shown, and any pole-zero 
% pairs that are canceled in y/x are reported.
%
% The output parameter gopt is an upper bound for the H-infinity norm.
%
% See also DSSHINF

% Author: H. Kwakernaak, 1997-1998
% Copyright 1998 by PolyX Ltd.
% Modified by J. Jezek, 13-Aug-2001, arg chec

% Checks

tolcncl = 1e-4; tolstable = 1e-8; tolspf = 1e-8; tollr = 1e-8;
tol = [tolcncl tolstable tolspf tollr];
show = '';

if nargin < 10
   error('Not enough input arguments.');
elseif nargin>=11,
   if ischar(p11)
      show = p11;
   elseif isnumeric(p11)
      tol = p11;
   else
      error('Invalid 11th argument.');
   end;
end;

if nargin == 12
   if ischar(p12)
      show = p12;
   elseif isnumeric(p12)
      tol = p12;
   else
      error('Invalid 12th argument.');
   end;
end;

if isempty(show), show = 0;
elseif strcmp(show,'show'), show = 1;
else error('Invalid command option.');
end;

if ndims(tol)>2 | ~all(size(tol)==[1 4])
   error('Invalid tolerance; must be a vector of length 4.');
end
if ~isreal(tol) | any(tol<0) | any(tol)>1,
   error('Invalid tolerance.');
end;
tolcncl = tol(1); tolstable = tol(2); 
tolspf = tol(3); tollr = tol(4);

eval(['d = pol(d); n = pol(n); m = pol(m);', ...
      'a1 = pol(a1); b1 = pol(b1); a2 = pol(a2); b2 = pol(b2);'], ...
   'error(peel(lasterr));');
sizes = [size(d) size(n) size(m) size(a1) size(b1) size(a2) size(b2)];
if ~all(sizes==1),
   error('Scalar polynomials only.');
end
if deg(d) ~= deg(m)
   error('The polynomials d and m have different degrees.')
end
if deg(m) > 0 & ~isstable(m,tolstable)
   error('The polynomial m is not Hurwitz.')
end
if deg(b1) > 0 & ~isstable(b1,tolstable)
   error('The polynomial b1 is not Hurwitz.')
end
if deg(b2) > 0 & ~isstable(b2,tolstable)
   error('The polynomial b2 is not Hurwitz.')
end

% Define the various polynomial matrices

D1 = [b1 0; 0 b2; 0 0];
D2 = [a1  ; 0   ; d  ];
N1 = [ 0;  0; -m ];
N2 = [ 0; a2; -n ];

% Left-to-right conversion

Q = diag([b1 b2 m]);
[Del,Lam] = lmf2rmf([N2 D2],Q,tollr);
DelLam = [Del; Lam];

% Scale Del and Lam to improve the numerical reliability

s1 = norm(DelLam(:,1),1); s2 = norm(DelLam(:,2),1); 
Del = Del*[1/s1 0; 0 1/s2]; Lam = Lam*[1/s1 0; 0 1/s2];

% Prepare the search

if ~isnumeric(gmin) | length(gmin)~=1 | ~isreal(gmin) | ...
   ~isnumeric(gmax) | length(gmax)~=1 | ~isreal(gmax) | ...   
   gmin > gmax,
      error('Invalid lower or upper bound of gamma.');
end
if ~isnumeric(accuracy) | length(accuracy)~=1 | ...
      ~isreal(accuracy) | accuracy<0,
   error('Invalid accuracy.');
end;

if show
   fprintf('\n\tgamma\ttest result\n')
   fprintf('\t-----\t-----------\n')
end

% Check existence and stability for gamma = gmax

Jgamma = eye(3); Jgamma(3,3)= -gmax^2;
DelDel = Del'*Jgamma*Del;
in = inertia(DelDel,tolstable);
if in(2) > 0
   disp('mixeds: No solution exists at the upper bound'); return
end
[Gam,J] = spf(DelDel,'ext','nnc',tolspf);
xy = rmf2lmf([1 0]*Gam,Lam,tollr);
x = xy(1,1); y = xy(1,2);
phi = d*x+n*y;
if ~isstable(phi,tolstable)
   display('mixeds: No stabilizing solution exists at the upper bound'); 
   y = []; x = []; gopt = []; return
end

% Retain stabilizing solution

xopt = x; yopt = y; gopt = gmax;
if show, fprintf('\t%g\tstable\n',gmax)
end

% Check existence and stability for gamma = gmin

Jgamma = eye(3); Jgamma(3,3)= -gmin^2;
DelDel = Del'*Jgamma*Del;
in = inertia(DelDel,tolstable);
if in(2) == 0
   [Gam,J] = spf(DelDel,'ext','nnc',tolspf);
   xy = rmf2lmf([1 0]*Gam,Lam,tollr);
   x = xy(1,1); y = xy(1,2);
   phi = d*x+n*y;
   if isstable(phi,tolstable)
      display('dessrch: A stabilizing solution exists at the lower bound'); 
      gopt = gmin; return
   else
      if show, fprintf('\t%g\tunstable\n',gmin),end
   end
else
   if show, fprintf('\t%g\tno solution\n',gmin),end
end

% Binary search

while gmax-gmin > accuracy
   gamma = (gmin+gmax)/2;
   Jgamma = eye(3); Jgamma(3,3)= -gamma^2;
   DelDel = Del'*Jgamma*Del;
   in = inertia(DelDel,tolstable);
   if in(2) > 0
      gmin = gamma;
      if show, fprintf('\t%g\tno solution\n',gamma), end
   else
      [Gam,J] = spf(DelDel,'ext','nnc',tolspf);
      xy = rmf2lmf([1 0]*Gam,Lam,tollr);
      x = xy(1,1); y = xy(1,2);
      phi = d*x+n*y;
      if ~isstable(phi,tolstable)
         gmin = gamma;
      if show, fprintf('\t%g\tunstable\n',gamma), end
      else
         gmax = gamma;
      if show, fprintf('\t%g\tstable\n',gamma),end
         xopt = x; yopt = y; gopt = gmax;
      end
   end
end

% Cancel any common roots of xopt and yopt

rootsy = roots(yopt);
xo = xopt{:}; degxo = deg(xopt);
yo = yopt{:}; degyo = deg(yopt);
for i = 1:length(rootsy)
   z = rootsy(i);

   % Divide xo(s) by s-z

   if degxo > 0
      q = zeros(1,degxo);
      q(degxo) = xo(degxo+1);
      for j = degxo-1:-1:1
         q(j) = xo(j+1)+z*q(j+1);
      end

      % If the remainder is small then cancel
      % the factor s-z both in xo and in yo

      if abs(xo(1)+z*q(1)) < tolcncl*norm(xo,1)
         xo = q; degxo = degxo-1;
         p = zeros(1,degyo);
         p(degyo) = yo(degyo+1);
         for j = degyo-1:-1:1
            p(j) = yo(j+1)+z*p(j+1);
         end
         yo = p; degyo = degyo-1;
         if show
            disp(sprintf('\nCancel root at %g\n',z));
         end
      end  
   end
end

x = pol(real(xo),degxo);
y = pol(real(yo),degyo);

%end .. mixeds

