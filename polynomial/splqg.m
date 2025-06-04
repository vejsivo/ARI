function [y,x,regpoles,obspoles] = splqg(d,n,p,q,rho,mu)
%SPLQG  Polynomial solution of the SISO LQG problem
%
% The function call
%    [y,x,regpoles,obspoles] = SPLQG(d,n,p,q,rho,mu)
% results in the solution of the continuous-time SISO LQG problem 
% defined by:
% -Response of the measured output to the control input:
% 	   y = P(s)u,  P(s) = n(s)/d(s)
% -Response of the controlled output to the control input:
% 	   z = Q(s)u,  Q(s) = p(s)/d(s)
% -Response of the measured output to the disturbance input:
% 	   y = R(s)v,  R(s) = q(s)/d(s)
% In state space form:
%    x' = Ax + Bu + Gv
%    z  = Dx 
%    y  = Cx + w
%
%    P(s) = C*inv(sI-A)*B  
%    Q(s) = D*inv(sI-A)*B 
%    R(s) = C*inv(sI-A)*G
% The scalar white state noise v has intensity 1 and the white
% measurement noise w has intensity mu.
%
% The compensator C(s) = y(s)/x(s) minimizes the steady-state 
% value of
%    E[z^2(t) + rho*u^2(t)]
% The output argument regpoles contains the regulator poles and 
% obspoles contains the observer poles. Together they are the 
% closed-loop poles.

% Author: H. Kwakernaak, 1997
% Copyright 1998 by PolyX Ltd.



% Initialization

global PGLOBAL;
eval('PGLOBAL.FORMAT;', 'painit;'); 
verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

% Checks

if nargin<6,
   error('Not enough input arguments.');
end;
eval('d = pol(d); n = pol(n); p = pol(p); q = pol(q);', ...
   'error(peel(lasterr));');
if length(d)~=1 | length(n)~=1 | length(p)~=1 | length(q)~=1,
   error('Scalar polynomials only.');
end;

v = {d.v, n.v, p.v, q.v};
u = ''; t = 0;
for i = 1:4,
   w = v{i};
   if ~(isempty(w) | strcmp(w,'s') | strcmp(w,'p')),
      error('Continuous time polynomials only.');
   end;
   if isempty(u),
      u = w;
   elseif ~t & ~isempty(w) & ~strcmp(u,w),
      warning('Inconsistent variables.');
      t = 1;
   end;
end;
d.v = u; n.v = u; p.v = u; q.v = u;

if deg(n) >= deg(d)   
   error('n/d is not proper.'); 
elseif deg(p) >= deg(d)
   error('p/d is not proper.'); 
elseif deg(q) >= deg(d)
   error('q/d is not proper.'); 
end

if ~isa(rho,'double') | length(rho)~=1 | ~isreal(rho) | rho<0,
   error('Invalid 5th argument.');
end;
if ~isa(mu,'double') | length(mu)~=1 | ~isreal(mu) | mu<0,
   error('Invalid 6th argument.');
end;

% Spectral factorizations and computation of 
% the closed-loop polynomial phi

regpol = d'*d+p'*p/rho;
obspol = d'*d+q'*q/mu;
regpoles = roots(regpol);
obspoles = roots(obspol);
regpoles = regpoles(find(real(regpoles)<0));
obspoles = obspoles(find(real(obspoles)<0));
phir = real(poly(regpoles)); phir = pol(fliplr(phir),length(phir)-1);
phio = real(poly(obspoles)); phio = pol(fliplr(phio),length(phio)-1);
phi = phir*phio;

% Solve for the compensator

[x,y] = axbyc(d,n,phi,'miny');

%end .. splqg
