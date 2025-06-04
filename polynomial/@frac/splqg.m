function [C,regpoles,obspoles] = splqg(P,Q,R,rho,mu)
%SPLQG  Polynomial solution of the SISO LQG problem
%
% The function call
%    [C,regpoles,obspoles] = SPLQG(P,Q,R,rho,mu)
% results in the solution of the continuous-time SISO LQG problem 
% defined by:
% -Response of the measured output to the control input:
% 	   y = P(s)u
% -Response of the controlled output to the control input:
% 	   z = Q(s)u
% -Response of the measured output to the disturbance input:
% 	   y = R(s)v
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
% The compensator C(s) minimizes the steady-state 
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

% Checks

if nargin~=5,
   error('Wrong number of input arguments. (5 expected)');
end;
eval('P = frac(P); Q = frac(Q); R = frac(R);', ...
   'error(peel(lasterr));');
if length(P)~=1 | length(Q)~=1 | length(R)~=1 ,
   error('Scalar polynomials only.');
end;

d=P.den; hcoef=d.c(d.d+1);
n=P.num;
d=d/hcoef; n=n/hcoef;

x=Q.den; hcoef=x.c(x.d+1);
p=Q.num;
p=p/hcoef; x=x/hcoef;
if x~=d,
   error('Transfer functions should have the same poles.');
end;

x=R.den; hcoef=x.c(x.d+1);
q=R.num;
q=q/hcoef; x=x/hcoef;
if x~=d,
   error('Transfer functions should have the same poles.');
end;

[x,y,regpoles,obspoles]=splqg(d,n,p,q,rho,mu);

C=x/y;

%end .. splqg
