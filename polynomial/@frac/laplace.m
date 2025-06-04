function [H,Q] = laplace(F,degree,T,tau,tol)
%LAPLACE    Inverse Laplace transform of fraction
%
% For fraction F in variable 's' or 'p', the command
%    [H,Q] = LAPLACE(F,DEGREE,T,TAU)
% where scalars T,TAU>=0 are sampling period and sampling
% phase, returns polynomial H in 'z^-1', corresponding
% to the inverse Laplace transform% of the strictly proper
% part of F, sampled with period T and phase TAU, i.e.
% in sampling points  k*T+TAU, k = 0,1,2,...DEGREE.
% In Q, the polynomial part of F (corresponding to Dirac
% impulses in time domain) is returned. Q is nonzero only
% for TAU = 0, as the Dirac pulses sit in TIME = 0 .
%
% TAU can also be, instead of scalar, a vector. In such
% a case, the results H,Q form cell vectors.
%
% The arguments DEGREE,T,TAU may be omitted or given by [],
% the defaults being  DEGREE = MAX(DEG(F.NUM),DEG(F.DEN)),
%  T = 1, TAU = 0 . A tolerance TOL may be specified
% as an optional argument.
%
% See also RDF/SAMP, LDF/SAMP, MDF/SAMP, SDF/SAMP,
% RDF/LAURENT, LDF/LAURENT, MDF/LAURENT, SDF/LAURENT.

%      Author: J.Jezek, 14-Jul-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 16-Jan-2001 $
%                    $ Date 25-Jul-2002 $
%                    $ Date 14-Oct-2002 $
%                    $ Date 28-Feb-2003 $

global PGLOBAL;

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;

if ~isa(F,'frac'),
   error('Invalid 1st argument.');
end;

if ni>=2 & ~isempty(degree),
   if ~isa(degree,'double') | length(degree)~=1 | ...
         ~isreal(degree) | degree<0 | floor(degree)~=degree,
      error('Invalid degree.');
   end;
else degree = max(deg(F.num),deg(F.den));
end;

if ni>=3 & ~isempty(T),
   if ~isa(T,'double') | length(T)~=1 | ~isreal(T) | T<0,
      error('Invalid sampling period.');
   end;
elseif ~isempty(F.h) & isfinite(F.h) & F.h>0,
   T = F.h;
else
   T = 1;
end;

if ni>=4 & ~isempty(tau),
   if ~isa(tau,'double') | ndims(tau)~=2 | all(size(tau)~=1) | ...
         ~isreal(tau) | any(tau<0),
      error('Invalid sampling phase.');
   end;
   tau = sort(tau);
else tau = 0;
end;

if ni==5 & ~isempty(tol),
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
else tol = PGLOBAL.ZEROING;
end;

Fcl = class(F);
if strcmp(Fcl,'frac'),
   error('Invalid 1st argument.');
end;

Fvar = F.v;
if ~isempty(Fvar) & ~strcmp(Fvar,'s') & ~strcmp(Fvar,'p'),
   error('Invalid variable symbol; must be continuous-time.');
end;

[sF1,sF2] = size(F);
Hc = zeros(sF1,sF2*(degree+1));
[A,B,C,D] = abcd(F,tol); D = pol(D);
Dz = pol(zeros(size(D)));
R = expm(A*T);

ltau = length(tau); stau = size(tau);
if ltau>1,
   H = cell(stau); Q = H;
end;

for k = 1:ltau,
   tauk = tau(k);
   
   X = expm(A*tauk);
   HH = pol(Hc,degree);
   HH{0} = C*X*B;
   for i = 1:degree,
      X = X*R;
      HH{i} = C*X*B;
   end;
   HH.v = 'z^-1'; HH.h = T;
   
   if tauk==0, QQ = D;
   else QQ = Dz;
   end;
   
   if ltau==1,
      H = HH; Q = QQ;
   else
      H{k} = HH; Q{k} = QQ;
   end;
end;

%end .. @frac/laplace

      