function [N,D] = dss2rmf(a,b,c,d,e,tol)
%DSS2RMF Converts a descriptor state space system to a right matrix fraction
%
% The command
%    [N,D] = DSS2RMF(a,b,c,d,e[,tol])
% with a, b, c, d and e all constant matrices characterizing a descriptor
% system
%    ex' = ax + bu,
%      y = cx + du,
% calculates a right coprime polynomial fraction representation
%     N D^-1
% such that
%     N D^-1 = c (xe-a)^-1 b + d,
% Here x is 's' (continuous-time) or 'z' (discrete-time) and
% the denominator matrix D is column reduced if N D^-1 is proper.

% Author: R.C.W. Strijbos, November 12, 1998.
% Copyright 1998 by Polyx, Ltd.
% Modified by J.Jezek, 13-Aug-2001,  arg checking

global PGLOBAL;
eval('PGLOBAL.VARIABLE;', 'painit;');

var = PGLOBAL.VARIABLE;

switch nargin
case 5
   tol = [];
case 6
   if ~isa(tol,'double') | length(tol) > 1  | ~isreal(tol) | ...
         tol<0 | tol>1,
      error('Invalid tolerance.')
   end
otherwise
   error('Not enough input arguments.');
end

if ~isnumeric(a)|~isnumeric(b)|~isnumeric(c)|~isnumeric(d)|~isnumeric(e) 
   error('Invalid 1st - 5th argument; should be constant matrices.');
end

[rA,cA]=size(a);[rB,cB]=size(b);[rC,cC]=size(c);
[rD,cD]=size(d);[rE,cE]=size(e);
if rA~=cA | rA~=rB | rA~=rE | rA~=cC | cB~=cD | rC~=rD | rE ~= cE
   error('Matrices of inconsistent dimensions.');
end
 
if isempty(a)
   N = pol(d);
   D = pol(eye(cB));
else
   Nr = c;
   Dr = shift(pol(e),1)-a;
   if rank(Dr) ~= rA
      error('Matrix se-a is singular.');
   end
   [Nl,Dl] = rmf2lmf(Nr,Dr,tol);
   Nl = Nl*b + Dl*d;
   [N,D] = lmf2rmf(Nl,Dl,tol);
end
if strcmp(var,'z^-1') | strcmp(var,'d')
   [N,D] = reverse(N,D,'r',tol);
end

%end .. dss2rmf
