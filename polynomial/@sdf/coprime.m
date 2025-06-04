function Q = coprime(M,tol)
%COPRIME   Make scalar-den fraction coprime
%
% For scalar-denominator fraction M, the command  Q = COPRIME(M)
% returns scalar fraction Q which is equal to M but
% modified so that its numerator matrix and denominator scalar
% are coprime, i.e. have no scalar common divisor.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also SDF/REDUCE.

%       Author:  J. Jezek  26-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 26-Apr-2000 $
%                     $ Date 22-Jun-2001 $
%                     $ Date 06-Oct-2002 $
%                     $ Date 14-Oct-2002 $

global PGLOBAL;

if nargin==1 | isempty(tol),
   tol = PGLOBAL.ZEROING;
elseif ~isa(tol,'double') | length(tol)~=1 | ...
      ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;

if strcmp(M.frac.c,'cop') & M.frac.tc==tol,     % quick exit
   Q = M; return;
end;

[G,R] = pgcd(vertcat(M.frac.den,M.frac.num(:)),[],tol);
Qd = R(1);
Qn = reshape(R(2:length(R)),M.frac.s);
Q = sdf(Qn,Qd);
props(Q,'cop',tol,M.frac.p,M.frac.tp);

%end .. @sdf/coprime
