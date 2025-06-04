function An = mpower(A,n,tol)
%MPOWER(^)   Matrix power of fraction
%              AN = A^N   or   AN = MPOWER(A,N)
%
% For square fraction A and for scalar integer N,
% the command returns fraction AN, the N-th power of A.
% When N is negative, A must be nonsingular.
%
% An optional input argument TOL may specify a zeroing tolerance
% to be used instead the standard one.
%
% See also RDF/POWER, LDF/POWER, MDF/POWER, SDF/POWEr,
%          RDF/MTIMES, LDF/MTIMES, MDF/MTIMES, SDF/MTIMES,
%          RDF/INV, LDF/INV, MDF/INV, SDF/INV.

%       Author:  J. Jezek  27-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 25-Jul-2002 $
%                     $ Date 14-Oct-2002 $

global PGLOBAL;

na = nargin;
if na<2,
   error('Not enough input arguments.');
end;
if ~isa(n,'double') | any(size(n)~=1) | ~isreal(n) | floor(n)~=n,
   error('Invalid power; must be scalar integer.');
end;
if na==2 | isempty(tol),
   tol = PGLOBAL.ZEROING;
elseif ~isa(tol,'double') | length(tol)~=1 | ...
      ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;

Acl = class(A);
if strcmp(Acl,'frac'),
   error('Function ''^'' not defined for variables of class ''frac''.');
end;

As = A.s;
if As(1)~=As(2),
   error('Matrix is not square.');
end; 

if isempty(A),
   An = A; An.c = 1; An.r = 1; return;
elseif n==0,
   An = A; An.num = pol(eye(As(1))); An.den = An.num;
   An.c = 1; An.r = 1;  An.p = 1; return;
end;

if n<0,
   eval('A = inv(A,tol);','error(peel(lasterr));');
   n = -n;
end;

bin = [];
while n~=0,
   bin = [bin, rem(n,2)];		
   n = floor(n/2);
end;
Z = A;
q = 1;
while bin(q)==0,
   Z = mtimes(Z,Z,tol);
   q = q+1;
end;
F = Z;
for k = q+1 : size(bin,2),
   Z = mtimes(Z,Z,tol);
   if bin(k)~=0,
      F = mtimes(F,Z,tol);
   end;
end;
An = F;

%end .. @frac/mpower
