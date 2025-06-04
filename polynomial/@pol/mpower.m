function [An,C] = mpower(A,n,tol)
%MPOWER (^)  Matrix power of polynomial
%              A^N
%
% AN = A^N or AN = MPOWER(A,N) is A to the N-th power where N is 
% scalar integer and A is square. If N is nonnegative, the result is
% a polynomial. If N is negative and A nonsingular, the result is
% a scalar-denominator fraction (barring the simple cases  z^N  or
% zi^N resulting in polynomials). The class of the result can be
% obtained in an optional output argument.
%
% The command works with zeroing activated through the global variable
% PGLOBAL.ZEROING. The command AN = MPOWER(A,N,TOL) works with zeroing
% specified by the input argument TOL (nonnegative scalar).
%
% See also POL/POWER, POL/MTIMES, POL/INV.

%       Author(s):  M. Hromcik, M. Sebek 16-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 02-Jul-1998 10:28:34   $
%       $Revision: 3.0 $  $Date: 05-Apr-2000   J.Jezek  $
%                         $Date: 23-May-2000   J.Jezek  $
%
% Effect on other properties:
% C.u: UserData are deleted.
% C.version: set to 3.0.

global PGLOBAL;

na = nargin;
if na<2,
   error('Not enough input arguments.');
end;

if ~isa(n,'double') | any(size(n)~=1) | ~isreal(n) | floor(n)~=n,
   error('Invalid power; must be integer scalar.');
end;

if na<3 | isempty(tol),
   tol = PGLOBAL.ZEROING;
elseif ~isa(tol,'double') | length(tol)~=1 | ...
      ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;


absn = abs(n);
nstore = absn;
A = pol(A);
me = [];

if n<0 & all([all(A.s == 1), isequal(A.c, cat(3,0,1))]),
   if strcmp(A.v,'z'),
      An = shift(pol([0 1],1,'z^-1'), -n-1);
      An.h = A.h;
      C = class(An); return;
   elseif strcmp(A.v,'z^-1'),
      An = shift(pol([0 1],1,'z'), -n-1);
      An.h = A.h;
      C = class(An); return;
   end;
end;  

As = A.s;
if As(1)~=As(2),
 error('Matrix is not square.');
end; 

if isempty(A),
   An = A;
elseif absn==0,
   An = pol(eye(As(1)));
elseif n<0,
   eval('A = inv(A,tol);','error(peel(lasterr));');
   An = mpower(A,absn,tol);
 
else 
 bin=[];
 while absn~=0,
  bin = [bin, rem(absn,2)];		
  absn = floor(absn/2);
 end; 	%while
 Z=A;
 q=1;
 while bin(q)==0,
  Z=mtimes(Z,Z,0);		% no zeroing yet
  q=q+1;
 end; 	%while
 F=Z;
 for k=q+1 : size(bin,2),
  Z=mtimes(Z,Z,0);		% no zeroing yet
  if bin(k)~=0,
   F=mtimes(F,Z,0);
  end;   %if
 end; 	%for
 An = F;
 
 % zeroing
 me = min(abs(nonzeros(A.c)));
 if ~isempty(me),
    me = me^nstore;
    An = pzer(An,tol*me);
 end;
 
end;	%if else

C = class(An);

%end .. @pol/mpower



