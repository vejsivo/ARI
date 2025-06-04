function r = isunimod(F,arg2)
%ISUNIMOD  Test if fraction is unimodular
%
% For square fraction F, the command
%  R = ISUNIMOD(F[,TOL])
% returns 1 if F is convertible to two-sided polynomial
% and the determinant of F is nonzero number and 0 otherwise.

% An optional tolerance TOL may be included. Its default value is the 
% global zeroing tolerance.
%
% See also POL/ISUNIMOD, TSP/ISUNIMOD.

%       Author:  J. Jezek  30-Sep-2002
%       Copyright (c) 2002 by Polyx, Ltd.

Fcl = class(F);
if strcmp(Fcl,'frac') | ~isa(F,'frac'),
   error('Invalid 1st argument.');
end;

[m,n] = size(F);
if m~=n,
   error('Matrix is not square.');
end;

F.v = 'z'; r = logical(1);
eval('Ft = tsp(F);', 'r = logical(0);');
if ~r, return;
end;
if nargin==1,
   eval('r = isunimod(Ft);', 'error(peel(lasterr));');
else
   eval('r = isunimod(Ft,arg2);', 'error(peel(lasterr));');
end;

%end .. @frac/isunimod

