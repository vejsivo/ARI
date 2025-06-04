function r = ismonomod(F,arg2)
%ISMONOMOD  Test if fraction is monomodular
%
% For square fraction F, the command
%  R = ISMONOMOD(F[,TOL])
% returns 1 if F is convertible to two-sided polynomial
% and the determinant of F is monomial  K*z^k
% with K nonzero, k any integer; 0 otherwise.

% An optional tolerance TOL may be included. Its default value is the 
% global zeroing tolerance.
%
% See also POL/ISMONOMOD, TSP/ISMONOMOD.

%       Author:  J. Jezek  30-Sep-2002
%       Copyright (c) 2002 by Polyx, Ltd.

Fcl = class(F);
if strcmp(Fcl,'frac') | ~isa(F,'frac'),
   error('Invalid 1st argument.');
end;

F.v = 'z'; r = logical(1);
eval('Ft = tsp(F);', 'r = logical(0);');
if ~r, return;
end;
if nargin==1,
   eval('r = ismonomod(Ft);', 'error(peel(lasterr));');
else
   eval('r = ismonomod(Ft,arg2);', 'error(peel(lasterr));');
end;

%end .. @frac/ismonomod

