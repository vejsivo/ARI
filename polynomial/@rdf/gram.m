function G = gram(F,tol)
%GRAM   Gramian of a stable right-den fraction
%
% The command  GRAM(F)  computes the Gramian of the stable
% right-den fraction  F = N/D  with N and D right coprime.
%
% A tolerance TOL may be specified as an additional input
% argument. Its default value is the global zeroing tolerance.
%
% See also POL/GRAM.

%    Author: J.Jezek, 17-Jul-2001
%    Copyright(c) 2001 by Polyx, Ltd.

global PGLOBAL;

if nargin==1 | isempty(tol),
   tol = PGLOBAL.ZEROING;
elseif ~isa(tol,'double'),
   error('Invalid tolerance.');
end;

eval('G = gram(F.frac.n,F.frac.d,''r'',tol);', ...
   'error(peel(lasterr));');

%end .. @rdf/gram
