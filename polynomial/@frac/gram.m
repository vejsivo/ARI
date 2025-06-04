function G = gram(F,tol)
%GRAM   Gramian of a stable fraction
%
% The command  GRAM(F)  computes the Gramian of the stable
% fraction F (matrix-den fraction or scalar-den fraction).
%
% A tolerance TOL may be specified as an additional input
% argument. Its default value is the global zeroing tolerance.
%
% See also RDF/GRAM, LDF/GRAM, POL/GRAM.

%    Author: J.Jezek, 17-Jul-2001
%    Copyright(c) 2001 by Polyx, Ltd.
%    $ Revision $  $ Date 25-Jul-2002 $

global PGLOBAL;

if nargin==1 | isempty(tol),
   tol = PGLOBAL.ZEROING;
elseif ~isa(tol,'double'),
   error('Invalid tolerance.');
end;

Fcl = class(F);
if strcmp(Fcl,'frac'),
   error('Invalid 1st argument.');
end;

[s1,s2] = size(F);
if s1>s2, F - rdf(F);
else F = ldf(F);
end;

eval('G = gram(F,tol);', 'error(peel(lasterr));');

%end .. @frac/gram

