function [A,B,C,D] = abcd(F,tol)
%ABCD   State space realization of a right-den fraction
%
% For right-denominator fraction F, the command
%    [A,B,C,D] = ABCD(F)
% returns the (generalized) contoller form realization A,B,C,D,
% that is    F(x) = C*(xI-A)^-1*B + D(x)
% where x is 's','p' (continuous time) or 'z','q' (discrete time)
% or         F(x) = C*x*(I-A*x)^-1*B + D(x^-1)
% where x is 'z^-1','d' (discrete time).
% The resulting A,B,C are numerical matrices, D is a numerical
% matrix or a polynomial in 's', 'p', 'z' or 'q'.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% For backward conversion, see RDF/RDF.

%     Author: J. Jezek, 10-Mar-2000
%     Copyright (c) 2000 by Polyx, Ltd.
%     $ Revision $  $ Date 25-Apr-2000 $
%                   $ Date 12-Oct-2000 $
%                   $ Date 14-Oct-2002 $

global PGLOBAL;

if nargin==1,
   tol = PGLOBAL.ZEROING;
else
   if ~isempty(tol) & ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
end;

eval('F = reduce(F,tol);','error(peel(lasterr));');
eval('[A,B,C,D] = rmf2ss(F.frac.num,F.frac.den,tol);', ...
   'error(peel(lasterr));');

%end .. @rdf/abcd

