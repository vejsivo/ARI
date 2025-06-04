function [A,B,C,D,E] = abcde(F,tol)
%ABCDE   Descriptor state space realization of a right-den fraction
%
% For right-denominator fraction F, the command
%    [A,B,C,D,E] = ABCDE(F)
% returns the descriptor state space realization A,B,C,D,E,
% that is    F(x) = C*(xE-A)^-1*B + D
% where x is 's','p'(continous-time) or 'z','q' (discrete-time)
% or         F(x) = C*x*(E-A*x)^-1*B + D
% where x is 'z^-1,'d' (discrete time).
% The resulting A,B,C,D,E are numerical matrices.
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
eval('[A,B,C,D,E] = rmf2dss(F.frac.num,F.frac.den,tol);', ...
   'error(peel(lasterr));');

%end .. @rdf/abcde

