function [A,B,C,D,E] = abcde(F,tol)
%ABCDE   Descriptor state space realization of scalar-den fraction
%
% For scalar-denominator fraction F, the command
%    [A,B,C,D,E] = ABCDE(F)
% returns the (generalized) controller or observer form
% realization A,B,C,D, that is
%    F(x) = C*(xE-A)^-1*B + D
% where x is 's','p' (continous-time) or 'z','q' (discrete-time),
% or F(x) = C*x*(E-A*x)^-1*B + D
% where x is 'z^-1','d' (discrete time).
% The resulting A,B,C,D,E are numerical matrices.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% For backward conversion, se SDF/SDF.

%     Author: J. Jezek, 10-Mar-2000
%     Copyright (c) 2000 by Polyx, Ltd.
%     $ Revision $  $ Date 26-Apr-2000 $
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

if size(F,1)<size(F,2),
   eval('F = ldf(F,tol);','error(peel(lasterr));');
   [A,B,C,D,E] = lmf2dss(F.num,F.den,tol);
else
   eval('F = rdf(F,tol);','error(peel(lasterr));');
   [A,B,C,D,E] = rmf2dss(F.num,F.den,tol);
end;

%end .. @sdf/abcde
