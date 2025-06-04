function G = deriv(F,arg2,arg3,arg4)
%DERIV    Derivative of constant
%
% The command  G = DERIV(F)
% computes the derivative of constant F (standard Matlab matrix).
% The result is zero; this macro exists only for completeness.
%
% The command  G = DERIV(F,N)
% computes the N-th derivative.
%
% See also POL/DERIV, TSP/DERIV, RDF/DERIV, LDF/DERIV,
% SDF/DERIV, MDF/DERIV.

%       Author:  J. Jezek  13-Apr-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 15-Jun-2000 $

if nargin==0,
   error('Not enough input arguments.');
end;
if ~isa(F,'double') | ndims(F)>2,
   error('Invalid 1st argument.');
end;

G = zeros(size(F));

%end .. deriv
