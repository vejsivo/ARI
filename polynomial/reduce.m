function G = reduce(F,tol)
%REDUCE    Reduce constant
%
% The command  G = REDUCE(F)  or  G = REDUCE(F,TOL)
% returns  G = F .
%
% This macro exists only for completeness.
% See also LDF/REDUCE, RDF/REDUCE, MDF/REDUCE, SDF/REDUCE.

%      Author:  J. Jezek  27-Apr-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 30-Sep-2002 $

if nargin==0,
   error('Not enough input arguments.');
end;
if ~isa(F,'double') | ndims(F)>2,
   error('Invalid 1st argument.');
end;

G = F;

%end .. reduce
