function Rs = prod(R,dim,tol);
%PROD   Element-wise product of right-den or left-den fraction
%
% For vectors,  PROD(X) is the product of the elements of X. 
% For matrices, PROD(X) is a row vector with the product over each column. 
%
% PROD(X,DIM) works along the dimension DIM.
%
% PROD(X,DIM,TOL) or PROD(X,[],TOL) works with zeroing specified by the
% input tolerance TOL.

%       Author:  J. Jezek 22-Feb-2003
%       Copyright (c) 2003 by Polyx, Ltd.
%       $ Revision $  $ Date 28-Feb-2003 $

global PGLOBAL;

Rcl = class(R);
if strcmp(Rcl,'frac'),
   error('Invalid 1st argument.');
end;

ni = nargin;
if ni<3,
   tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
end;
if ni<2,
   dim = [];
else
   if ~isa(dim,'double'),
      error('Invalid dimension.');
   end;
end;

Rm = mdf(R);
Rs = 0;
eval('Rs = prod(Rm,dim,tol);','error(peel(lasterr));');
if isa(R,'rdf'), Rs = rdf(Rs);
elseif isa(R,'ldf'), Rs = ldf(Rs);
end;

%end .. @frac/prod

