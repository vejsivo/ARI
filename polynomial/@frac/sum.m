function Rs = sum(R,dim,tol);
%SUM   Element-wise sum of right-den or left-den fraction
%
% For vectors,  SUM(X) is the sum of the elements of X. 
% For matrices, SUM(X) is a row vector with the sum over each column. 
%
% SUM(X,DIM) works along the dimension DIM.
%
% SUM(X,DIM,TOL) or SUM(X,[],TOL) works with zeroing specified by the
% input tolerance TOL.

%       Author:  J. Jezek 22-Feb-2003
%       Copyright (c) 2003 by Polyx, Ltd.

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
eval('Rs = sum(Rm,dim,tol);','error(peel(lasterr));');
if isa(R,'rdf'), Rs = rdf(Rs);
elseif isa(R,'ldf'), Rs = ldf(Rs);
end;

%end .. @frac/sum
