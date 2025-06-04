function Rs = prod(R,dim,tol);
%PROD   Element-wise product of scalar-den fraction
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

eval('R = sdf(R);','error(peel(lasterr));');

ni = nargin;
if ni==1,
  dim = []; tol = PGLOBAL.ZEROING;
elseif ni==2 | isempty(tol),
  tol = PGLOBAL.ZEROING;
else
  if ~isa(tol,'double'),
     error('Invalid tolerance.');
  end;
end;

Rsi = R.frac.s;
Rd = R.frac.num.d;
if isempty(dim),
   if isempty(Rd),
      Rs = sdf(prod(zeros(Rsi)));
      return;
   end;
   dim = find(Rsi~=1);
   if isempty(dim),
      dim = 1;
   else
      dim = dim(1);
   end;
end;   
if ~isa(dim,'double') | length(dim)~=1 | (dim~=1 & dim~=2),
   error('Invalid dimension.');
end;
if isempty(Rd),
   Rs = sdf(prod(zeros(Rsi),dim));
   return;
end;

Rs = 0;
eval('Rs = sdf(prod(R.frac.num,dim,tol),power(R.frac.den,Rsi(dim),tol));', ...
   'error(peel(lasterr));');

if strcmp(R.frac.p,'prop'),
   props(Rs,'prop');
end;
if strcmp(R.frac.r,'red'),
   props(Rs,'red');
end;
if strcmp(PGLOBAL.COPRIME,'cop'), Rs = coprime(Rs,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'), Rs = reduce(Rs);
else Rs = smreduce(Rs);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), Rs = defract(Rs);
end;

%end .. @sdf/prod

