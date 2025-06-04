function Ts = prod(T,dim,tol)
%PROD   Element-wise product of two-sided polynomial
%
% For  vectors, PROD(X) is the product of the elements of X. 
% For matrices, PROD(X) is a row vector with the product over each column. 
%
% PROD(X,DIM) works along the dimension DIM.
%
% PROD(X,DIM,TOL) or PROD(X,[],TOL) works with zeroing specified by the
% input tolerance TOL .

%       Author:  J.Jezek  11-8-99
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 13-Oct-1999  12:00:00  $
%                         $Date: 24-May-2000  11:30:00  $
%                         $Date: 28-Feb-2003            $

global PGLOBAL;

eval('T = tsp(T);','error(peel(lasterr));');

na = nargin;
if na<2, dim = []; end;
if na<3 |isempty(tol), tol = PGLOBAL.ZEROING; end;
if ~isa(tol,'double') | length(tol)~=1 | ...
      ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;

Tsi = T.s;
if isempty(dim),
   if Tsi(1)==0 | Tsi(2)==0,
      Ts = tsp(prod(zeros(Tsi)));
      return;
   end;
   dim = find(Tsi~=1);
   if isempty(dim),
      dim = 1;
   else
      dim = dim(1);
   end;
end;
if ~isa(dim,'double') | length(dim)~=1 | (dim~=1 & dim~=2),
   error('Invalid dimension.');
end;

Ps = prod(T.p,dim,tol);
Ts = tsp(Ps); Ts.h = T.h;
Ts = shift(Ts,T.o*Tsi(dim));

%end .. @tsp/prod

