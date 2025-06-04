function B = repmat(A,M,N)
%REPMAT   Replicate and tile two-sided polynomial
%
% B = REPMAT(A,M,N)  or  B = REPMAT(A,[M N])  replicates and tiles 
% the tsp matrix A to produce the M-by-N block tsp matrix B.
%
% REPMAT(A,N) means REPMAT(A,N,N).

%       Author:  J. Jezek  23-Dec-1999
%       Copyright(c) 1999 by Polyx, Ltd.
%       $ Revision $   $ Date 24-May-2000 $
%                      $ Date 31-Oct-2000 $

if nargin<2,
   error('Not enough input arguments.');
end;
if ~isa(M,'double'),
   error('Invalid 2nd argument.');
end;

PP = 0;
if nargin==2,
   eval('PP = repmat(A.p,M);','error(peel(lasterr));');
else
   if ~isa(N,'double'),
      error('Invalid 3rd argument.');
   end;
   eval('PP = repmat(A.p,M,N);','error(peel(lasterr));');
end;
B = tsp(PP); B.h = A.h;
B = shift(B,A.o);
  
%end .. @tsp/repmat
