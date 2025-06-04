function B = repmat(A,M,N)
%REPMAT   Replicate and tile polynomial
%
% B = REPMAT(A,M,N)  or  B = REPMAT(A,[M N])  replicates and tiles 
% the polynomial matrix A to produce the M-by-N block polynomial matrix B.
%
% REPMAT(A,N) means REPMAT(A,N,N).

%       Author:  J. Jezek  23-Dec-1999
%       Copyright(c) 1999 by Polyx, Ltd.

if nargin<2,
   error('Not enough input arguments.');
end;
lM = length(M);
if ~isa(M,'double') | lM<1 | lM>2 | ~isreal(M) | ...
      any(M<0) | any(floor(M)~=M),
   error('Invalid 2nd argument.');
end;
if nargin==2,
   if lM==1,
      MN = [M M];
   else
      MN = M;
   end;
else
   if ~isa(N,'double') | length(N)~=1 | ~isreal(N) | ...
         N<0 | floor(N)~=N | lM>1,
      error('Invalid 3rd argument.');
   end;
   MN = [M N];
end;

if all(MN~=0), B.d = A.d;
else B.d = [];
end
B.s = A.s .* MN;
B.c = repmat(A.c,MN);
B.v = A.v;
B.h = A.h;
B.u = [];
B.version = 3.0;

B = class(B,'pol');

%end .. @pol/repmat
