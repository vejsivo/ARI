function B = repmat(A,M,N)
%REPMAT   Replicate and tile left-den fraction
%
% B = REPMAT(A,M,N)  or  B = REPMAT(A,[M N])  replicates and tiles 
% the left-den fraction A to produce the M-by-N block left-den fraction B.
%
% REPMAT(A,N) means REPMAT(A,N,N).

%       Author:  J. Jezek  04-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 21-Apr-2000 $
%                     $ Date 02-Nov-2000 $
%                     $ Date 30-Sep-2002 $
%                     $ Date 14-Oct-2002 $

global PGLOBAL;

if nargin<2,
   error('Not enough input arguments.');
end;
if ~isa(M,'double'),
   error('Invalid 2nd argument.');
end;

NN = 0;
if nargin==2,
   eval('NN = repmat(A.frac.num,M);','error(peel(lasterr));');
   if length(M)==1, N = M;
   else N = M(2); M = M(1);
   end;
else
   if ~isa(N,'double'),
      error('Invalid 3rd argument.');
   end;
   eval('NN = repmat(A.frac.num,M,N);','error(peel(lasterr));');
end;

CC = cell(1,M);
for i = 1:M,
   CC{i} = A.frac.den;
end;
DD = blkdiag(CC{1:M});
B = ldf(DD,NN);

if M>0 & N>0,
   props(B,A.frac.p,A.frac.tp);
end;
if strcmp(PGLOBAL.COPRIME,'cop'), B = coprime(B);
end;
if strcmp(PGLOBAL.REDUCE,'red'), B = reduce(B);
else B = smreduce(B);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), B = defract(B);
end;

%end .. @ldf/repmat
