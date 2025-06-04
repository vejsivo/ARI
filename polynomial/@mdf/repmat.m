function B = repmat(A,M,N)
%REPMAT   Replicate and tile matrix-den fraction
%
% B = REPMAT(A,M,N)  or  B = REPMAT(A,[M N])  replicates and tiles 
% the matrix-den fraction A to produce the M-by-N block matrix-den
% fraction B.
%
% REPMAT(A,N) means REPMAT(A,N,N).

%       Author:  J. Jezek  23-Dec-1999
%       Copyright(c) 1999 by Polyx, Ltd.
%       $ Revision $  $ Date 26-Apr-2000 $
%                     $ Date 06-Nov-2000 $
%                     $ Date 30-Sep-2002 $
%                     $ Date 14-Oct-2002 $

global PGLOBAL;

if nargin<2,
   error('Not enough input arguments.');
end;
if ~isa(M,'double'),
   error('Invalid 2nd argument.');
end;

B = 0;
if nargin==2,
   eval('B = mdf(repmat(A.frac.num,M),repmat(A.frac.den,M));', ...
        'error(peel(lasterr));'); N = M;
else
   if ~isa(N,'double'),
      error('Invalid 3rd argument.');
   end;
   eval('B = mdf(repmat(A.frac.num,M,N),repmat(A.frac.den,M,N));', ...
      'error(peel(lasterr));');
end;
  
if M>0 & N>0,
   props(B,A.frac.p,A.frac.tp,A.frac.c,A.frac.tc,A.frac.r);
end;
if strcmp(PGLOBAL.COPRIME,'cop'), B = coprime(B);
end;
if strcmp(PGLOBAL.REDUCE,'red'), B = reduce(B);
else B = smreduce(B);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), B = defract(B);
end;

%end .. @mdf/repmat
