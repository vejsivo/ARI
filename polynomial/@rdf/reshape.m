function Q = reshape(P,M,N)
%RESHAPE        Change size of a right-den fraction
%
% For right-den fraction P, the command Q = RESHAPE(P,M,N)
% returns a M-by-N right-den fraction whose elements are taken
% columnwise from P. An error results if P has not M*N
% elements.
%
% Q = RESHAPE(P,[M N])  is the same thing.

%      Author:  J. Jezek  08-Feb-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 25-Apr-2000 $
%                    $ Date 01-Nov-2000 $
%                    $ Date 30-Sep-2002 $

global PGLOBAL;

eval('P = rdf(P);','error(peel(lasterr));');
P = sdf(P);

ni = nargin;
if ni==1,
   error('Not enough input arguments.');
end;
if ~isa(M,'double'),
   error('Invalid 2nd argument.');
end;
Q = 0;
if ni==2,
   eval('Q = reshape(P,M);','error(peel(lasterr));');
else
   if ~isa(N,'double'),
      error('Invalid 3rd argument.');
   end;
   eval('Q = reshape(P,M,N);','error(peel(lasterr));')
end;

Q = rdf(Q);

props(Q,P.p,P.tp);
if strcmp(PGLOBAL.COPRIME,'cop'), Q = coprime(Q);
end;
if strcmp(PGLOBAL.REDUCE,'red'), Q = reduce(Q);
else Q = smreduce(Q);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), Q = defract(Q);
end;

%end .. @rdf/reshape
