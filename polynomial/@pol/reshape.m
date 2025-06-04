function Q = reshape(P,M,N)
%RESHAPE        Change size of polynomial
%
% For polynomial matrix P, the command  Q = RESHAPE(P,M,N)  returns
% a M-by-N polynomial matrix whose elements are taken columnwise
% from P. An error results if P has not M*N elements.
%
% Q = RESHAPE(P,[M N])  is the same thing.

%      Author:  J. Jezek  05-Jan-2000
%      Copyright(c) by Polyx, Ltd.
%      Revision   $ Date $   $ 25-Feb-2000  $
%                 $ Date $   $ 06-Feb-2002  $

eval('P = pol(P);','error(peel(lasterr));');
d1 = P.d+1;
if isempty(d1) | d1<0, d1 = 0;
end;

ni = nargin;
if ni==1,
   error('Not enough input arguments.');
end;
if ~isa(M,'double'),
   error('Invalid 2nd argument.');
end;
if ni==2,
   if length(M)~=2,
      error('Invalid 2nd argument or missing 3rd one.');
   end;
   if M(1)*M(2)~=P.s(1)*P.s(2),
      error('To RESHAPE the number of elements must not change.');
   end;
   eval('P.c = reshape(P.c,M(1),M(2),d1);','error(peelf(lasterr));');
   P.s = M;
else
   if ~isa(N,'double'),
      error('Invalid 3rd argument.');
   end;
   if M*N~=P.s(1)*P.s(2),
      error('To RESHAPE the number of elements must not change.');
   end;
   eval('P.c = reshape(P.c,M,N,d1);','error(peelf(lasterr));');
   P.s = [M N];
end;
Q = P;
Q.u = [];
Q.version = 3.0;

%end .. @pol/reshape
