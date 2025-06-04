function A = sum(R,dim,tol);
%SUM   Element-wise sum of matrix-den fraction
% 
% For a matrix-den vector fraction R (i.e., a matrix-den fraction
% with one row or one column), SUM(R) is the sum of the elements of R. 
%
% For other matrix-den fractions, SUM(R) is a row vector with the sum 
% over each column. 
%
% SUM(R,DIM) sums along the dimension DIM, where DIM is 1 or 2.
%
% SUM(R,DIM,TOL) or SUM(R,[],TOL) works with zeroing specified by 
% the input tolerance TOL.

%      Author:  J. Jezek  03-Jan-2000
%      Copyright (c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 26-Apr-2000 $
%                    $ Date 25-May-2001 $
%                    $ Date 30-Sep-2002 $
%                    $ Date 14-Oct-2002 $
%                    $ Date 28-Feb-2003 $

global PGLOBAL;

ni = nargin;
if ni==1,
  dim = []; tol = PGLOBAL.ZEROING;
elseif ni==2 | isempty(tol),
  tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;
if ~isempty(dim) & (~isa(dim,'double') | (dim~=1 & dim~=2)),
   error('Invalid dimension.');
end;

Rtp = R.frac.tp; Rtc = R.frac.tc;
Rs = R.frac.s; Rs1 = Rs(1); Rs2 = Rs(2);
Rd = R.frac.num.d;
if isempty(dim),
   if isempty(Rd),
      A = mdf(sum(zeros(Rs)));
      return;
   end;
   dim = find(Rs~=1);
   if isempty(dim),
      dim = 1;
   else
      dim = dim(1);
   end;
end;

if Rs1==0 | Rs2==0,
   A = mdf(sum(zeros(Rs),dim));
   return;
end;
if Rs1==1 &  dim==2,
   n = Rs2;
elseif Rs2==1 & dim==1,
   n = Rs1;
else
   if dim==1,
      A = mdf(R.frac.num(1,:),R.frac.den(1,:));
      props(A,'prop',Rtp,'cop',Rtc,'red');
      for j = 1:Rs2,
         RR = sum(mdf(R.frac.num(:,j),R.frac.den(:,j)),[],tol);
         A.frac.num(j) = RR.frac.num; A.frac.den(j) = RR.frac.den; 
         [RR,A] = logic(RR,A);
      end;
   else
      A = mdf(R.frac.num(:,1),R.frac.den(:,1));
      props(A,'prop',Rtp,'cop',Rtc,'red');
      for i = 1:Rs1,
         RR = sum(mdf(R.frac.num(i,:),R.frac.den(i,:)),[],tol);
         A.frac.num(i) = RR.frac.num; A.frac.den(i) = RR.frac.den;
         [RR,A] = logic(RR,A);
      end;
   end;
   return;
end;

[D,DD] = plcm(R.frac.den,[],tol);
A = mdf(sum(R.frac.num.*DD,[],tol),D);

if strcmp(R.frac.p,'prop'),
   props(A,'prop',Rtp);
end;
if strcmp(PGLOBAL.COPRIME,'cop'),
   A = coprime(A,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   A = reduce(A);
else 
   A = smreduce(A);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'),
   A = defract(A);
end;

function [RR,A] = logic(RR,A)

if strcmp(RR.frac.p,'prop?'),
   props(A,'prop?',[]);
end;
if ~strcmp(A.frac.p,'prop?'),
   if strcmp(RR.frac.p,'nprop'),
      props(A,'nprop',Rtp);
   end;
end;

if strcmp(RR.frac.c,'cop?'),
   props(A,'cop?',[]);
end;
if ~strcmp(A.frac.c,'cop?'),
   if strcmp(RR.frac.c,'ncop'),
      props(A,'ncop',Rtc);
   end;
end;

if strcmp(RR.frac.r,'red?'),
   props(A,'red?');
end;
if ~strcmp(A.frac.r,'red?'),
   if strcmp(RR.frac.r,'nred'),
      props(A,'nred');
   end;
end;

%end .. logic

%end .. @mdf/sum
