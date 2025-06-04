function [Q,R] = laurent(F,degree,T,tau,tol);
%LAURENT      Laurent series of right-den fraction
%
% If F is a right-den fraction in variables 'z' or 'q' then
% the command   Q = LAURENT(F,DEGREE)
% returns (possibly two-sided) polynomial Q ... the first
% terms of the Laurent series of F about the point z=q=Inf
% (from positive powers of 'z' down to the power  z^-DEGREE ).
%
% If F is a right-den fraction in variables 'z^-1' or 'd' then
% the command   Q = LAURENT(F,DEGREE)
% returns (possibly two-sided) polynomial Q ... the first
% terms of the Laurent series of F about the point z^-1=d=0
% (from positive powers of 'z' down to the power  z^-DEGREE ).
%
% If F is a right-den fraction in variables 's' or 'p' then
% the command   Q = LAURENT(F,DEGREE)  or  Q = LAURENT(F)
% returns polynomial Q ... nonnegative powers of the Laurent
% series of F about the point s = p = Inf.
% The command  [Q,R] = LAURENT(F,DEGREE)
% returns polynomial Q ... positive powers, and  R ... three-
% dimensional array R(:,:,DEGREE+1) of nonpositive powers
%      R(s) = R0 + R1*s^-1 + ... + RDEGREE*s^-DEGREE
% with  RI = R(:,:,I+1)  for  I = 0,1,...DEGREE .
%
% The default for DEGREE is  MAX(DEG(F.NUM),DEG(F.DEN)) ,
%
% For discrete-time fraction F, optional arguments T,TAU may
% be used, defaults being  T = F.h, TAU = 0 . The resulting
% Q.h is set to T. If TAU is vector, the resulting Q is a
% cell vector, all cells containing the same. For continuous
% time fractions F, the arguments T,TAU, if given, are not used.
%
% A tolerance TOL may be specified as an optional (fifth)
% input argument.

%        Author:  J.Jezek  16-Feb-2000
%        Copyright(c) 2000 by Polyx, Ltd.
%        $ Revision $  $ Date 18-Jul-2000 $
%                      $ Date 16-Jan-2001 $
%                      $ Date 24-Sep-2002 $
%                      $ Date 14-Oct-2002 $
%                      $ Date 28-Feb-2003 $

global PGLOBAL;

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;

if ~isa(F,'rdf'),
   error('Invalid 1st argument.');
end;
Fdegree = max(deg(F.frac.num),deg(F.frac.den));
Fh = F.frac.h; Ftau = 0;
Ftol = PGLOBAL.ZEROING;

if ni>=2 & ~isempty(degree),
   if ~isa(degree,'double') | length(degree)~=1 | ~isfinite(degree) |...
         ~isreal(degree) | degree<0 | floor(degree)~=degree,
      error('Invalid degree.');
   end;
   Fdegree = degree;
end;

if ni>=3 & ~isempty(T),
   if ~isa(T,'double') | length(T)~=1 | ~isreal(T) | T<0,
      error('Invalid sampling period.');
   end;
   Fh = T;
end;

if ni>=4 & ~isempty(tau),
   if ~isa(tau,'double') | ndims(tau)~=2 | all(size(tau)~=1) | ...
         ~isreal(tau) | any(tau<0),
      error('Invalid sampling phase.');
   end;
   Ftau = tau;
end;

if ni==5 & ~isempty(tol),
   if isa(tol,'char'),
      tol = str2num(tol);
      if isempty(tol),
         error('Invalid tolerance.');
      end;
   end;
   if ~isa(tol,'double') | length(tol)~=1 | ~isreal(tol) | ...
         tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
   Ftol = tol;
end;

if isempty(F.frac.v) | strcmp(F.frac.v,'s') | strcmp(F.frac.v,'p'),
   if nargout<=1,
      Q = longrdiv(F.frac.num,F.frac.den,Fdegree,Ftol);
   else
      [Q,R] = longrdiv(F.frac.num,F.frac.den,Fdegree,Ftol);
      if isa(R,'pol'),
         RR = R{0};
         [sR1,sR2] = size(RR);
         R = zeros(sR1,sR2,Fdegree+1);
         R(:,:,1) = RR;
      end;      
   end;
   
else
   if nargout>1,
      error('Too many output arguments.');
   end;
   
   [QQ,RR] = longrdiv(F.frac.num,F.frac.den,Fdegree,Ftol);
   if QQ==0, QQ = pol(RR);
   elseif deg(RR)<=0,
      QQ = QQ+RR;
   else QQ.v = 'z'; RR.v = 'z^-1';
      QQ = tsp(QQ+RR);
   end;
   QQ.h = Fh;
   
   ltau = length(Ftau);
   if ltau==1,
      Q = QQ;
   else
      Q = cell(size(Ftau));
      for k = 1:ltau,
         Q{k} = QQ;
      end;
   end;
end;

%end .. @rdf/laurent
