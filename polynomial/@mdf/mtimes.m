function R = mtimes(P,Q,tol);
%MTIMES(*)     Matrix multiply matrix-den fractions
%           R = P*Q   or   R = MTIMES(P,Q)
%
% For matrix-den fractions P,Q the command returns their matrix product.
% The number of columns of P must be equal to the number of rows of Q,
% unless one of matrices P,Q is scalar. In such a case, this scalar
% is multiplied with every entry of the second matrix.
%
% The variable symbols of P,Q should be the same. When not, a warning
% is issued and the symbols are changed to the standard one. However, 
% if one symbol is 'z' and the second 'z^-1' then the symbols play
% a role. The resulting symbol is taken form P, no warning being issued.
%
% An optional input argument TOL can specify the zeroing tolerance
% to be used instead of the standard one.
%
% See also MDF/TIMES.

%       Author:  J. Jezek  03-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 26-Apr-2000 $
%                     $ Date 29-May-2000 $
%                     $ Date 02-Aug-2000 $
%                     $ Date 24-Jan-2002 $
%                     $ Date 30-Sep-2002 $
%                     $ Date 14-Oct-2002 $

global PGLOBAL;

ni = nargin;
if ni<2,
   error('Not enough input arguments.');
elseif ni==2 | isempty(tol),
   tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;

eval('P = mdf(P);','error(''Invalid 1st argument.'');');
eval('Q = mdf(Q);','error(''Invalid 2nd argument.'');');

[tv,Rv,P,Q] = testvf(P,Q);
if ~tv, warning('Inconsistent variables.');
end;
[th,Rh,P,Q] = testhf(P,Q,Rv);
if ~th, warning('Inconsistent sampling periods.');
end;

Ps = P.frac.s; Ps1 = Ps(1); Ps2 = Ps(2);
Qs = Q.frac.s; Qs1 = Qs(1); Qs2 = Qs(2);
if Ps2~=Qs1,
   if all(Ps==1) | all(Qs==1),
      R = times(P,Q,tol);
   else
      error('Matrices of inconsistent dimensions.');
   end;
else
   RR = pol(ones(Ps1,Qs2));
   R = mdf(RR,RR);
   Rp = 'prop'; Rc = 'cop'; Rr = 'red';
   Rtp = tol; Rtc = tol;
   for i = 1:Ps1,
      for j=1:Qs2,
         PP = mdf(P.frac.num(i,:),P.frac.den(i,:));
         QQ = mdf(Q.frac.num(:,j),Q.frac.den(:,j));
         SS = mdf(sum(PP.*QQ.',[],tol));
         R.frac.num(i,j) = SS.frac.num;
         R.frac.den(i,j) = SS.frac.den;
         if ~strcmp(Rp,'prop?') & ~strcmp(SS.frac.p,'prop'),
            Rp = SS.frac.p; Rtp = SS.frac.tp;
         end;
         if ~strcmp(Rc,'cop?') & ~strcmp(SS.frac.c,'cop'),
            Rc = SS.frac.c; Rtc = SS.frac.tc;
         end;
         if ~strcmp(Rr,'red?') & ~strcmp(SS.frac.r,'red'),
            Rr = SS.frac.r;
         end;
      end;
   end;
   R.frac.h = Rh;
   R.frac.num.h = Rh; R.frac.den.h = Rh;
   props(R,Rc,Rtc,Rp,Rtp,Rr);
end;

if strcmp(PGLOBAL.COPRIME,'cop'), R = coprime(R,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'), R = reduce(R);
else R = smreduce(R);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), R = defract(R);
end;

%end .. @mdf/mtimes
