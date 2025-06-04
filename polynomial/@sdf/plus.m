function R = plus(P,Q,tol);
%PLUS(+)     Add scalar-den fractions
%           R = P+Q   or   R = PLUS(P,Q)
%
% For scalar-den fractions P,Q, the command returns scalar-den fraction R
% whose (i,j)-th entry is  P(i,j)+Q(i,j) .
%
% Matrices P,Q must have the same dimensions, unless one of them is
% scalar. In such a case, this scalar is added with every entry
% of the second matrix.
%
% The variable symbols of P,Q should be the same. When not, a warning
% is issued and the symbols are changed to the standard one. However, 
% if one symbol is 'z' and the second 'z^-1' then the symbols play
% a role. The resulting symbol is taken form P, no warning being issued.
%
% An optional input argument TOL can specify the zeroing tolerance
% to be used instead of the standard one.
%
% See also FRAC/UPLUS, FRAC/MINUS.

%       Author:  J. Jezek  26-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 26-Apr-2000 $
%                     $ Date 29-May-2000 $
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
   if ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
end;

eval('P = sdf(P);','error(''Invalid 1st argument.'');');
eval('Q = sdf(Q);','error(''Invalid 2nd argument.'');');

[tv,Rv,P,Q] = testvf(P,Q);
if ~tv, warning('Inconsistent variables.');
end;
[th,Rh,P,Q] = testhf(P,Q,Rv);
if ~th, warning('Inconsistent sampling periods.');
end;
[td,P,Q] = testdf(P,Q);
if ~td, error('Matrices not of the same dimensions.');
end;

[X Y] = axby0(P.frac.den,-Q.frac.den,tol);
NN = plus(times(P.frac.num,X,tol),times(Q.frac.num,Y,tol),tol);
DD = times(P.frac.den,X,tol);
R = sdf(NN,DD);

if strcmp(P.frac.p,'prop') & strcmp(Q.frac.p,'prop') & ...
      P.frac.tp==Q.frac.tp,
   props(R,'prop',P.frac.tp);
end;
if strcmp(PGLOBAL.COPRIME,'cop'),
   R = coprime(R,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   R = reduce(R);
else
   R = smreduce(R);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'),
   R = defract(R);
end;

%end .. @sdf/plus
