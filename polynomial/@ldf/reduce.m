function G = reduce(F,tol)
%REDUCE   Reduce left-den fraction
%
% The command  G = REDUCE(F) brings the fraction  F = D^-1*N
% to the form with denominator row reduced and in the 
% row echelon form.
%
% If the variable symbol is 'z','q','s' or 'p'
%    then the row leading matrix plays a role,
% if the variable symbol is 'z^-1' or 'd'
%    then the row trailing matrix.
%
% The reduced form is unique among all left-denominator fractions
% that are equal each to other and left coprime.
%
% See also ROWRED, ECHELON, LDF/COPRIME.

%        Author:  J. Jezek  07-Jan-2000
%        Copyright(c) 2000 by Polyx, Ltd.
%        $ Revision $   $ Date 21-Apr-2000 $
%                       $ Date 06-Feb-2001 $
%                       $ Date 06-Oct-2002 $
%                       $ Date 14-Oct-2002 $

global PGLOBAL;

if nargin==1 | isempty(tol),
   tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;

if strcmp(F.frac.r,'red'),    % quick exit
   G = F; return;
end;

if isempty(F.frac.den),
   G = F; props(G,'red'); return;
end;

Fc = F.frac.c; Ftc = F.frac.tc;
Fp = F.frac.p; Ftp = F.frac.tp;
if strcmp(F.frac.v,'z^-1') | strcmp(F.frac.v,'d'),
   F = reverse(F); rev = 1;
else
   rev = 0;
end;

[D,rk,U] = rowred(F.frac.den,tol);
N = mtimes(U,F.frac.num,tol);
[D,U] = echelon(D,'row',tol);
N = mtimes(U,N,tol);
G = ldf(D,N);
G.frac.h = F.frac.h;

props(G,'red');
if rev, G = reverse(G);
end;
props(G,'red',Fc,Ftc,Fp,Ftp);

%end .. @ldf/reduce
