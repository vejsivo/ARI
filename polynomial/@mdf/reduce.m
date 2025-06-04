function Q = reduce(R,tol)
%REDUCE    Reduce matrix-den fraction
%
% For matrix-den fraction R, the command  Q = REDUCE(R)  returns
% matrix-den fraction Q, which is equal to R but modified:
% if the variable symbol is 'z','q','s' or 'p'
%    then all leading coefficients of the denominators are 1
% if the variable symbol is 'z^-1' or 'd'
%    then all trailing coefficients of the denominators are 1.
%
% The reduced form is unique among all matrix-den fractions that
% are equal each to other and coprime.
%
% Possible optional input argument TOL is ignored.
%
% See also MDF/COPRIME.

%        Author:  J. Jezek  07-Jan-2000
%        Copyright(c) 2000 by Polyx, Ltd.
%        $ Revision $  $ Date 26-Apr-2000 $
%                      $ Date 06-Oct-2002 $
%                      $ Date 14=Oct-2002 $

if nargin==2,
   if ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
end;      % possible TOL argument ignored

Q = R;
if strcmp(Q.frac.r,'red'),
   return;
end;

Rc = R.frac.c; Rtc = R.frac.tc;
Rp = R.frac.p; Rtp = R.frac.tp;
if ~isempty(Q),
   if strcmp(Q.frac.v,'z^-1') | strcmp(Q.frac.v,'d'),
      M = tcoef(Q.frac.den,'ent');
   else
      M = lcoef(Q.frac.den,'ent');
   end;
   Q.frac.num = Q.frac.num./M; Q.frac.den = Q.frac.den./M;
end;
props(Q,'red',Rc,Rtc,Rp,Rtp);

%end .. @mdf/reduce
