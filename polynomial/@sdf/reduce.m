function Q = reduce(M,tol)
%REDUCE    Reduce scalar-den fraction
%
% For scalar-den fraction M, the command  Q = REDUCE(M) returns
% scalar-den fraction Q, which is equal to M but modified:
% if the variable symbol is 'z','q','s' or 'p'
%    then the leading coefficient of the denominator is 1
% if the variable symbol is 'z^-1' or 'd'
%    then the trailing coefficient of the denominator is 1.
%
% The reduced form is unique among all scalar-den fractions
% that are equal each to other and coprime.
%
% Possible optional input argument TOL is ignored.
%
% See also SDF/COPRIME.

%        Author:  J. Jezek  24-Jan-2000
%        Copyright(c) 2000 by Polyx, Ltd.
%        $ Revision $  $ Date 26-Apr-2000 $
%                      $ Date 06-Feb-2001 $
%                      $ Date 06-Oct-2002 $
%                      $ Date 14-Oct-2002 $

if nargin==2 & ~isa(tol,'double'),
      error('Invalid tolerance.');
end;      % possible TOL argument ignored
   
if strcmp(M.frac.r,'red'),   % quick exit
   Q = M; return;
end;

Mc = M.frac.c; Mtc = M.frac.tc;
Mp = M.frac.p; Mtp = M.frac.tp;
if strcmp(M.frac.v,'z^-1') | strcmp(M.frac.v,'d'),
   MM = 1/tcoef(M.frac.den);
else
   MM = 1/lcoef(M.frac.den);
end;
Q = sdf(M.frac.num*MM, M.frac.den*MM);
Q.frac.h = M.frac.h;
props(Q,'red',Mc,Mtc,Mp,Mtp);

%end .. @sdf/reduce
