function Q = reduce(M)
%SMREDUCE    Small reduce scalar-den fraction
%
% For scalar-den fraction M, the command  Q = SMREDUCE(M) returns
% scalar-den fraction Q, which is equal to M but modified:
% if the variable symbol is 'z','q','s' or 'p'
%    then the leading coefficient of the denominator is 1
% if the variable symbol is 'z^-1' or 'd'
%    then the trailing coefficient of the denominator is 1.
%
% This macro exists only for completeness.
% See also SDF/REDUCE.

%        Author:  J. Jezek  30-Sep-2002
%        Copyright(c) 2002 by Polyx, Ltd.
%        $ Revision $  $ Date 14=Oct-2002 $

if strcmp(M.frac.r,'red'),   % quick exit
   Q = M; return;
end;

Mc = M.frac.c; Mp = M.frac.p;
Mtp = M.frac.tp; Mtc = M.frac.tc;
if strcmp(M.frac.v,'z^-1') | strcmp(M.frac.v,'d'),
   MM = 1/tcoef(M.frac.den);
else
   MM = 1/lcoef(M.frac.den);
end;
Q = sdf(M.frac.num*MM, M.frac.den*MM);
Q.frac.h = M.frac.h;
props(Q,'red',Mc,Mtc,Mp,Mtp);

%end .. @sdf/smreduce
