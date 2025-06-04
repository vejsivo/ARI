function G = smreduce(F)
%SMREDUCE    Small reduce left-den fraction
%
% For left-den fraction F, the command
%   G = SMREDUCE(F)  returns left-den fraction G,
% which is equal to F but modified:
% if the variable symbol is 'z','q','s' or 'p'
%    then the leading coefficient matrix of the denominator
%    is identity matrix,
% if the variable symbol is 'zi' or 'd'
%    then the trailing coefficient matrix of the denominator
%    is identity matrix.
%
% This modification takes place only when the leading/trailing
% matrix is nonsingular.
%
% See also LDF/REDUCE.

%     Author:  J. Jezek, 30-Sep-2002
%     Copyright(c) 2002 by Polyx, Ltd.
%     $ Revision $  $ Date 14-Oct-2002 $

G = F;
if strcmp(F.frac.r,'red') | isempty(F), return;   % quick exit
end;

if strcmp(F.frac.v,'z^-1') | strcmp(F.frac.v,'d'),
   M = tcoef(F.frac.den);
else
   M = lcoef(F.frac.den);
end;

[n,m] = size(F);
if rank(M)==n,
   DN = M\horzcat(F.frac.den,F.frac.num);
   G.frac.den = DN(1:n,1:n);
   G.frac.num = DN(1:n,n+1:n+m);
end;

%end .. @ldf/smreduce


