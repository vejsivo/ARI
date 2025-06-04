function Fr = reverse(F);
%REVERSE  Reverse variable of matrix-den fraction
%
% Given matrix-den fraction F(VAR), the command
%    FR = REVERSE(F)
% returns matrix-den fraction FR such that  FR(VAR) = F(1/VAR) .
%
% If the variable of F is 'p','s','q',or 'd', so is the variable of FR.
% However, if it is 'z' or 'z^-1', it is changed to 'z^-1' or 'z'.
% So, in the case of 'z' or 'z^-1', it holds  FR==F .

%       Author:  J. Jezek  28-Dec-1999
%       Copyright(c) 1999 by Polyx, Ltd.
%       $ Revision $  $ Date 26-Apr-2000 $
%                     $ Date 12-Jul-2000 $
%                     $ Date 06-Feb-2001 $
%                     $ Date 03-Jan-2002 $
%                     $ Date 14-Oct-2002 $
%                     $ Date 22-Nov-2002  bug:coprime  $
%                     $ Date 28-Feb-2003  warnings  $

Fv = F.frac.v;
if strcmp(Fv,'z'), symb = 'z^-1';
elseif strcmp(Fv,'z^-1'), symb ='z';
else symb = Fv;
end;
var = pol([0 1],1,symb);

N = F.frac.num; D = F.frac.den;
degDN = max(deg(N,'ent'),deg(D,'ent'));
Fs1 = F.frac.s(1); Fs2 = F.frac.s(2);
saved_w = warning; warning off;
NN = repmat(var,size(N)); DD = repmat(var,size(D));
for i = 1:Fs1,
   for j = 1:Fs2,
      NN(i,j) = rev(N(i,j),degDN(i,j),var);
      DD(i,j) = rev(D(i,j),degDN(i,j),var);
   end;
end;
warning(saved_w);
Fr = mdf(NN,DD);
h = F.frac.h; Fr.frac.h = h;
Fr.frac.num.h = h; Fr.frac.den.h = h;

if strcmp(symb,'z') | strcmp(symb,'z^-1'),
   props(Fr,F.frac.p,F.frac.tp,F.frac.r);
end;

%end .. @mdf/reverse
