function Fr = reverse(F);
%REVERSE  Reverse variable of right-den fraction
%
% Given the right-den fraction F(VAR), the command
%    FR = REVERSE(F)
% returns a right-den fraction FR such that  FR(VAR) = F(1/VAR) .
%
% If the variable of F is 'p','s','q',or 'd', so is the variable of FR.
% However, if it is 'z' or 'z^-1', it is changed to 'z^-1' or 'z'.
% In such a case, it holds  F==FR  and if F is known to be (column)
% reduced, so is FR.

%       Author:  J. Jezek  28-Dec-1999
%       Copyright(c) 1999 by Polyx, Ltd.
%       $ Revision $  $ Date 25-Apr-2000 $
%                     $ Date 06-Feb-2001 $
%                     $ Date 30-Sep-2002 $
%                     $ Date 14-Oct-2002 $
%                     $ Date 22-Nov-2002  bug: coprime  $

Fv = F.frac.v;
if strcmp(Fv,'z'), symb = 'z^-1';
elseif strcmp(Fv,'z^-1'), symb ='z';
else symb = Fv;
end;
var = pol([0 1],1,symb);

Fs1 = F.frac.s(1); Fs2 = F.frac.s(2);
%DN = [F.frac.d;F.frac.n];
DN = vertcat(F.frac.den,F.frac.num);
DDNN = repmat(var,size(DN));
if strcmp(F.frac.r,'red'),
   degcolDN = deg(DN,'col');
   for i = 1:Fs2,
      DDNN(:,i) = rev(DN(:,i),degcolDN(i),var);
   end;
else
   degDN = deg(DN);
   DDNN = rev(DN,degDN,var);
end;
DD = DDNN(1:Fs2,:); NN = DDNN(Fs2+1:Fs2+Fs1,:);
Fr = rdf(NN,DD);
Fr.frac.h = F.frac.h;

if strcmp(symb,'z') | strcmp(symb,'z^-1'),
   props(Fr,F.frac.r,F.frac.p,F.frac.tp);
end;

%end .. @rdf/reverse
