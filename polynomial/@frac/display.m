function display(F)
%DISPLAY   Display fraction

%    Author: J. Jezek, 2000
%    Copyright(c) 2000 by Polyx, Ltd.
%    $ Revision $  $ Date 14-Oct-2002 $

disp('Fraction');
numerator = F.num, denominator = F.den

if strcmp(F.c,'cop'), disp('coprime');
elseif strcmp(F.c,'ncop'), disp('noncoprime');
end;
if strcmp(F.r,'red'), disp('reduced');
elseif strcmp(F.r,'nred'), disp('nonreduced');
end;
if strcmp(F.p,'prop'), disp('proper');
elseif strcmp(F.p,'nprop'), disp('improper');
end;
   
%end .. @frac/display
