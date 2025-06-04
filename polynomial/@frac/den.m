function D = den(F,R)
%DEN    Denominator of fraction
%
% For fraction F, the command   D = DEN(F)  delivers
% polynomial D, which is denominator of F.
%
% For fraction F and polynomial R, the command
%   DEN(F,R)  sets the denominator of F to R.
%
% See also FRAC/NUM.

%       Author: J.Jezek, 02-Feb-2002
%       Copyright(c) 2002 by  Polyx, Ltd.
%       $ Revision $  $ Date: 15-Jul-2002 $
%                     $ Date: 14-Oct-2002 $

if nargin==1,
   D = F.den;
else
   if isempty(inputname(1)),
      error('1st argument must be a named variable.');
   end;
   eval('D = pol(R); F.den = D;','error(peel(lasterr));');
   assignin('caller',inputname(1),F);
end;

%end .. @frac/den
