function N = num(F,R)
%NUM    Numerator of fraction
%
% For fraction F, the command   N = NUM(F)  delivers
% polynomial N, which is numerator of F.
%
% For fraction F and polynomial R, the command
%    NUM(F,R)  sets the numerator of F to R.
%
% See also FRAC/DEN.

%       Author: J.Jezek, 02-Feb-2002
%       Copyright(c) 2002 by  Polyx, Ltd.
%       $ Revision $  $ Date 14-Oct-2002 $

if nargin==1,
   N = F.num;
else
   if isempty(inputname(1)),
      error('1st argument must be a named variable.');
   end;
   eval('N = pol(R); F.num = N;','error(peel(lasterr));');
   assignin('caller',inputname(1),F);
end;

%end .. @frac/num

