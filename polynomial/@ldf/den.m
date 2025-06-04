function D = den(F,R)
%DEN    Denominator of left-den fraction
%
% For left-den fraction F, the command  D = DEN(F)
% delivers polynomial D, which is denominator of F.
%
% For left-den fraction F and polynomial R, the command
%   DEN(F,R)  sets the denominator of F to R.
%
% See also @LDF/NUM.

%      Author: J.Jezek, 15-Jul-2002
%      Copyright(c) 2002 by Polyx, Ltd.
%      $ Revision $  $ Date 14-Oct-2002 $
%                    $ Date 23-Feb-2003 $

if nargin==1,
   D = F.frac.den;
else
   if isempty(inputname(1)),
      error('1st argument must be a named variable.');
   end;
   eval('F = ldf(F); D = pol(R);','error(peel(lasterr));');
   S.type = '.';
   S.subs = 'den';
   eval('F = subsasgn(F,S,R);','error(peel(lasterr));');
   assignin('caller',inputname(1),F);
end;

%end .. @ldf/den
