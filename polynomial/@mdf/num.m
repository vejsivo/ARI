function N = num(F,R)
%NUM    Numerator of matrix-den fraction
%
% For matrix-den fraction F, the command  N = NUM(F)
% delivers polynomial N, which is numerator of F.
%
% For matrix-den fraction F and polynomial R, the command
%   NUM(F,R)  sets the numerator of F to R.
%
% See also @MDF/DEN.

%      Author: J.Jezek, 15-Jul-2002
%      Copyright(c) 2002 by Polyx, Ltd.
%      $ Revision $  $ Date 14-Oct-2002 $
%                    $ Date 28-Feb-2003 $

if nargin==2,
   if isempty(inputname(1)),
      error('1st argument must be a named variable.');
   end;
   eval('F = mdf(F); N = pol(R);','error(peel(lasterr));');
   S.type = '.';
   S.subs = 'num';
   eval('F = subsasgn(F,S,N);','error(peel(lasterr));');
   assignin('caller',inputname(1),F);
end;
N = F.frac.num;

%end . @mdf/num