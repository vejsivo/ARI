function Q = smreduce(R)
%SMREDUCE    Small reduce matrix-den fraction
%
% For matrix-den fraction R, the command  Q = SMREDUCE(R)  returns
% matrix-den fraction Q, which is equal to R but modified:
% if the variable symbol is 'z','q','s' or 'p'
%    then all leading coefficients of the denominators are 1
% if the variable symbol is 'z^-1' or 'd'
%    then all trailing coefficients of the denominators are 1.
%
% Tis macro exists only for completeness.
% See also MDF/REDUCE.

%        Author:  J. Jezek  30-Sep-2002
%        Copyright(c) 2002 by Polyx, Ltd
%        $ Revision $  $ Date 14-Oct-2002 $

Q = R;
if strcmp(Q.frac.r,'red'),
   return;
end;

Rc = R.frac.c; Rp = R.frac.p;
Rtp = R.frac.tp; Rtc = R.frac.tc;
if ~isempty(Q),
   if strcmp(Q.frac.v,'z^-1') | strcmp(Q.frac.v,'d'),
      M = tcoef(Q.frac.den,'ent');
   else
      M = lcoef(Q.frac.den,'ent');
   end;
   Q.frac.num = Q.frac.num./M; Q.frac.den = Q.frac.den./M;
end;
props(Q,'red',Rc,Rtc,Rp,Rtp);

%end .. @mdf/smreduce
