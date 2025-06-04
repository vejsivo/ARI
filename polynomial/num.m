function N = num(F,R)
%NUM    Numerator of constant
%
% For constant F, the command   N = NUM(F)  delivers
% polynomial N = F. So, when F is considered as
% a fraction, N is its numerator.
%
% For constant F and polynomial R, the command
%    NUM(F,R)  sets  F = R . So, when F is considered
% as a fraction, R is now its numerator.
%
% This macro exists only for completeness.
%
% See also DEN, RDF/NUM, LDF/NUM, MDF/NUM, SDF/NUM.

%       Author: J.Jezek, 02-Feb-2002
%       Copyright(c) 2002 by  Polyx, Ltd.
%       $ Revision $  $ Date: 15-Jul-2002 $

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
eval('N = pol(F);','error(peel(lasterr));');

if ni==2,
   if isempty(inputname(1)),
      error('1st argument must be a named variable.');
   end;
   eval('N = pol(R);','error(peel(lasterr));');
   F = declass(N);
   assignin('caller',inputname(1),F);
end;

%end .. num
