function N = num(F,R)
%NUM    Numerator of polynomial
%
% For polynomial F, the command  N = NUM(F)  delivers
% polynomial N = F. So, when F is considered as
% a fraction, N is its numerator.
%
% For polynomial F and polynomial R, the command
%     NUM(F,R)  sets  F = R. So, when F is considered
% as a fraction, R is now its numerator.
%
% This macro exists only for completeness.
%
% See also POL/DEN, RDF/NUM, LDF/NUM, MDF/NUM, SDF/NUM.

%       Author: J.Jezek, 02-Feb-2002
%       Copyright(c) 2002 by  Polyx, Ltd.
%       $ Revision $  $ Date: 15-Jul-2002 $

eval('N = pol(F);','error(peel(lasterr));');

if nargin==2,
   if isempty(inputname(1)),
      error('1st argument must be a named variable.');
   end;
   eval('F = pol(R);','error(peel(lasterr));');
   assignin('caller',inputname(1),F);
   N = F;
end;

%end .. @pol/num
