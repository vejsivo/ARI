function D = den(F)
%DEN    Denominator of constant
%
% For constant F, the command   D = DEN(F)  delivers
% polynomial D = 1. So, when F is considered as a
% fraction, D is its denominator.
%
% This macro exists only for completeness.
%
% See also NUM, RDF/DEN, lDF/DEN, MDF/DEN, SDF/DEN.

%       Author: J.Jezek, 02-Feb-2002
%       Copyright(c) 2002 by  Polyx, Ltd.
%       $ Revision $  $ Date: 15-Jul-2002 @

eval('N = pol(F);','error(peel(lasterr));');
D = pol(1);

%end .. den

