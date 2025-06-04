function D = den(F)
%DEN    Denominator of polynomial
%
% For polynomial F, the command   D = DEN(F)  delivers
% polynomial D = 1. So, when F is considered as a
% fraction, D is its denominator.
%
% This macro exists only for completeness.
%
% See also POL/NUM, RDF/DEN, LDF/DEN, MDF/DEN, SDF/DEN. 

%       Author: J.Jezek, 02-Feb-2002
%       Copyright(c) 2002 by  Polyx, Ltd.
%       $ Revision $  $ Date: 15-Jul-2002 $

D = pol(1);

%end .. @pol/den
