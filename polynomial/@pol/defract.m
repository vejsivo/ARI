function [Q,C] = defract(F)
%DEFRACT    Defract polynomial
%
% The command  Q = DEFRACT(F)  returns polynomial Q = F.
% The command  [Q,C] = DEFRACT(F)  returns also the class  C = 'pol' .
%
% This macro exists only for completeness.
% See also RDF/DEFRACT, LDF/DEFRACT, MDF/DEFRACT, SDF/DEFRACT.

%       Author:  J. Jezek  27-Apr-2000
%       Copyright(c) 2000 by Polyx, Ltd.

Q = F; C = class(F);

%end .. @pol/defract

