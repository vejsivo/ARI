function F = sdf(T,var)
%SDF    Convert two-sided polynomial to scalar-den fraction
%
% The command  F = SDF(T)  converts two-sided polynomial T to
% scalar-denominator fraction in variable 'z' or 'z^-1',
% decided standardly by PGLOBAL.DISCRVAR or by optional input
% argument VAR, e.g.  F = SDF(T,'z^-1') .

%       Author:  J. Jezek  10-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision 3.0 $  $ Date 28-Apr-2000 $
%                         $ Date 24-May-2000 $

global PGLOBAL;

if nargin==1 | isempty(var),
   var = PGLOBAL.DISCRVAR;
elseif isa(var,'pol'),
   [vs1,vs2,vd] = size(var);
   if all([vs1,vs2,vd]==1) & all(var.c(:,:)==[0,1]),
      var = var.v;
   else
      error('Invalid variable symbol.');
   end;
end;
if ~isa(var,'char') | (~strcmp(var,'z') & ...
      ~strcmp(var,'z^-1') & ~strcmp(var,'zi')),
   error('Invalid variable symbol.');
end;
if strcmp(var,'zi'),
   var = 'z^-1';
end;

if T.o>=0,
   F = sdf(shift(T.p,T.o,'z'));
else
   zz = z; zz.h = T.h;
   F = sdf(T.p,zz^abs(T.o));
end;

if strcmp(var,'z^-1'),
   F = reverse(F);
end;

%end .. @tsp/sdf
