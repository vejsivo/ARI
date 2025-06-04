function C = coef(A,arg2,arg3)
%COEF    Coefficient of fraction
%
% For fraction A, the command  C = COEF(A,k) returns K-th
% coefficient in Laurent series of A (at k-th power of A.v) .
% It has the same effect as  A{k} .
% 
% For fractions in 'z' or in 'z^-1', it is possible 
% to specify the variable for coefficient as an optional
% argument. So, with  A = z/(z-.8),  the commands
%    A{-2},  COEF(A,-2),  COEF(A,2,zi),  COEF(A,-2,z)
% yield 0.64, whereas
%    COEF(A,-2,zi),  COEF(A,2,z)  yield 0.
% With  A = 1/(1-.8*z^-1), the commands
%    A{2},  COEF(A,2),  COEF(A,2,zi),  COEF(A,-2,z)
% yield 0.64, whereas
%    COEF(A,2,z),  COEF(A,-2,zi)  yield 0.
%
% See also: POL/COEF, POL/SUBSREF, RDF/SUBSREF, LDF/SUBSREF,
%           MDF/SUBSREF, SDF/SUBSREF

%       Author:  J. Jezek, 22-Jul-2002
%       Copyright(c) 2002 by Polyx, Ltd.

if nargin<2,
   error('Not enough input arguments.');
end;

isk = 0; var = '';
if isa(arg2,'double'), k = arg2; isk =1;
elseif isa(arg2,'char'), var = arg2;
else
   eval('arg2 = pol(arg2);', ...
      'error(''Invalid 2nd argument.'');');
   [vs1,vs2,vd] = size(arg2);
   if all([vs1,vs2,vd]==1) & all(arg2.c(:,:)==[0,1]),
      var = arg2.v;
   else error('Invalid 2nd argument.');
   end;
end;

if nargin==3,
   if isa(arg3,'double'), k = arg3; isk = 1;
   elseif isa(arg3,'char'), var = arg3;
   else
      eval('arg3 = pol(arg3);', ...
         'error(''Invalid 3rd argument.'');');
      [vs1,vs2,vd] = size(arg3);
      if all([vs1,vs2,vd]==1) & all(arg3.c(:,:)==[0,1]),
         var = arg3.v;
      else error('Invalid 3rd argument.');
      end;
   end;
end;

if ~isa(A,'frac'),
   eval('A = sdf(A);','error(peel(lasterr));');
end;
if ~isk,
   error('Subscript for coefficient is missing.');
end;
if ~isempty(var),
   if strcmp(var,'zi'), var = 'z^-1';
   end;
   if ~isempty(A.v) & ~strcmp(A.v,var),
      if (strcmp(A.v,'z') & strcmp(var,'z^-1')) | ...
            (strcmp(A.v,'z^-1') & strcmp(var,'z')),
         k = -k;
      else
         error('Invalid variable for coefficient.');
      end;
   end;
end;

S.type = '{}';
S.subs = cell(1,1);
S.subs{1} = k;

eval('C = subsref(A,S);','error(peel(lasterr));');

%end .. @frac/coef
