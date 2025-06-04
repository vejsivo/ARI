function C = coef(A,arg2,arg3)
%COEF    Coefficient of polynomial
%
% The command  C = COEF(A,K)  returns K-th coefficient
% of polynomial A. It has the same effect as  C = A{K} .
%
% For polynomials in 'z' or in 'z^-1', it is possible 
% to specify the variable for coefficient as an optional
% argument. So, with A = 3+4z, the commands
%    A{1},  COEF(A,1),  COEF(A,1,z),  COEF(A,-1,zi)  yield 4,
% whereas  COEF(A,1,zi),  COEF(A,-1,z)  yield 0.
% With A = 3+4z^-1, the commands
%    A{1},  COEF(A,1),  COEF(A,1,zi),  COEF(A,-1,z)  yield 4,
% whereas  COEF(A,1,z),  COEF(A,-1,zi)  yield 0.
%
% See also: POL/SUBSREF

%      Author:  J. Jezek, 19-Apr-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 25-May-2000 $
%                    $ Date 24-Jun-2001 $
%                    $ Date 22-Jul-2002 $

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

eval('A = pol(A);','error(peel(lasterr));');
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

%end .. @pol/coef
