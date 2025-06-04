function C = coef(A,arg2,arg3)
%COEF    Coefficient of two-sided polynomial
%
% The command  C = COEF(A,K)  returns K-th coefficient 
% of two-sided polynomial A. It has the same effect as 
%    C = A{K} .
%
% It is also possible to specify the variable for 
% coefficient. So, with  A = 2z^-1 + 3 + 4z, the commands
%   A{1},  COEF(A,1),  COEF(A,1,z),  COEF(A,-1,zi)
% yield 4, whereas
%   COEF(A,-1), COEF(A,-1,z),  COEF(A,1,zi)
% yield 2.
%
% See also: TSP/SUBSREF

%       Author:  J. Jezek, 22-Jul-2002
%       Copyright(c) 2002 by Polyx, Ltd.
%       $ Revision $  $ Date 28-Feb-2003 $

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

eval('A = tsp(A);','error(peel(lasterr));');
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

eval('A = tsp(A);','error(peel(lasterr));');
if ~isk,
   error('Subscript for coefficient is missing.');
end;
if ~isempty(var),
   if strcmp(var,'zi'), var = 'z^-1';
   end;
   if ~(strcmp(var,'z') | strcmp(var,'z^-1')),
      error('Invalid variable for coefficient.');
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

Ao = A.o;
if ~isempty(Ao), k = k-Ao;
end;
eval('C = coef(A.p,k);','error(peel(lasterr));');

%end .. @tsp/coef
