function r = eig(A, arg2, arg3, arg4);
%EIG   Find the poles of a fraction.
%
% P = EIG(F) where F is a fraction object returns 
%            the vector of (finite) poles of F. 
%
% For possible further arguments, see POL/ROOTS.

%       Author(s): M. Hromcik
%       Copyright (c) 2000 by Polyx, Ltd.
%       $Date: 12-Oct-2000            Version 3, M. Hromcik 
%       $ Revision $  $ Date 17-Jul-2001  J.Jezek  arg checking  $
%                     $ Date 25-Jul-2002  J.Jezek  class 'frac'  $
%                     $ Date 14-Oct-2002  J.Jezek  $

global PGLOBAL;

all = 0;
met = 'det';
tol = PGLOBAL.ZEROING;

ni = nargin;
if ni>=2,
   if ~isempty(arg2);
      if ischar(arg2),
         if strcmp(arg2,'all'),
            all = 1;
         else met = arg2;
         end;
      elseif isnumeric(arg2),
         tol = arg2;
      else
         error('Invalid 2nd argument.');
      end;
   end;
end;

if ni>=3,
   if ~isempty(arg3);
      if ischar(arg3),
         if strcmp(arg3,'all'),
            all = 1;
         else met = arg3;
         end;
      elseif isnumeric(arg3),
         tol = arg3;
      else
         error('Invalid 3rd argument.');
      end;
   end;
end;
      
if ni>=4,
   if ~isempty(arg4);
      if ischar(arg4),
         if strcmp(arg4,'all'),
            all = 1;
         else met = arg4;
         end;
      elseif isnumeric(arg4),
         tol = arg4;
      else
         error('Invalid 4th argument.');
      end;
   end;
end;

As = A.s;
Acl = class(A);
if strcmp(Acl,'frac'),
   error('Invalid 1st argument.');
end;

if isa(A, 'mdf'), 
   r = [];
   for i=1:As(1),
      for j=1:As(2),
         eval('rr = roots(A.den(i,j),all,met,tol);', ...
            'error(peel(lasterr));');
         r = [r; rr];
      end;
   end;
else		% sdf, ldf or rdf
   eval('r = roots(A.den,all,met,tol);', ...
      'error(peel(lasterr));');
end;

%end .. @frac/eig
