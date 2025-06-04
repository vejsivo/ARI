function H = h2norm(F,arg2);
%H2NORM  H2-norm of a stable polynomial matrix fraction
%
% The command
%    H2NORM(F)
% where F is a polynomial fraction of any kind computes its H2-norm.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
% 
% See also FRAC/HINFNORM, FRAC/NORM.

%       Author(s): M. Hromcik, M. Sebek 22-10-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 27-Oct-1998 10:28:34   $
%       $Revision: 3.0 $  $Date: 23-Sep-1999 $
%       $Revision: 3.1 $  $Date: 23-Nov-1999 $
%
%		  $ Version 3 - 22-Nov-2000, Martin Hromcik $
%       $             25-Jul-2002, Jan Jezek      $
%       $             14-Oct-2002, Jan Jezek      $           

global PGLOBAL;
tol = PGLOBAL.ZEROING;

switch nargin,
  case 1,
  case 2,
   if ~isempty(arg2),
      if isnumeric(arg2) & all(size(arg2)==1) & ...
            isreal(arg2) & arg2>=0 & arg2<=1,
         tol = arg2;
      else error('Invalid 2nd argument.');
      end;
   end;
end;  	% switch ...

Fcl = class(F);
if strcmp(Fcl,'frac'),
   error('Invalid 1st argument.');
end;

if isa(F,'rdf'),
   N = pol(F.num); D = pol(F.den);
   eval('H = h2norm(N,D,''r'',tol);', ...
      'error(peel(lasterr));');
elseif isa(F,'ldf'),
   N = pol(F.num); D = pol(F.den);
   eval('H = h2norm(N,D,''l'',tol);', ...
      'error(peel(lasterr));');
else
   [s1,s2] = size(F);
   if s1>s2,   
      F = rdf(F);
      N = pol(F.num); D = pol(F.den);
      eval('H = h2norm(N,D,''r'',tol);', ...
         'error(peel(lasterr));');
   else
      F = ldf(F);
      N = pol(F.num); D = pol(F.den);
      eval('H = h2norm(N,D,''l'',tol);', ...
         'error(peel(lasterr));');
   end;
   
end;
  
%end .. @frac/h2norm
