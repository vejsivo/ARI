function H = norm(F,arg2);
%NORM  Norms of polynomial matrix fractions
%
% For a polynomial matrix fraction F of any kind, the commands
%    NORM(F)
%    NORM(F,Inf)
% compute its H-infinity norm. The commmand
%    NORM(F,2)
% returns the H-2 norm of F.
%
% See also FRAC/HINFNORM, FRAC/H2NORM.

%	Author(s): M. Hromcik, M. Sebek 16-2-98
%	Copyright (c) 1998-2000 by PolyX, Ltd.
%       $Revision: 4.0 $  $Date: 23-Oct-1998 10:28:34   $
%
%		  $ Version 3 - 22-Nov-2000, Martin Hromcik $
%       $             06-Jul-2001, Jan Jezek      $
%       $             25-Jul-2002, Jan Jezek      $

switch nargin,
case 1,
   eval('H = hinfnorm(F);','error(peel(lasterr));');
case 2
   if isempty(arg2) | ...
         (isa(arg2,'double') & length(arg2)==1 & arg2==Inf) | ...
         (isa(arg2,'char') & (strcmp(arg2,'Inf') | strcmp(arg2,'inf'))),
      eval('H = hinfnorm(F);','error(peel(lasterr));');
   elseif (isa(arg2,'double') & length(arg2)==1 & arg2==2) | ...
         (isa(arg2,'char') & strcmp(arg2,'2')),
      eval('H = h2norm(F);','error(peel(lasterr));');
   else
      error('Invalid command option.');
   end;
end;
 
%end .. @frac/norm