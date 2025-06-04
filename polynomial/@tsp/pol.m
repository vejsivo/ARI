function P = pol(T,name)
%POL   Convert two-sided polynomial to polynomial
%
%    P = POL(T)  or  P = POL(T,'z')  or  P = POL(T,'z^-1')
% This command converts a two-sided polynomial matrix T to polynomial
% in z or in z^-1, if possible. If not possible, error.
%
% See also POL/POL, TSP/DECLASS, TSP/NNEG, TSP/NPOS.

%        Author: J. Jezek  11-10-99
%        Copyright (c) 1999 by Polyx, Ltd.

notconv = 0;
if nargin==1,
   if isempty(T) | tdeg(T)>=0,
      P = nneg(T); return;
   elseif deg(T)<=0,
      P = npos(T); return;
   else notconv =1;
   end;
else
   if isa(name,'pol'),
      [vs1,vs2,vd] = size(name);
      if all([vs1,vs2,vd]==1) & all(name.c(:,:)==[0,1]),
         name = name.v;
      end;
   end;
   if isa(name,'char'),
      if strcmp(name,'z'),
         if isempty(T) | tdeg(T)>=0,
            P = nneg(T); return;
         else notconv =1;
         end;
      elseif strcmp(name,'zi') | strcmp(name,'z^-1'),
         if isempty(T) | deg(T)<=0,
            P = npos(T); return;
         else notconv =1;
         end;
      end;
   end;
end;        
if notconv,
   error('Argument is not convertible to polynomial.');
end;
error('Invalid 2nd argument.');

%end .. @tsp/pol

      
      
         