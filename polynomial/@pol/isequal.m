function c = isequal(varargin);
%ISEQUAL  Test if polynomials are equal
%
% ISEQUAL(A,B) is 1 if the two arrays are of the same size
% and contain the same values, and 0 otherwise.
%
% A difference between the polynomial variables of A and B causes a 
% warning message, but does not affect the result - for instance,
%   ISEQUAL(s,z)  is 1.    
% However, if one of the variables is 'z' and the other 'z^-1',
% then the variable names play a role in the comparison;
% no warning is issued in such a case. So,
%   ISEQUAL(z,z^-1)  is 0.
%
% ISEQUAL(A,B,C,...) is 1 if all the input arguments are of the
% same size and contain the same values.
%
% ISEQUAL(A,B,'tol',T) or ISEQUAL(A,B,C,...,'tol',T) uses 
% the prescribed tolerance T for comparison. Its default 
% value is the global zeroing tolerance.
%
% See also POL/EQ.

%       Author(s): M. Hromcik, M. Sebek 21-6-01
%       Copyright (c) 1998 by Polyx, Ltd.
%       Modified by J. Jezek, 17-Sep-2001
%                   J. Jezek, 28-Mar-2003  warning

global PGLOBAL;

ni = nargin;
if ni <= 1, error('Not enough input arguments.');
end;

flg = varargin{ni-1};
if isa(flg,'char') & strcmp(flg,'tol'), 
   tol = varargin{ni};
   if ~isempty(tol),
      if ~isa(tol,'double') | length(tol)~=1 | ...
            ~isreal(tol) | tol<0 | tol>1,
         error('Invalid tolerance.');
      end;
   end;
   ni = ni-2;
   if ni <= 1, error('Not enough input arguments.');
   end;
else 
   tol = PGLOBAL.ZEROING;
end;

c = logical(1); vOK = logical(1); hOK = logical(1);
savedw = warning; warning off;
eval('x = pol(varargin{1}); [xs1,xs2] = size(x);', ...
   'warning(savedw); c = logical(0);');
if c,
   for i = 2:ni,
      d = logical(1);
      eval('y = pol(varargin{i});', ...
         'warning(savedw); c = logical(0);');
      c = c & d;
      if d,
         c = c & all(size(y)==[xs1,xs2]);
         [tv,v,x,y] = testvp(x,y);
         vOK = vOK & tv;
         [th,h,x,y] = testhp(x,y,v);
         hOK = hOK & th;
         eval('c = c & all(all(eq(x,y,tol)));', ...
            'warning(savedw); c = logical(0);');
      end;
   end;
end;
warning(savedw);
if ~vOK,
   warning('Inconsistent variables.');
elseif ~hOK,
   warning('Inconsistent sampling periods.');
end;

%end .. @pol/isequal
