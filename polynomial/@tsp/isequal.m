function c = isequal(varargin);
%ISEQUAL  Test if two-sided polynomials are equal
%
% ISEQUAL(A,B) is 1 if the two arrays are of the same size
% and contain the same values, and 0 otherwise.
%
% ISEQUAL(A,B,C,...) is 1 if all the input arguments are of the
% same size and contain the same values.
%
% ISEQUAL(A,B,'tol',T) or ISEQUAL(A,B,C,...,'tol',T) uses 
% the prescribed tolerance T for comparison. Its default 
% value is the global zeroing tolerance. 

%      Author: J. Jezek, 18-Sep-2001
%      Copyright(c) 2001 by Polyx, Ltd.
%      $ Revision $  $ Date 28-Feb-2003  warning  $

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

c = logical(1); hOK = logical(1);
eval('x = tsp(varargin{1}); [xs1,xs2] = size(x);', ...
   'c = logical(0);');
if c,
   for i = 2:ni,
      d = logical(1);
      eval('y = tsp(varargin{i});', 'd = logical(0);');
      c = c & d;
      if d,
         c = c & all(size(y)==[xs1,xs2]);
         savedw = warning; warning off;
         [th,h,x,y] = testht(x,y);
         hOK = hOK & th;
         eval('c = c & all(all(eq(x,y,tol)));', ...
            'warning(savedw); c = logical(0);');
         warning(savedw);
      end;
   end;
end;
if ~hOK,
   warning('Inconsistent sampling periods.');
end;

%end .. @tsp/isequal
