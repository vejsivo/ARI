function T = vertcat(varargin)
%VERTCAT      Vertical concatenate two-sided polynomials
%                   T = [ A;B;C ... ]
%
% The command
%    T = [A;B]
% returns the vertical concatenation of the tsp matrices A and B.
% A and B must have the same number of columns. Any number of matrices
% may be concatenated with one pair of brackets. Horizontal and
% vertical concatenation may be combined together as for standard
% MATLAB matrices. The command
%    T = VERTCAT(A,B)
% works similarly.
%
% See also TSP/HORZCAT.

%      Author: J. Jezek  11-8-99
%      Copyright (c) 1999 by Polyx, Ltd.
%      $Revision: 3.0 $  $Date: 29-Sep-1999  13:00:00  $
%                        $Date: 16-Jun-2000  17:00:00  $
%                        $Date: 31-Oct-2000  13:00:00  $
%                        $Date: 25-Jan-2002            $
%                        $Date: 28-Feb-2003            $

% last argument may be TOL
ni = length(varargin);
if ni>=2,
   lastarg = varargin{ni};
   if isa(lastarg,'char'),
      ni = ni-1; tol = str2double(lastarg);
      if ~isempty(lastarg),
         if isnan(tol) | length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
            error('Invalid last argument.');
         end;
      end;     
   end;
end;

h = []; warn = 0;
v = zeros(1,ni); var = cell(1,ni);
for i=1:ni,
   U = varargin{i};
   eval('U = tsp(U);', ...
      'eval(''U = mdf(U);'',''error(''''Invalid argument.'''');'');');
   if isa(U,'mdf'),         % mdf case
      varargin{i} = U;
      eval('T = vertcat(varargin{1:ni});','error(peel(lasterr));');
      eval('T = tsp(T);',';');
      return;
   end;
   Uh = U.h;
   if isempty(h), h = Uh;
   elseif ~isempty(Uh) & ~isnan(Uh) & ~isnan(h) & Uh~=h,
      warn = 1; U.p.h = NaN;
   end;   
   var{i} = U.p; var{i}.v = 'z'; v(i) = U.o;
end;

if warn,
   warning('Inconsistent sampling periods.');
   h = NaN;
end;

oo = min(v);
for i=1:ni,
   var{i} = shift(var{i},v(i)-oo,'z');
end;
Q = 0;
eval('QQ = vertcat(var{1:ni});','error(peel(lasterr));');
T = tsp(QQ); T.h = h;
T = shift(T,oo);

%end .. @tsp/vertcat
