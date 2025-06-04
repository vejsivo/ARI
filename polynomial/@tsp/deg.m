function [D,C] = deg(A,varargin)
%DEG  Degree of two-sided polynomial
%
% DEG(A)       default, the same as DEG(A,'mat').
% DEG(A,'mat') returns the degrees of matrix A.
% DEG(A,'ent') returns the matrices of degrees of A entries.
% DEG(A,'row') returns the column vectors of row degrees of A.
% DEG(A,'col') returns the row vectors of column degrees of A.
%
% By an optional input argument VAR, it is possible to specify
% that the degree is to be understood by the highest power
% of 'z' or 'z^-1'.
%
% The second output argument in all cases returns the
% corresponding leading matrix of cofficients, the same as the
% first and argument in function LCOEF.
%
% See also TSP/LCOEF, TSP/TDEG, TSP/TCOEF.

%       Author:  J. Jezek  11-8-99
%       Copyright (c) by Polyx, Ltd.
%       $Revision: 3.0 $  $ Date 29-Sep-1999  13:00:00  $
%                         $ Date 17-Nov-2002  z,z^-1    $

if ~isa(A,'tsp'),
   error('Some argument but not 1st is invalidly tsp.');
end;
recip = logical(0);

string = 'mat';
li = length(varargin);
if li>0,
   for i = 1:li,
      arg = varargin{i};
      if isa(arg,'pol'),
         [vs1,vs2,vd] = size(arg);
         if all([vs1,vs2,vd]==1) & all(arg.c(:,:)==[0 1]),
            arg = arg.v;
         else
            error(['Invalid ',nth(i+1),' argument.']);
         end;
      end;
      if ischar(arg),
         if strcmp(arg,'zi'), arg = 'z^-1';
         end;
         I = strmatch(arg,{'s';'p';'z^-1';'d';'z';'q'},'exact');
         if ~isempty(I),
            if strcmp(arg,'z'),
            elseif strcmp(arg,'z^-1'),
               recip = logical(1);
            else
               error('Invalid variable symbol.');
            end; 
         else
            string = arg;
         end;
      else
         error(['Invalid ',nth(i+1),' argument.']);
      end;
   end;
end;        
  
if recip,
   if nargout<2,
      eval('D = tdeg(A,string);','error(peel(lasterr));');
   else
      eval('[D,C] = tdeg(A,string);','error(peel(lasterr));');
   end;
   D = -D;
   return;
end;

eval('[D,C] = deg(A.p,string);', 'error(peel(lasterr));');
if ~isempty(A.o), D = D + A.o;
end;

%end .. @tsp/deg



