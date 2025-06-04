function [D,C] = tdeg(A,varargin)
%TDEG  Trailing degree of two-sided polynomial
%
% TDEG(A)       default, the same as TDEG(A,'mat').
% TDEG(A,'mat') returns the trailing degrees of matrix A.
% TDEG(A,'ent') returns the matrices of trailing degrees of A entries.
% TDEG(A,'row') returns the column vectors of trailing row degrees of A.
% TDEG(A,'col') returns the row vectors of trailing column degrees of A.
%
% By an optional input argument VAR, it is possible to specify
% that the trailing degree is to be understood by the lowest power
% of 'z' or 'z^-1'.

% The second output argument in all cases returns the
% corresponding trailing matrix of cofficients, the same as the
% first and argument in function TCOEF.
%
% See also TSP/TCOEF, TSP/DEG, TSP/LCOEF.

%       Author:  J. Jezek  11-8-99
%       Copyright (c) by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 29-Sep-1999  13:00:00  $
%                         $Date: 13-Oct-1999  12:00:00  $
%                         $Date: 17-Nov-2002  z,z^-1    $

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
   
no = nargout;
if recip,
   if no<2,
      eval('D = deg(A,string);','error(peel(lasterr));');
   else
      eval('[D,C] = deg(A,string);','error(peel(lasterr));');
   end;
   D = -D;
   return;
end;
      
Am = mirror(A);
eval('[D,C] = deg(Am.p,string);', 'error(peel(lasterr));');
if ~isempty(Am.o), D = -(D + Am.o);
end;

%end .. @tsp/tdeg


