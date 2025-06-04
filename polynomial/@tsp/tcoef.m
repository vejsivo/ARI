function [T,D] = tcoef(A,varargin)
%TCOEF  Trailing coefficients and trailing degrees of a tsp matrix A
%
% TCOEF(A)       default, the same as TCOEF(A,'mat')
% TCOEF(A,'mat') returns the trailing coefficients matrix of A.
% TCOEF(A,'ent') returns the trailing scalar coefficients of entries of A.
% TCOEF(A,'row') returns the row trailing coefficient matrix of A.
% TCOEF(A,'col') returns the column trailing coefficient matrix of A.
%
% By an optional input argument VAR, it is possible to specify
% that the trailing coefficient is to be understood by the highest power
% of 'z' or 'z^-1'.
%
% The second output argument in all cases returns the
% corresponding matrix or vector of trailing degrees, the same as the first
% argument in function TDEG.
%
% See also TSP/TDEG, TSP/LCOEF, TSP/DEG.

%       Author:  J. Jezek  11-8-99
%       Copyright (c) by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 29-Sep-1999  13:00:00  $
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

if recip,
   if nargout<2,
      eval('T = lcoef(A,string);','error(peel(lasterr));');
   else
      eval('[T,D] = lcoef(A,string);','error(peel(lasterr));');
   end;
   D = -D;
   return;
end;      

Am = mirror(A);
eval('[T,D] = lcoef(Am.p,string);', 'error(peel(lasterr));');
if ~isempty(Am.o), D = -(D + Am.o); end;

%end .. @tsp/tcoef


