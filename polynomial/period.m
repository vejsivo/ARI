function h = period(P,number)
%PERIOD  Display/modify sanpling period of object
%
% The command
%    PERIOD(P,NUMBER)
% sets the sampling period of P equal to the value NUMBER. The command
%    H = PERIOD(P)
% returns the current value of the sampling period of P.

%       Author:  J. Jezek  03-Jul-2002
%       Copyright(c) 2002 by Polyx, Ltd.

ni = nargin;
if ~ni, error('Not enough input arguments.');
end;

I1 = strmatch(class(P),{'pol';'tsp';'frac';'rdf';'ldf';'mdf';'sdf';'double'},'exact');
if isempty(I1),
   error('1st agument is not a polynomial object.');
end;

if ni==1,
   if I1==8, h = [];
   else h = P.h;
   end;
else,
   if isempty(inputname(1)),
      error('1st argument must be a named variable.');
   end;
   if I1==8,
      if ~isa(number,'double') | (~isempty(number) & ...
            length(number)~=1 | ~isreal(number) | number<0 | isinf(number)),
         error('Invalid sampling period.');
      end;
      h = [];
   else
      eval('P.h = number;','error(peel(lasterr));');
      % Assign P in caller's workspace
      assignin('caller',inputname(1),P);
      h = P.h;
   end;
end;

%end .. period


