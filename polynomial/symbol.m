function Sy = symbol(P,string)
%SYMBOL  Display/modify variable symbol of object
%
% The commmand
%    SYMBOL(P,STRING) 
% sets the variable symbol of P equal to the value STRING. The commmand
%   SY = SYMBOL(P)     
% returns the current value of the variable symbol of P.

%       Author(s):  S. Pejchova, M. Sebek 25-9-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 17-Mar-2000 16:28:50   $
%                         $Date: 08-Aug-2001  J.Jezek   $
%                         $Date: 04-Jul-2002  J.Jezek   $

ni = nargin;
if ~ni, error('Not enough input arguments.');
end;

I1=strmatch(class(P),{'pol';'tsp';'frac';'rdf';'ldf';'mdf';'sdf';'double'},'exact'); 
if isempty(I1), 
   error('1st argument is not a polynomial object.'); 
end;

if ni==1,
   if I1==8, Sy = '';
   else Sy = P.v;
   end;
else,
   if isempty(inputname(1)),
      error('1st argument must be a named variable.')
   end;
   if I1==8,
      test1 = 0;
      if isa(string,'pol'),
         [vs1,vs2,vd]=size(string);
         if all([vs1,vs2,vd]==1)&(~any(string.c(:,:)-[0,1])),
            string = string.v;
         end; 
      end;
      if isstr(string),
         if isempty(string), test1 = 1;
         else
            if strcmp(string,'zi'), string='z^-1';
            end;
            I2 = strmatch(string,{'s';'p';'z^-1';'d';'z';'q'},'exact');
            if ~isempty(I2), test1 = 1;
            end;
         end;
      end;
      if ~test1, error('Invalid variable symbol.');
      end;
      Sy = '';
      
   else
      eval('P.v = string;','error(peel(lasterr));');   
      % Assign P in caller's workspace
      assignin('caller',inputname(1),P);
      Sy = P.v;
   end;
end;

%end .. symbol
