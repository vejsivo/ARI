function [varargout] = old2new(varargin)
%OLD2NEW  Converts the old polynomial format to a POL object
%
% Syntax
%    Q = OLD2NEW(P)
% See also NEW2OLD.

%       Author(s):  S. Pejchova, M. Sebek 31-3-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 16-Apr-1998 11:22:34   $

% Effect on other properties:
% Variables are set to equal to the current value of the
%       global variable symbol.
% UserData are deleted.

ni = nargin;  no = nargout;
if ni<1,
   error('Not enough input arguments.');
end; 
for i=1:ni,
   Ao=varargin{i}; An = Ao;
   if isa(Ao,'double') & ndims(Ao)==2,
      [rA,cA]=size(Ao);
      if (rA>1)&(cA>1)&(~any(Ao(2:rA-1,cA)))&...
         (~any(Ao(rA,1:cA-1)))&isnan(Ao(rA,cA)),
         degA=Ao(1,cA);
         Ao=Ao(1:rA-1,1:cA-1);
         eval(['An=pol(Ao,degA);'],'error(lasterr)');
      end;
   end;
   if no>=i,
      varargout{i}=An;
   else,
      name=inputname(i);
      assignin('caller',name,An);
   end;
end;
if no>ni,
   for i=1:no-ni,
      varargout(ni+i)=[];
   end;
end;

%end .. old2new
