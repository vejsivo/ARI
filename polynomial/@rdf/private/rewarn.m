function rewarn(a)
%REWARN     Restore warning
%
% The command  REWARN(A)  restores the warning
% status from variable A.

%     Author:  J.Jezek  30-Nov-2000
%     Copyrigh(c) 2000 by Polyx,Ltd.

if nargin>=1 & isa(a,'char'),
   if strcmp(a,'on'), warning on;
   elseif strcmp(a,'off'), warning off;
   elseif strcmp(a,'backtrace'), warning backtrace;
   elseif strcmp(a,'debug'), warning debug;
   end;
end;

%end .. @rdf/private/rewarn

