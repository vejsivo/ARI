function U = userdata(P,argm)
%USERDATA  Display/modify user data of object
%
% The commmand
%    USERDATA(P,ARGM)  
% sets the user data of P equal to the value ARGM. The commmand
%    U = USERDATA(P)  
% returns the current value of the user data of P

%       Author(s):  S. Pejchova, M. Sebek 25-9-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 17-Mar-2000 15:23:50   $
%                         $Date: 08-Aug-2001  J.Jezek   $

ni = nargin;
if ~ni, error('Not enough input arguments.'); end;
I=strmatch(class(P),{'pol';'tsp';'frac';'rdf';'ldf';'mdf';'sdf'},'exact');
if isempty(I), 
   error('1st argument is not a polynomial object.'); 
end;

if ni==1,
   U=P.u;
else,
  if isempty(inputname(1)),
      error('1st argument must be a named variable.')
  end;
  P.u = argm;    
  % Assign P in caller's workspace
  assignin('caller',inputname(1),P);
  U=P.u;
end;

%end .. userdata
