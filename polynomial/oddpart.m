function [Aodd,Aeven] = oddpart(A,string)
%ODDPART  Odd part of a polynomial matrix
%
% Given the polynomial matrix
%   A = A{0} + A{1}*v + A{2}*v^2 + ... + A{d}*v^d,
% the command
%    Aodd = ODDPART(A)  
% returns the odd part of A This contains the terms 
%    Aodd = A{1}*v + A{3}*v^3 + ...
%                  + A{2*floor(d/2)+1}*v^(2*floor(d/2)+1)
% The command
%    Aodd = ODDPART(A,'SQZ')  
% returns the odd part in the  "squeezed" form:
%    Aodd = A{1} + A{3}*v + A{5}*v^2 + ...
%                + A{2*floor(d/2)+1}*v^(floor((d-1)/2))
%
% The second output argument Aeven in both cases returns the 
% even part of the polynomial matrix A.
%
%  See also EVENPART.

%       Author(s):  S. Pejchova 28-07-2000
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 04-Aug-2000 10:32:34  S.Pejchova  $
%                         $Date: 08-Aug-2002  J.Jezek  arg check   $

if ~nargin,
   error('Not enough input arguments.');
end;
if nargin==1,
   eval('[Aodd,Aeven] = oddpart(pol(A));', ...
      'error(peel(lasterr));');
else
   eval('[Aodd,Aeven] = oddpart(pol(A),string);', ...
      'error(peel(lasterr));');
end;

%end .. oddpart
