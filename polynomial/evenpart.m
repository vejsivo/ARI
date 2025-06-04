function [Aeven,Aodd] = evenpart(A,string)
%EVENPART Even part of a polynomial matrix
%
% Given the polynomial matrix
%   A = A{0} + A{1}*v + A{2}*v^2 + ... + A{d}*v^d
% the command
%    Aeven = EVENPART(A)  
% returns the even part of A. This consists of the terms
%    Aeven = A{0} + A{2}*v^2 + A{4}*v^4 + ...
%                 + A{2*floor(d/2)}*v^(2*floor(d/2))
% The command
%    Aeven = EVENPART(A,'SQZ')  
% returns the even part 'squeezed' to the form 
%    Aeven = A{0} + A{2}*v + A{4}*v^2 + ...
%                 + A{2*floor(d/2)}*v^(floor(d/2))
% The second output argument Aodd in both cases returns
% the odd part of A. 
%
% See also ODDPART.

%       Author(s):  S. Pejchova 28-07-2000
%       Copyright (c) 2000 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 04-Aug-2000 10:18:34  S.Pejchova  $
%                         $Date: 08-Aug-2001  J.Jezek  arg check   $

if ~nargin,
   error('Not enough input arguments.');
end;
if nargin==1,
   eval('[Aeven,Aodd] = evenpart(pol(A));', ...
      'error(peel(lasterr));');
else
   eval('[Aeven,Aodd] = evenpart(pol(A),string);', ...
      'error(peel(lasterr));');
end;

%end .. evenpart
