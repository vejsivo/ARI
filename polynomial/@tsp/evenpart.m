function [Aeven,Aodd] = evenpart(A,string)
%EVENPART Even part of a two-sided polynomial matrix
%
% Given the two-sided polynomial matrix
%   A = A{-k}*z^(-k) + ... + A{-2}*z^(-2) + A{-1}*z + 
%       + A{0} + A{1}*z + A{2}*z^2 + ... + A{m}*z^m
% the command
%    Aeven = EVENPART(A)  
% returns the even part of A. This consists of the terms
%    Aeven = A{-2*fix(k/2)}*z^(-2*fix(k/2)) + ... + A{-2}*z^(-2) + 
%            + A{0} + A{2}*z^2 + A{4}*z^4 + ... + A{2*fix(m/2)}*z^(2*fix(m/2))
% The command
%    Aeven = EVENPART(A,'sqz')  
% returns the even part 'squeezed' to the form 
%    Aeven = A{-2*fix(k/2)}*z^(-fix(k/2))+ ...+ A{-2}*z^(-1) + 
%            + A{0} + A{2}*z + A{4}*z^2 + ... + A{2*fix(m/2)}*z^(fix(m/2))
%
% The second output argument Aodd  returns the odd or 'squeezed' odd part of A. 
%
% See also TSP/ODDPART.

%       Author(s):  S. Pejchova 04-08-2000
%       Copyright (c) 2000 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 04-Aug-2000 14:57:34  S.Pejchova  $

if ~nargin,
   error('Not enough input arguments.');
elseif ~isa(A,'tsp'),
   error('Invalid 1st argument.');
end;
o_s = A.o; o_n = o_s;
if nargin==1,
   [Aeven,Aodd] = evenpart(A.p);
elseif isa(string,'char') & strcmp(string,'sqz'),
   [Aeven,Aodd] = evenpart(A.p,string);
   o_n = o_s/2;
else
   error('Invalid command option.');
end;
if mod(o_s,2), 
   Axx = Aeven;  Aeven = Aodd;  Aodd = Axx;
end;
Aeven = shift(Aeven,fix(o_n),'z');
Aodd = shift(Aodd,floor(o_n),'z');

% end ../evenpart.m
