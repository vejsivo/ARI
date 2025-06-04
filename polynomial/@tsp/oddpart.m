function [Aodd,Aeven] = oddpart(A,string)
%ODDPART Odd part of a two-sided polynomial matrix
%
% Given the two-sided polynomial matrix
%   A = A{-k}*z^(-k) + ... + A{-2}*z^(-2) + A{-1}*z + 
%       + A{0} + A{1}*z + A{2}*z^2 + ... + A{m}*z^m
% the command
%    Aodd = ODDPART(A)  
% returns the odd part of A. This contains the terms 
%    Aodd = A{-2*fix((k-1)/2)-1}*z^(-2*fix((k-1)/2)-1)+ ... + A{-1}*z^(-1) + 
%           + A{1}*z + A{3}*z^3 + ... + A{2*fix((m-1)/2)+1}*z^(2*fix((m-1)/2)+1)
% The command
%    Aodd = ODDPART(A,'SQZ')  
% returns the odd part in the  "squeezed" form:
%    Aodd = A{-2*fix((k-1)/2)-1}*z^(-fix((k-1)/2)-1)+ ... + A{-2}*z^(-1) + 
%           + A{1} + A{3}*z + A{5}*z^2 + ... + A{2*fix((m-1)/2)+1}*z^(fix((m-1)/2))
%
% The second output argument Aeven returns the even or 'squeezed' even part of A. 
%
%  See also TSP/EVENPART.

%       Author(s):  S. Pejchova 04-08-2000
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 07-Aug-2000 10:03:34  S.Pejchova  $

if ~nargin,
   error('Not enough input arguments.');
elseif ~isa(A,'tsp'),
   error('Invalid 1st argument.');
end;
o_s=A.o;  o_n = o_s;
if nargin==1,
   [Aodd,Aeven] = oddpart(A.p);
elseif isa(string,'char') & strcmp(string,'sqz'),
   [Aodd,Aeven] = oddpart(A.p,string);
   o_n = o_s/2;
else
   error('Invalid command option.');
end;
if mod(o_s,2), 
   Axx = Aeven;  Aeven = Aodd;  Aodd = Axx;
end;
Aodd = shift(Aodd,floor(o_n),'z');
Aeven = shift(Aeven,fix(o_n),'z');

%end .. @tsp/oddpart
