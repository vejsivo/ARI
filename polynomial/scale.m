function [Q,a] = scale(P,a)
%SCALE     Scale a polynomial matrix
%
% Given a polynomial matrix P of degree n and a scalar a the 
% function call
%    Q = SCALE(P,a)
% results in the polynomial matrix Q given by 
%    Q(s) = (a^n)*P(s/a)
% This leaves the leading coefficient matrix unaltered.
% Without the second input argument, the call
%   [Q,a] = scale(P)
% scales P automatically while taking the scaling factor as
%   a = [norm(Ph)/norm(Pl)]^(1/d)
% Pl is the coefficient matrix of the term with the lowest degree,
% Ph the coefficient matrix of the term with highest degree,
% and d the difference of the highest and lowest degrees.
%
% Automatic scaling makes the norm of the coefficient matrix of 
% the term with lowest power equal to that of the highest power.
%
% See also LINVT, REVERSE.

%      Author: H. Kwakernaak, 1992 and 1998
%      Copyright(c) 1998 by PolyX Ltd

% Preparation

if nargin<1,
   error('Not enough input arguments.');
end;
eval('P = pol(P);', 'error(peel(lasterr));');
if nargin == 2,
   if ~isa(a,'double') | length(a)~=1,
      error('Invalid scale.');
   end;
end;

if isempty(P)
   Q = P; 
   if nargin == 1;
      a = 1; 
   end
   return
elseif deg(P) <= 0
   Q = P; 
   if nargin == 1;
      a = 1; 
   end
   return
end

if nargin == 1
   i = 0;
   while norm(P{i}) == 0
      i = i+1;
   end
   d = deg(P)-i;
   if d > 0
      a = (norm(P{deg(P)})/norm(P{i}))^(1/d);
   else
      a = 1;
   end
end

Q = P;
if a == 1
   return
else
   for j = 0:deg(P)
      Q{j} = P{j}*a^(deg(P)-j);   
   end
end

%end .. scale
