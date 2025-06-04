function [Aeven,Aodd] = evenpart(A,string)
%EVENPART  Even part of a polynomial matrix
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
% See also POL/ODDPART.

%       Author(s):  S. Pejchova 16-10-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 16-Oct-1998 15:23:34   $
%       $Revision: 3.0 $  $Date: 28-Jul-2000 10:13:34  S.Pejchova  $

ni=nargin; k=0;
if ~ni,
   error('Not enough input arguments.');
elseif ~isa(A,'pol'),
   error('Invalid 1st argument.');
elseif ni==2,
   if ischar(string)&strcmp(lower(string),'sqz'), k=1;
   else, error('Invalid command option.');
   end;
end

Ad=A.d; Ac=A.c; [As1,As2]=size(A); var=A.v; Ah=A.h;

if isempty(Ad) | isinf(Ad), 
   Aeven=A; Aodd=A;  return;
elseif Ad==0,
   Aeven=A; Aodd=pol(zeros(As1,As2));  return;
end;

if ~k,
   Aeven=zeros(As1,As2,Ad+1); Aodd=Aeven;
   Aeven(:,:,1:2:end)=Ac(:,:,1:2:end);
   Aeven=pol(Aeven(:,:),Ad,var,Ah);
   if nargout==2,
      Aodd(:,:,2:2:end)=Ac(:,:,2:2:end);
      Aodd=pol(Aodd(:,:),Ad,var,Ah);
   end;
else,
   Aeven=Ac(:,:,1:2:end);
   Aeven=pol(Aeven(:,:),floor(Ad/2),var,Ah*2);
   if nargout==2,
      Aodd=Ac(:,:,2:2:end);
      Aodd=pol(Aodd(:,:),floor((Ad-1)/2),var,Ah*2);
   end;
end;

%end .. @pol/evenpart
