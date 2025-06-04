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
% See also POL/EVENPART.

%       Author(s):  S. Pejchova 16-10-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 16-Oct-1998 15:30:34   $
%       $Revision: 3.0 $  $Date: 04-Aug-2000 10:36:34  S.Pejchova  $
%                         $Date: 17-Jul-2001 J.Jezek sampling period $

ni=nargin; k=0;
if ~ni,
   error('Not enough input arguments.');
elseif ~isa(A,'pol'),
   error('Invalid 1st argument.');
elseif ni==2,
   if ischar(string)&strcmp(string,'sqz'), k=1;
   else, error('Invalid command option.');
   end;
end

Ad=A.d; Ac=A.c; [As1,As2]=size(A); var=A.v; Ah=A.h;

if isempty(Ad) | isinf(Ad), 
   Aodd=A; Aeven=A;  return;
elseif Ad==0,
   Aodd=pol(zeros(As1,As2)); Aeven=A;  return;
end;

if ~k,
   Aodd=zeros(As1,As2,Ad+1); Aeven=Aodd;
   Aodd(:,:,2:2:end)=Ac(:,:,2:2:end);
   Aodd=pol(Aodd(:,:),Ad,var,Ah);   
   if nargout==2,
      Aeven(:,:,1:2:end)=Ac(:,:,1:2:end);
      Aeven=pol(Aeven(:,:),Ad,var,Ah);
   end;
else,
   Aodd=Ac(:,:,2:2:end);
   Aodd=pol(Aodd(:,:),floor((Ad-1)/2),var,Ah*2);
   if nargout==2,
      Aeven=Ac(:,:,1:2:end);
      Aeven=pol(Aeven(:,:),floor(Ad/2),var,Ah*2);
   end;
end;

%end .. @pol/oddpart
