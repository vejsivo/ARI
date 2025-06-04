function [varargout] = new2old(varargin)
%NEW2OLD  Converts a POL object into the old format
%used in the Polynomial Toolbox versions 1.5 and 1.6
% 
% Syntax
%    P = NEW2OLD(Q)
% See also OLD2NEW.

%       Author(s):  S. Pejchova, M. Sebek 26-3-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 16-Apr-1998 11:51:34   $

% Effect on other properties:
% Outputs are standard Matlab matrices.

ni =nargin; no = nargout;
if ni<1,
   error('Not enough input arguments.');
end;
for i=1:ni,
   An=varargin{i};
   if isa(An,'pol'),
      Ad = An.d; As = An.s; [As1,As2]=size(An);  Ac = An.c;
      if isempty(Ad),
         Ao = [];
      elseif isinf(Ad),
         Ao = zeros(As1+1,As2+1);
         Ao(1,As2+1)=-Inf;
         Ao(As1+1,As2+1)=NaN;
      else,
         Ao = zeros(As1+1,1+As2*(Ad+1));
         Ao(1:As1,1:As2*(Ad+1)) = Ac(:,:);
         Ao(1,1+As2*(Ad+1)) = Ad;
         Ao(As1+1,1+As2*(Ad+1)) = NaN;
      end;
   else,
      Ao = An;
   end;
   if no>=i,
      varargout{i}=Ao;
   else,
      name=inputname(i);
      assignin('caller',name,Ao);
   end;
end;
if no>ni,
   for  i=1:no-ni,
      varargout{ni+i}=[];
   end;
end;

%end .. new2old
