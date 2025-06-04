function P=mono(n);
%MONO  Create the N-th power of the current global variable
%
% If the current global variable string is s then the command
%    P = MONO(N)
% creates the polynomial (monomial) P = s^N. 
%
% If N is a vector or matrix then the vector or matrix of
% monomials is returned.
%
% See also S, P, Z, Q, D, ZI, V.

%       Author(s):  S. Pejchova, M. Sebek 16-4-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 17-Apr-1998 10:50:34   $
%       $Revision: 3.0 $  $Date: 18-Jul-2000 11:50:34  S.Pejchova  $
%       $Revision: 3.1 $  $Date: 9-Oct-2002 15:11:55  Z. Hurak  $

% Effect on other properties:
% P.v: Variable is set equal to the current value of global
%       variable symbol.
% P.u: UserData are deleted.

if nargin ~=1,
   error('Not enough input arguments.');
elseif ndims(n)>2|~isa(n,'double'),
  error('Argument is not two-dimensional double.');
end;
if isempty(n), P=pol(zeros(size(n))); return; 
end;
n(n==-1)=-2;
n(isinf(n)&(n<0))=-1;
n=n+1;
if any(n(:)<0),
   error('Negative degree > -Inf');
elseif any(~isfinite(n(:))) | (any(round(n(:))-n(:))),   
   error('Degree is not integer.');
end;
degp=max(max(n));
if degp==0,
   P=pol(zeros(size(n)));
else,
   [s1,s2]=size(n);
   if s1==1&s2==1,
      Pc=[zeros(1,degp-1), 1];
   else,
      Pc=repmat(reshape(1:degp,[1 1 degp]),[s1 s2]);
      nc=repmat(n,[1 1 degp]);
      Pc=+(Pc(:,:)==nc(:,:));
   end;
      P=pol(Pc,degp-1);
end;

%end .. mono
