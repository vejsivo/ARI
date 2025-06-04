function H = hurwitz(A,k)
%HURWITZ  Create the Hurwitz matrix corresponding to a polynomial matrix
%
% The command
%    H = HURWITZ(A,K)  
% creates the constant Hurwitz matrix H of order K (K-by-K blocks)  
% corresponding to the polynomial matrix A. If
%     A = A{0} + A{1}*v + A{2}*v^2 + ... + A{d}*v^d
% then
%     H = [ A{d-1}  A{d-3}  A{d-5}  ...    ...   |  block row 1.
%         | A{d}    A{d-2}  A{d-4}  ...    ...   |  block row 2.
%         | 0       A{d-1}  A{d-3}  A{d-5} ...   |  ...
%         | 0       A{d}    A{d-2}  A{d-4} ...   |  ...
%         | 0       0       ...     ...    A{0}  ]  block row K.
%
% The default value of K is D.

%       Author(s):  S. Pejchova, M. Sebek 25-9-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 02-Nov-1998 14:41:34   $
%       $Revision: 3.0 $  $Date: 08-Aug-2000  S. Pejchova   $

% Effect on other properties:
% H is a standard Matlab matrix.

ni=nargin;
narginchk(1,2);
% error(nargchk(1,2,ni));	%REMOVED IN NEW MATLABS

eval('A = pol(A);','error(peel(lasterr));');
Ad=A.d; Ac=A.c; [As1,As2]=size(A);  H=[];
if ni==1,  k=Ad;  end;
if ~isa(k,'double') | any(size(k)~=1) | ~(isfinite(k)) | (k~=fix(abs(k))),
   error('The requested order must be a nonnegative integer.');
end;
if isempty(Ad) | isinf(Ad) | ~k, return; end;
if k>Ad, Ac(:,:,k+1)=zeros(As1,As2);
elseif k<Ad, Ac=Ac(:,:,Ad-k+1:end);
end;
if k==1, H=Ac(:,:,1); return; end;
P_E=flipdim(Ac(:,:,1:2:end),3);
P_O=flipdim(Ac(:,:,2:2:end),3);
H = zeros(k*As1,k*As2);

if mod(k,2),
   Ha1=P_E(:,:);                              %K is odd
   Hab=[Ha1; P_O(:,:)];
   if k>2,
      H((k-1)*As1+1:k*As1,floor(k/2)*As2+1:k*As2)=Ha1;
   end;
else,
   Hab=[P_O(:,:),zeros(As1,As2); P_E(:,:)];   %K is even
end;
l_h=size(Hab,2);
H(1:2*As1,1:l_h)=Hab;
for ii=1:floor(k/2)-1,
   H(2*(ii)*As1+1:2*(ii+1)*As1,As2*(ii)+1:As2*(ii)+l_h)=Hab;
end;

%end .. hurwitz
