function [B,U] = otrowred(A);
%OTRORED    Orthogonal trailing row reduction of
%           a polynomial matrix
%
% For square nonsingular polynomial matrix A in vaziable
% 'z', the command
%       [B,U] = OTROWRED(A)
% brings A(z) to B(z), where  B = U*A  is such that
% the absolute term B{0} is nonsingular upper triahgular.
% The resulting U is an orthogonal (in the complex case:
% unitary) polynomial matrix in  variable 'z^-1',
% i.e. it satisfies   CTRANSPOSE(U) * U = I .
%
% See also ROWRED, COLRED, SPF.

%       Author:  J. Jezek  13-Mar-2003
%       Copyright(c) 2003 by Polyx, Ltd.

if nargin<1,
   error('Not enough input arguments.');
end;
eval('A = pol(A);','error(peel(lasterr));');
Av = A.v;
if ~(isempty(Av) | strcmp(Av,'z')),
   error('Invalid variable symbol: must be ''z''.');
end;
[rA,cA] = size(A);
if rA~=cA,
   error('Matrix is not square.');
end;
if issingular(A),
   error('Matrix is singular.');
end;

B = A; rB = rA; cB = cA;
U = eye(rB);
if rB==0, return;
end;

no = nargout;
while(1),
   L = B{0};
   [Q,L] = qr(L); Qt = Q';
   B = Qt*B;
   if no==2, U = Qt*U;
   end;
   r = rank(L);
   if r==rB, return;
   end;
   Blow = B(r+1:rB,:);
   t = tdeg(Blow,'row');
   if any(t==Inf),
      error('Matrix is singular.');
   end;
   D = diag([ones(r,1);zi.^t]);
   if no==2, U = D*U;
   end;
   t = repmat(t,1,cB);
   B(r+1:rB,:) = shift(Blow,-t);
end;

%end .. @pol/otrowred


