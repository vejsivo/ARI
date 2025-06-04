function [B,H,F,bsizes] = bhf(M,tol)
%BHF   Convert polynomial into state space in block Hessenberg form
%
% For polynomial M, convertible to a strictly proper fraction,
% the function
%    [F,G,H,BSIZES] = BHF(M,TOL)
% converts M into state space realization in upper block Hessenberg form.
% The controllable part of the system is returned as the 
% system {F,G,H}, where F is an upper block Hessenberg matrix and G has 
% nonzero elements in its first M rows only, where M is the rank of B. 
% The vector BSIZES contains the sizes of the different blocks of the 
% matrix F.
%
% A tolerance TOL may be specified as an additional input argument.
%
% For strict properness, see FRAC/ISPROPER.
% See also FRAC/BHF, BHF.

%        Author: J. Jezek, 28-Jul-2001
%        Copyright (c) 2001 by Polyx, Ltd.

if nargin==3 & ~isempty(tol),
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
else
   tol = [];
end;

[A,B,C,D] = abcd(M,tol);
if ~isempty(D) & any(any(D~=0)),
   error('Argument is not convertible to a strictly proper fraction.');
end;
[F,G,H,bsizes] = bhf(A,B,C,tol);

%end .. @pol/bhf

   