function [F,G,H,bsizes] = bhf(M,tol)
%BHF   Convert fraction into state space in block Hessenberg form
%
% The function
%    [F,G,H,BSIZES] = BHF(M,TOL)
% converts a strictly proper fraction M into state space realization
% in upper block Hessenberg form.
% The controllable part of the system is returned as the 
% system {F,G,H}, where F is an upper block Hessenberg matrix and G has 
% nonzero elements in its first M rows only, where M is the rank of B. 
% The vector BSIZES contains the sizes of the different blocks of the 
% matrix F.
%
% A tolerance TOL may be specified as an additional input argument.
%
% See also BHF.

%        Author: J. Jezek, 28-Jul-2001
%        Copyright (c) 2001 by Polyx, Ltd.
%        $ Revision $  $ Date 25-Jul-2002 $

Mcl = class(M);
if strcmp(Mcl,'frac'),
   error('Invalid 1st argument.');
end;

if  nargin==1,
   if ~isproper(M,'stric')
      error('Fraction is not strictly proper.');
   end;
   [A,B,C,D] = abcd(M);
   [F,G,H,bsizes] = bhf(A,B,C);
else
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
   if ~isproper(M,'stric')
      error('Fraction is not strictly proper.');
   end;
   [A,B,C,D] = abcd(M,tol);
   [F,G,H,bsizes] = bhf(A,B,C,tol);
end;

%end .. @frac/bhf

   

