function Al = tril(A,k,tol);
%TRIL   Lower triangular part of polynomial
%
% TRIL(P) is the lower triangular part of X.
% TRIL(P,K) is the elements on and below the K-th diagonal
% of P .  K = 0 is the main diagonal, K > 0 is above the
% main diagonal and K < 0 is below the main diagonal.
%
% TRIL(P,K,TOL) works with zeroing specified by the input 
% tolerance TOL.
%
% See also POL/TRIU.

%       Author(s): M. Hromcik, M. Sebek 16-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 06-Mar-1998 10:10:34   $
%       $          3.0 $  $Date: 12-Jul-2000 14:45:00   $                  
%                         $Date: 24-Jun-2001  J.Jezek   $
% Effect on other properties:
% Al.u: UserData are deleted.

global PGLOBAL;

switch nargin 
   case 1, 
     	k = 0; 
	   tol = PGLOBAL.ZEROING;
   case 2,
     if isempty(k), 
      Al = pol([]);
      warning('K-th diagonal input must be an integer scalar.');
      return;
     end;
     if ~isnumeric(k) | length(k)>1,
        error('Invalid 2nd argument; must be numeric scalar.');
     end;
     tol = PGLOBAL.ZEROING;
   case 3
     
     if ~isnumeric(tol) | length(tol)>1,
        error('Invalid 3rd argument; must be numeric scalar.');
     end;
     tol = PGLOBAL.ZEROING;
     
     if ~isnumeric(k) | length(k)>1,
        error('Invalid 2nd argument; must be numeric scalar.');
     end;
     tol = PGLOBAL.ZEROING;
      
     if isempty(k), 
      Al = pol([]);
      warning('K-th diagonal input must be an integer scalar.');
      return;
     end;
          
end; 	%switch

Al.d = A.d;
Al.s = A.s;
Ac = A.c;
Alc = zeros(size(A.c));
for i=1:A.d+1,
  Alc(:,:,i) = tril(Ac(:,:,i),k);
end;
Al.c = Alc;
Al.v = A.v;
Al.h = A.h;
Al.u = [];
Al.version = 3.0;

Al = class(Al,'pol');  

% clearing:
Al = pclear(Al);

%end .. @pol/tril
  