function Au = triu(A,k,tol);
%TRIU   Upper triangular part of polynomial
%
%   TRIU(P) is the lower triangular part of X.
%   TRIU(P,K) is the elements on and above the K-th diagonal
%   of P .  K = 0 is the main diagonal, K > 0 is above the
%   main diagonal and K < 0 is below the main diagonal.
%   TRIU(P,K,TOL) works with zeroing specified by the input 
%   tolerance TOL.
%
%  See also POL/TRIL.

%       Author(s): M. Hromcik, M. Sebek 16-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 06-Mar-1998 10:11:34   $
%       $          3.0 $  $Date: 12-Jul-2000 14:45:00   $
%                         $Date: 24-Jun-2001   J.Jezek  $

% Effect on other properties:
% Au.u: UserData are deleted.

global PGLOBAL;

switch nargin 
   case 1, 
     	k = 0; 
	   tol = PGLOBAL.ZEROING;
   case 2,
     if isempty(k), 
      Au = pol([]);
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
      Au = pol([]);
      warning('K-th diagonal input must be an integer scalar.');
      return;
     end;
          
end; 	%switch

Au.d = A.d;
Au.s = A.s;
Ac = A.c;
Auc = zeros(size(A.c));
for i=1:A.d+1,
  Auc(:,:,i) = triu(Ac(:,:,i),k);
end;
Au.c = Auc;
Au.v = A.v;
Au.h = A.h;
Au.u = [];
Au.version = 3.0;

Au = class(Au,'pol');  

% clearing:
Au = pclear(Au);

%end .. @pol/triu
  