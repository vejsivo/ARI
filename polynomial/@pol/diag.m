function D = diag(A,k,tol);
%DIAG  Diagonal matrix or diagonal of polynomial
%
% Let Av be a polynomial vector with N components. Then
%    DIAG(Av,K)  
% is a square polynomial matrix of dimensions (N+ABS(K))-by-(N+ABS(K))
% with the elements of Av on the K-th diagonal. K = 0 is the main diagonal, 
% K > 0 is above the main diagonal, and K < 0 is below the main diagonal.
% DIAG(Av) is the same as DIAG(Av,0).
%
% Let Am be a polynomial matrix. Then
%    DIAG(Am,K)  
% is a column polynomial vector formed from the elements of the K-th 
% diagonal of Am. DIAG(Am) is the main diagonal of Am. 
%
% DIAG(DIAG(Am)) is a diagonal matrix.
%
% See also POL/TRIL, POL/TRIU.

%       Author(s): S. Pejchova, M. Sebek 20-3-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 23-Feb-1999 14:31:11   $
%       $Revision: 3.0 $  $Date: 11-Aug-1999 12:00:00   J. Jezek  $
%                         $Date: 26-May-2000 15:00:00   J. Jezek  $
%                         $Date: 13-Oct-2002 J.Jezek warning $
%                         $Date: 22-Oct-2002  S. Pejchova -warning row 47. $  
%                         $Date: 28-Feb-2003 J.Jezek tol $
%
% Effect on other properties:
% D.u: UserData are deleted.

if nargin == 1, k = 0;
elseif ~isa(k,'double') | length(k)~=1 | ~isreal(k) | round(k)~=k,
   error('Invalid 2nd argument; must be scalar integer.');
end;
if nargin==3,
   if ~isa(tol,'double'),
      error('Invalid 3rd argument.');
   end;
end;
      
A = pol(A);
Ac = A.c; Ad = A.d; var = A.v;
D.d = Ad; Dc = [];
eval('Dc = diag(zeros(A.s(1),A.s(2)),k);','error(lasterr);');
D.s = size(Dc);
if ~isempty(Dc),
   if Ad > 0, Dc = repmat(Dc,[1 1 Ad+1]); end;  
   sw=warning; warning off;
   for i=1:Ad+1,
       Dc(:,:,i) = diag(Ac(:,:,i),k);
   end;
   warning(sw);
else, D.d=[]; var='';
end;
D.c = Dc;
D.v = var;
D.h = A.h;
D.u = [];
D.version = 3.0;
D = class(D,'pol');
D = pclear(D);

%end .. @pol/diag
