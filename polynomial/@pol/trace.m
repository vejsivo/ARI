function At = trace(A, tol);
%TRACE  Sum of the diagonal elements of polynomial
%   
% TRACE(P) is the sum of the diagonal elements of P
% with zeroing activated through the global variable
% PGLOBAL.ZEROING.
%
% TRACE(P,TOL) works with zeroing specified by the input 
% tolerance TOL. 

%       Author(s): M. Hromcik, M. Sebek 16-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 27-Feb-1998 10:28:34   $
%       $Revision: 3.0 $  $Date: 16-Jun-2000 12:00:00  J. Jezek  $
%                         $Date: 28-Feb-2003           J. Jezek  $

% Effect on other properties:
% At.u: UserData are deleted.

global PGLOBAL;

if nargin == 1 | isempty(tol), 
   tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;

% The following statements are peculiar to be
% compatible with the peculiar standard MATLAB
% function TRACE for non-square matrices
Ad = A.d;
if isempty(Ad), Ad = -Inf; end;
Ac = A.c;
[Ac1,Ac2,Ac3] = size(Ac);
if Ac3==0, Ac = zeros(Ac1,Ac2,1); end;
Ats = size(trace(Ac(:,:,1)));
Atc = zeros(Ats);
At.d = Ad;
At.s = Ats;
if ~isfinite(Ad), Ad = -1; end
if Ad>=0,
  for i=1:Ad+1,
    Atc(:,:,i) = trace(Ac(:,:,i));
  end;
end; 
At.c = Atc;
At.v = A.v;
At.h = A.h;
At.u = [];
At.version = 3.0;

At = class(At,'pol');

me =  min(abs(nonzeros(A.c)));
% zeroing
if ~isempty(me),
   At = pzer(At,tol*me);
end;

%end .. @pol/trace
  