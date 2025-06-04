function [r,rts] = isprime(P,arg2,arg3)
%ISPRIME  Test whether a polynomial matrix is left or right prime
%
% The commands
%    ISPRIME(P)
%    ISPRIME(P,'l')
% return 1 if the polynomial matrix is left prime and 0 if it is not.
% In the form
%    ISPRIME(P,'r')
% it tests whether the matrix is right prime. In the form
%    [R,RTS] = ISPRIME(P)
% the optional output argument RTS returns the roots of P. 
% -If P is prime then RTS is empty. 
% -If P is not left prime and RTS is returned as NaN then P does not 
%  have full row rank. 
% -If P is not right prime and RTS is returned as NaN then P does not 
%  have full column rank.
%
% As last input argument an optional tolerance TOL may be included. This
% is used to test whether P has any roots. Its default value is the global
% zeroing tolerance.
%
% See also GLD, GRD.

%       Author(s): H. Kwakernaak, M. Hromcik
%       Copyright (c) by PolyX, Ltd.
%       Modified by M. Hromcik, 21-June-2001: logical outputs
%                   J. Jezek,   02-July-2001: cosmetics

% Initializations

global PGLOBAL

eval('P = pol(P);', 'error(peel(lasterr));');

if nargin == 1
   option = 'l'; 
   tol = PGLOBAL.ZEROING;
elseif nargin == 2
   if ischar(arg2)
      option = arg2;
      tol = PGLOBAL.ZEROING;
   elseif isnumeric(arg2)
      option = 'l';
      tol = arg2;
   else
      error('Invalid 2nd argument.');
   end
elseif nargin == 3
   if ~ischar(arg2),
      error('Invalid 2nd argument.');
   end;
   if ~isnumeric(arg3),
      error('Invalid 3rd argument.');
   end;
   option = arg2; tol = arg3;
end

if length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;
if strcmp(option,'r')
   P = P.';
elseif ~strcmp(option,'l')
   error('Invalid command option.');
end

n = size(P,1); m = size(P,2);

% Compute a set of roots of the "squared down" matrix

detP1 = det(rand(n,n)*P*rand(m,n));
if norm(detP1,1) < tol*norm(P,1)
   r = logical(0); rts = NaN; return
end
if deg(pzer(detP1,tol)) == 0
   r = logical(1); rts = []; return
end
rts1 = roots(detP1);


% Test if any of these are the actual roots of P

detP2 = det(rand(n,n)*P*rand(m,n));
rts = [];
for i = 1:length(rts1)
   if abs(polyval(detP2,rts1(i))) < tol*norm(detP1,1)
      rts = [rts; rts1(i)];
   end
end


% Finalize

r = (length(rts) == 0);

%end .. @pol/isprime
