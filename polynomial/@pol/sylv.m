function S = sylv(A,varargin)
%SYLV  Create the Sylvester matrix corresponding to a polynomial matrix
%
% The command
%    S = SYLV(A,K,STRING)  
% creates the Sylvester resultant matrix S of order K corresponding
% to the polynomial matrix A. If
%    A = A0 + A1*v + A2*v^2 + ... + Ad*v^d
% then
%    S = [ A0 A1 A2 ... Ad  0  0  ...  0 |   block row 1.
%        | 0  A0 A1  A2 ... Ad 0  ...  0 |   block row 2.
%        | 0  0  ...    ... ...   ...  0 |   ...
%        | 0  0  ...    ... ...A0 ...  Ad]   block row K+1.
%
% K is the number of zeros blocks in each block row. 
% Its default value is K = D.
%
% The input argument STRING may be:
%    'row' or missing - row Sylvester matrix
%    'col'            - column Sylvester matrix

%       Author(s):  S. Pejchova, M. Sebek 25-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 06-Oct-1998 09:45:34   $
%       $Revision: 3.0 $  $Date: 08-Aug-2000  S. Pejchova   $

% Effect on other properties:
% S is a standard Matlab matrix.

ni = nargin; k=[];  string=0;
narginchk(1,3);
% error(nargchk(1,3,ni));	%REMOVED IN NEW MATLABS

for ii=1:ni-1,
   args=varargin{ii};,
   if ischar(args),
      switch args,
      case 'row', 
      case 'col', string = 1;
      otherwise,  error('Invalid command option.');
      end;
   elseif isa(args,'double')
      if all(size(args)==1) & isreal(args) & ...
            isfinite(args) & args==fix(abs(args)) & args>=0,
         k = args;
      else
         error('Invalid order; must be a nonnegative integer.');
      end;
   else
      error(['Invalid ',nth(ii+1),' argument.']);
   end;
end;

[As1,As2,Ad]=size(A);  Ac=A.c;
if ~size(Ac,3),
   Ad=0; Ac=zeros(As1,As2);  string = 0;
end;
if isempty(k),  k = Ad;  end;
if ~string,
   RAo = Ac(:,:);
   S = zeros((k+1)*As1,(Ad+k+1)*As2);
   S(1:As1,1:(Ad+1)*As2) = RAo;
   for d = 1:k,
      S((d*As1+1):(d+1)*As1,(d*As2+1):(As2*(Ad+d+1))) = RAo;
   end;
else,
   RAo = reshape(permute(Ac,[1 3 2]),[(Ad+1)*As1, As2]);
   S = zeros((Ad+k+1)*As1,(k+1)*As2);
   S(1:(Ad+1)*As1,1:As2) = RAo;
   for d = 1:k,
      S((d*As1+1):(As1*(Ad+d+1)),(d*As2+1):(d+1)*As2) = RAo;
   end;
end;

%end .. @pol/sylv
