function A = spf(B,varargin)
%SPF  Two-sided polynomial spectral factorization
%
% If B is a para-Hermitian two-sided polynomial matrix,
% that is,   B'(z^-1) = B(z) , positive definite on the unit circle,
% the command
%    A = spf(B)
% returns a stable polynomial A(z), satisfying
%    A'(z^-1)*A(z) = B(z) .
% The result A(z) has its absolute coefficient upper triangular
% with positive diagonal entries.
%
% The factorization is performed with a Newton-Raphson iterative
% scheme using the macro AXXAB.
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
% The iterative scheme is stopped when norm(A'*A-B) is less than
% TOL*norm(B)*100. The tolerance TOL is also used in the macro AXXAB.
%
% The commmand
%    SPF(B,'syl') 
% is based on a Sylvester matrix algorithm. This the default method.
% The commmand
%    SPF(B,'red') 
% is based on polynomial reductions.

%        Author:  J. Jezek  11-8-1999
%        Copyright (c) 1999 by Polyx, Ltd.
%        $ Revision $  $ Date 24-May-2000 $

eval('B = tsp(B);','error(peel(lasterr));');
if B~=B',
   error('Polynomial matrix is not para-Hermitian.');
end;
Bd = deg(nneg(B));
Bnn = nneg(shift(B,Bd));
symbol(Bnn,'q');

ni = length(varargin);
if ni>=1,
   for i = 1:ni;
      arg = varargin{i};
      if ~isempty(arg),
         if ~isa(arg,'char') & ~isa(arg,'double'),
            error(['Invalid ',nth(i+1),' argument.']);
         end;
      end;
   end;
end;

eval('A = spf(Bnn,varargin{:});','error(peel(lasterr));');
symbol(A,'z'); A.h = B.h;

%end .. @tsp/spf





