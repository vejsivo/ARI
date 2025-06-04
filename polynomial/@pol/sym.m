function Asym = sym(Apol, var)
%SYM  Convert a polynomial matrix to symbolic format
%
% The commands
%    AS = SYM(AP)
%    AS = SYM(AP,VAR)
% convert a polynomial matrix AP of the Polynomial Toolbox
% into a symbolic matrix AS of the Symbolic Toolbox.
%
% Entries in AS are polynomial expressions in the indeterminate VAR.
% If VAR is unspecified then the variable symbol of symbolic matrix AS
% is the same as that of polynomial matrix AP.

% Author: D. Henrion, October 14, 1998.
% Copyright 1998 by Polyx, Ltd.

global PGLOBAL;

eval('Apol = pol(Apol);', 'error(peel(lasterr));');

if nargin == 1,
 var = symbol(Apol);
elseif ~isa(var, 'char'),
 error('Invalid 2nd argument; must be a string.');
end;

As = Apol.s;
Asym = sym(zeros(As));

Ad = Apol.d;

if Ad >= 0, % non-zero polynomial matrix

 % Build the symbolic matrix from entries returned by Symbolic Toolbox
 % macro POLY2SYM.

 for row = 1:As(1),
  for col = 1:As(2),
   coef = Apol.c;
   entry = coef(row, col, :);
   if isempty(var),
    entry = poly2sym(fliplr(reshape(entry, 1, length(entry))));
   else
    entry = poly2sym(fliplr(reshape(entry, 1, length(entry))), var);
   end;
   if length(entry) == 1,
    Asym(row, col) = entry;
   end;
  end;
 end;

end;

%end .. @pol/sym




