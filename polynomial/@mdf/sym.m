function Asym = sym(Afrac, var)
%SYM  Convert a matrix-den matrix fraction to symbolic format.
%
% The commands
%    AS = SYM(AF)
%    AS = SYM(AF,VAR)
% convert a matrix-den matrix fraction AF of the Polynomial Toolbox
% into a symbolic matrix AS of the Symbolic Toolbox.
%
% Entries in AS are ratios of polynomials in the indeterminate VAR.
% If VAR is unspecified then the variable symbol of symbolic matrix AS
% is the same as that of fraction matrix AF.
%
% For backward conversion, see MDF/MDF.

% Author: D. Henrion, July 28, 2000.
% Copyright 2000 by Polyx, Ltd.
% Modified by D. Henrion, July 31, 2000.
%          by J. Jezek,   Nov  06, 2000.
%          by J. Jezek,   Oct  14, 2002.

if nargin == 1,
 var = Afrac.frac.v;
elseif ~isa(var, 'char'),
 error('Invalid 2nd argument; must be a string.');
end;

den = Afrac.frac.den;
num = Afrac.frac.num;

densym = 0; numsym = 0;
eval('densym = sym(den,var);','error(peel(lasterr));');
eval('numsym = sym(num,var);','error(peel(lasterr));');
eval('Asym = numsym ./ densym;','error(peel(lasterr));');

%end .. @mdf/sym
