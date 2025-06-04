function Asym = sym(Atsp, var)
%SYM  Convert a two-sided matrix polynomial to symbolic format.
%
% The commands
%    AS = SYM(AT)
% convert a two-sided matrix polynomial AT of the Polynomial Toolbox
% into a symbolic matrix AS of the Symbolic Toolbox.

% Author: D. Henrion, July 28, 2000.
% Copyright 2000 by Polyx, Ltd.

Atsp = tsp(Atsp);
Asym = expand(sym(Atsp.p)*power(sym('z'),Atsp.o));

%end .. @tsp/sym


