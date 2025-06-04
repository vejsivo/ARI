function R = times(P,Q,arg3,arg4);
%TIMES(.*)    Element-wise multiply right-den fractions
%           R = P.*Q   or   R = TIMES(P,Q)
%
% For right-den fractions  P,Q  the command returns right-den fraction R
% whose (i,j)-th entry is  P(i,j)*Q(i,j)
%
% Matrices P,Q must have the same dimensions, unless one of them is
% scalar. In such a case, this scalar is multiplied with every entry
% of the second matrix.
%
% The variable symbols of P,Q should be the same. When not, a warning
% is issued and the symbols are changed to the standard one. However, 
% if one symbol is 'z' and the second 'z^-1' then the symbols play
% a role. The resulting symbol is taken form P, no warning being issued.
%
% An optional input argument TOL can specify the zeroing tolerance
% to be used instead of the standard one.
%
% The result is standardly of class 'rdf'. However, it is less time
% consuming to reach the result in class 'mdf'. For that, an optional
% input argument CLASS may be used. It must be a character string,
% possible values 'mdf', 'sdf', 'rdf' or 'ldf'.
%
% See also RDF/MTIMES.

%       Author:  J. Jezek  03-Feb-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 25-Apr-2000 $
%                     $ Date 01-Nov-2000 $

global PGLOBAL;

tol = PGLOBAL.ZEROING; type = 'rdf';
na = nargin;
if na<2,
   error('Not enough input arguments.');
elseif na>2 & ~isempty(arg3),
   if isa(arg3,'double'), tol = arg3;
   elseif isa(arg3,'char'), type = arg3;
   else error('Invalid 3rd argument.');
   end;
end;
if na>3 & ~isempty(arg4),
   if isa(arg4,'double'), tol = arg4;
   elseif isa(arg4,'char'), type = arg4;
   else error('Invalid 4th argument.');
   end;
end;

cop = PGLOBAL.COPRIME; PGLOBAL.COPRIME = 'cop'; R = 0;
eval('P = mdf(P); Q = mdf(Q);','error(peel(lasterr));');
eval('R = times(P,Q,tol);','error(peel(lasterr));');
PGLOBAL.COPRIME = cop;

if strcmp(type,'rdf'), R = rdf(R,tol);
elseif strcmp(type,'ldf'), R = ldf(R,tol);
elseif strcmp(type,'sdf'), R = sdf(R,tol);
elseif ~strcmp(type,'mdf'),
   error('Invalid command option.');
end;

%end .. @rdf/times
