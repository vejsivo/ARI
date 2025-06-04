function R = minus(P,Q,arg3,arg4);
%MINUS(-)     Subtract fractions
%           R = P-Q   or   R = MINUS(P,Q)
%
% For fractions P,Q, the command returns fraction R
% whose (i,j)-th entry is  P(i,j)-Q(i,j) .
%
% Matrices P,Q must have the same dimensions, unless one of them is
% scalar. In such a case, this scalar is subtracted with every entry
% of the second matrix.
%
% The variable symbols of P,Q should be the same. When not, a warning
% is issued and the symbols are changed to the standard one. However, 
% if one symbol is 'z' and the second 'z^-1' then the symbols play
% a role. The resulting symbol is taken form P, no warning being issued.
%
% An optional numerical input argument can specify the zeroing
% tolerance to be used instead of the standard one.
%
% When both operands are left-den or right-den fractions, an
% optional 'ldf' or 'rdf' input argument can specify the
% required class of the result.
% 
% See also FRAC/UMINUS, RDF/PLUS, LDF/PLUS, MDF/PLUS, SDF/PLUS.

%       Author:  J. Jezek  26-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 25-Jul-2002 $

na = nargin;
if na<2,
   error('Not enough input arguments.');
end;

Pcl = class(P); Qcl = class(Q);
if strcmp(Pcl,'frac') | strcmp(Qcl,'frac'),
   error('Function ''minus'' not defined for variables of class ''frac''.');
end;

if na==2,
   eval('R = plus(P,-Q);','error(peel(lasterr));');
else
   if ~isa(arg3,'double') & ~isa(arg3,'char'),
      error('Invalid 3rd argument.');
   end;
   if na==3,
      eval('R = plus(P,-Q,arg3);','error(peel(lasterr));');
   else
      if ~isa(arg4,'double') & ~isa(arg4,'char'),
         error('Invalid 4th argument.');
      end;
      eval('R = plus(P,-Q,arg3,arg4);','error(peel(lasterr));');
   end;   
end;

%end .. @frac/minus
