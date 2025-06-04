function Y = value(F,X)
%VALUE    Value of right-denominator fraction
%
% For right-den fraction F, the command
%  Y = VALUE(F,X)  returns matrix Y whose
% every entry Y(i,j) is the value of function
% F(i,j) of the argument X(i,j). Matrices F,X
% must have the same dimensions unless one of
% them is scalar. In such a case, the scalar
% operates with every entry of the other matrix.
%
% When F is empty or constant, the function does not
% depend on X. In such a case, X may be any scalar or
% any compatible matrix or it may be omitted.
%
% Matrix X or some entries of X may be infinite.
%
% See also POL/VALUE, POL/MVALUE, RDF/MVALUE.

%       Author: J.Jezek, 17-Dec-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 27-Apr-2003 $

ni = nargin;
if ni==2,
   if ~isa(X,'double') | ndims(X)~=2,
      error('Invalid 2nd argument.');
   end;
   eval('F = mdf(F); Y = value(F,X);', ...
      'error(peel(lasterr));');
else
   eval('F = mdf(F); Y = value(F);', ...
      'error(peel(lasterr));');
end;   

%end .. @rdf/value
