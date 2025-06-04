function H = mldivide(G,F,arg3,arg4)
%MLDIVIDE    Matrix left divide fractions
%               H = G\F
%
% The command   H = G\F  or  H = MLDIVIDE(G,F) is the matrix left
% division of fractions F and G. Any scalar may be divided by
% anything or into anything. Otherwise, G must be square nonsingular
% and its dimension must equal to the number of rows of F.
%
% An optional numerical input argument may specify the zeroing
% tolerance to be used.
%
% When both operands are left-den or right-den fractions, an
% optional 'rdf' or 'ldf' input argument may specify the required
% class of the result.
%
% See also FRAC/MRDIVIDE, RDF/RDIVIDE, LDF/RDIVIDE, MDF/RDIVIDE,
%          SDF/RDIVIDE, RDF/LDIVIDE, LDF/LIVIDE, MDF/LDIVIDE, SDF/LDIVIDE.

%         Author:  J. Jezek  18-Nov-1999
%         Copyright(c) 1999 by Polyx, Ltd.
%         $ Revision $   $ Date 25-Jan-2002 $
%                        $ Date 25-Jul-2002 $

na = nargin;
if na<2,
   error('Not enough input arguments.');
end;

Fcl = class(F); Gcl = class(G);
if strcmp(Fcl,'frac') | strcmp(Gcl,'frac'),
   error('Function ''\'' not defined for variables of class ''frac''.');
end;

if ~isa(F,'frac') & ~isa(F,'pol') & ~isa(F,'tsp') & ...
      ~(isa(F,'double') & ndims(F)==2),
   eval('F = pol(F);', ...
      'eval(''F = sdf(F);'',''error(''''Invalid 2nd argument.'''');'');');
end;
if ~isa(G,'frac') & ~isa(G,'pol') & ~isa(G,'tsp') & ...
      ~(isa(G,'double') & ndims(G)==2),
   eval('G = pol(G);', ...
      'eval(''G = sdf(G);'',''error(''''Invalid 1st argument.'''');'');');
end;

if na==2,
   eval('H = mtimes(inv(G),F);','error(peel(lasterr));');
else
   if ~isempty(arg3) & (~isa(arg3,'double') & ~isa(arg3,'char')),
      error('Invalid 3rd argument.');
   end;
   if na==3,
      eval('H = mtimes(inv(G),F,arg3);','error(peel(lasterr));');
   else
      if ~isempty(arg4) & (~isa(arg4,'double') & ~isa(arg4,'char')),
         error('Invalid 4th argument.');
      end;
      eval('H = mtimes(inv(G),F,arg3,arg4);','error(peel(lasterr));');
   end;
end;

%end .. @frac/mldivide

