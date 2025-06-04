function  T = cat(dim,varargin)
%CAT      Concatenation of fractions
%
%  CAT(1,A,B)  acts as  VERTCAT(A,B)  or  [A;B]
%  CAT(2,A,B)       as  HORZCAT(A,B)  or  [A,B]  or  [A B]
% 
%  As the concatenation operation for fractions includes
%  computation, namely solving a polynomial equation, the tolerance
%  TOL may be given as an optional input parameter. However, to be
%  recognizable from other ones, this argument must have a form
%  of character string, e.g.  CAT(1,A,B,'10^-6')

%  See also RDF/HORZCAT, LDF/HORZCAT, MDF/HORZCAT, SDF/HORZCAT,
%  RDF/VERTCAT, LDF/VERTCAT, MDF/VERTCAT, SDF/VERTCAT.

%         Author:  J. Jezek  26-Jan-2000
%         Copyright (c) 2000 by Polyx, Ltd.
%         $ Revision $  $ Date 25-Jul-2002 $

if nargin<2, error('Not enough input arguments.');  end;
if ~isa(dim,'double') | length(dim)~=1 | ...
      ~(dim==1 | dim==2),
   error('Invalid 1st argument; must be 1 or 2.');
end;

ni = length(varargin);
cl = class(varargin{1});
if strcmp(cl,'frac'),
   error('Dummy function for HORZCAT and VERTCAT.');
end;

ni = length (varargin);
switch dim,
case 1, 
   eval('T = vertcat(varargin{1:ni});','error(peel(lasterr));');
case 2,
   eval('T = horzcat(varargin{1:ni});','error(peel(lasterr));');
end;

%end .. @frac/cat
