function  P = cat(dim,varargin)
%CAT      Concatenation of polynomial matrices
%
%  CAT(1,A,B)  acts as  VERTCAT(A,B)  or  [A;B]
%  CAT(2,A,B)       as  HORZCAT(A,B)  or  [A,B]  or  [A B]
% 
%  See also POL/HORZCAT, POL/VERTCAT.

%         Author:  J. Jezek  11-8-99
%         Copyright 1999 by Polyx, Ltd.
%         $Revision 3.0$  $Date: 11-Oct-1999 12:00:00  $

if nargin<2, error('Not enough input arguments.');  end;
if ~isnumeric(dim), error('Invalid 1st argument; must be 1 or 2.'); end;

ni = length(varargin);
switch dim,
case 1,
   eval('P = vertcat(varargin{1:ni});','error(peel(lasterr));');
case 2,
   eval('P = horzcat(varargin{1:ni});','error(peel(lasterr));');
otherwise, error('Invalid 1st argument; must be 1 or 2.');
end;

%end .. @pol/cat
