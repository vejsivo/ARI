function S = sylv(A,k,string)
%SYLV  Create Sylvester matrix of polynomial object.
% 
% S = SYLV(A,K,STRING)  creates constant Sylvester resultant
% matrix S of order K  from polynomial matrix A.
% If
%        A = A0 + A1*v + A2*v^2 + ... + Ad*v^d
% then
%        S = [ A0 A1 A2 ... Ad  0  0  ...  0 |   block row 1.
%            | 0  A0 A1  A2 ... Ad 0  ...  0 |   block row 2.
%            | 0  0  ...    ... ...   ...  0 |   ...
%            | 0  0  ...    ... ...A0 ...  Ad]   block row K+1.
%
% K is the number of zeros blocks in each block row; its default 
% value is K = D.
% The input argument STRING can be:
% 'row' or missing - row Sylvester matrix
% 'col' - column Sylvester matrix

%       Author(s):  S. Pejchova, M. Sebek 25-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 30-Sep-1998 16:54:34   $
%       $Revision: 3.0 $  $Date: 08-Aug-2000  S. Pejchova   $
%                         $Date: 08-Aug-2001  J. Jezek  arg check $   
% Effect on other properties:
% S is a standard Matlab matrix.

ni = nargin;
if ni==0,
   error('Not enough input arguments.');
elseif ni==1,
   eval('S = sylv(pol(A));', 'error(peel(lasterr));');
else
   if ~isa(k,'double') & ~isa(k,'char'),
      error('Invalid 2nd argument.');
   end;
   if ni==2,
      eval('S = sylv(pol(A),k);','error(peel(lasterr));');
   else
      if ~isa(string,'char') & ~isa(string,'double'),
         error('Invalid 3rd argument.');
      end;
      eval('S = sylv(pol(A),k,string);','error(peel(lasterr));');
   end;   
end;

%end .. sylv
