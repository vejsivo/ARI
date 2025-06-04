function [C,D] = compan(A,varargin)
%COMPAN  Create the block companion matrix of a polynomial matrix
%
% Given a square polynomial matrix A
%     A = A{0} + A{1}*v + A{2}*v^2 + ... + A{N}*v^N
% with nonsingular leading coefficient matrix A{N}, the commands
%    [C,D] = COMPAN(A)
%    [C,D] = COMPAN(A,'bot')
% return the "bottom" block companion matrix
%    C = [     0           I           0     ...         0      |
%        |     0           0           I     ...         0      |
%        |    ...         ...         ...    ...        ...     |
%        |     0          ...         ...    ...         I      |
%        | -A{N}\A{0}    -A{N}\A{1}   ...    ...   -A{N}\A{N-1} ]
% The matrix D is an identity matrix of the same dimensions as C.
%
% Similarly, the commands
%    [C,D] = COMPAN(A,'top')
%    [C,D] = COMPAN(A,'right')
%    [C,D] = COMPAN(A,'left')
% return "top", "right" and "left" versions of the block companion
% matrix, respectively.
%
% If, in addition, the string 'sep' is used as an input argument
% then the resulting block companion matrix is returned as two 
% separate matrices, for instance
%    C = [   0     I    0  ...    0    |   D = [ I   0   ...   0  |
%        |   0     0    I  ...    0    |       | 0   I   ...   0  |
%        |  ...   ...  ... ...   ...   |       |... ...  ...  ... |
%        |   0    ...  ... ...    I    |       | 0  ...   I    0  |
%        | -A{0} -A{1} ... ... -A{N-1} ]       | 0  ...   0   A{N}]
% Here A{N} may be singular.

%       Author(s):  S. Pejchova  01-10-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 30-Nov-1998 14:41:34   $
%       $Revision: 3.0 $  $Date: 28-Jul-2000 09:53:34  S.Pejchova  $
%                         $Date: 19-Nov-2002  S. Pejchova  - warning row 115.$  

% Effect on other properties:
% C and D are standard Matlab matrices.

global PGLOBAL; in_tol=0;
eval('in_tol=PGLOBAL.ZEROING;', 'painit; in_tol=PGLOBAL.ZEROING;');
ni=nargin;
narginchk(1,3);
% error(nargchk(1,3,ni));	%REMOVED IN NEW MATLABS

A=pol(A); s1='bot'; s2=1; C=[];  D=[];
for ii=1:ni-1,
    argm=varargin{ii};
    switch argm,
     case {'bot','top','right','left'},
       s1=argm;
     case 'sep',
       s2=0;
     otherwise,
       error('Invalid command option.');       
    end;
end;

[rA, cA,degA] = size(A);

if rA~=cA
   error('Matrix is not square.');
elseif (isempty(degA) | degA<=0),
   return;
end

Ac=A.c;  LA=Ac(:,:,degA+1);
if s2 & (abs(det(LA))<=in_tol)
   warning('Singular leading coefficient matrix.');
end
w1=warning; warning off;
if degA==1                                           %degree is 1
   C=-Ac(:,:,1);  D=LA;
   if s2
      switch s1,
      case {'bot','top'}, C = D\C;
      case {'right','left'}, C = C/D;
      end;
      D=eye(rA);
   end;
else,
   CZC=zeros((degA-1)*cA,cA);                         %zeros column block
   CZR=zeros(rA,(degA-1)*cA);                         %zeros row block
   CE=eye((degA-1)*cA);                               %identity block
   D=eye(degA*cA); Ac=Ac(:,:,1:degA);
   switch s1,
   case 'bot',
      CA=Ac(:,:);                                     %'bottom' block row
      if s2, CA=LA\CA;
      else, D(end-rA+1:end,end-rA+1:end)=LA;
      end
      C=[CZC,CE;-CA];
   case 'top',
      CA=flipdim(Ac,3); CA=CA(:,:);                   %'top' block row
      if s2, CA=LA\CA;
      else, D(1:rA,1:rA)=LA;
      end
      C=[-CA;CE,CZC];
   case 'right',
      CA=reshape(permute(Ac,[1 3 2]),[degA*rA, rA]);   %'right' block column
      if s2, CA=CA/LA;
      else, D(end-rA+1:end,end-rA+1:end)=LA;
      end
      C=[[CZR;CE],-CA];     
   case 'left',                                       %'left' block column
      CA=reshape(permute(Ac(:,:,end:-1:1),[1 3 2]),[degA*rA, rA]);   
      if s2, CA=CA/LA;
      else, D(1:rA,1:rA)=LA;
      end
      C=[-CA,[CE;CZR]];
   end;
end; %if degA==1
warning(w1);

%end .. @pol/compan
