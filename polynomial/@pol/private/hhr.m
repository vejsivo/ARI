function[V,T,rankV] = hhr(U,tol)
%HHR
%
% This macro carries out a sequence of Householder transformations
% from the right to bring the constant input matrix in the form
%
%      | 0 0 .......0 * * * * |
%      | 0 0 .......0 0 * * * |
%      | 0 0 .......0 0 0 * * |
%      | 0 0 .......0 0 0 0 * |
%
% The command
%
%    [V,T,rankV] = HHR(U,tol)
%
% returns the output matrix V, the transformation matrix T with 
% V = U*T and the rank of the input matrix U.

% Author: Rens C.W. Strijbos, July 2, 1998.
% Copyright 1998 by Polyx, Ltd.

[ru,cu]=size(U);
   
T=eye(cu);
V=U';
rankV=0;

if nargin == 1
   tol = 1e5*eps;
end
rowoffset=0;
for i=ru:-1:1
   inter=cu-ru+i+rowoffset;
   u=V(1:inter,i);
   lengthu=length(u);
   k = norm(u);
   if k>tol & lengthu>1
      if u(inter)>0
      	k = -k;
      end
      w = u;
      w(inter) = w(inter)-k;
      W = w/norm(w);
      Tloc = eye(lengthu)-2*W*W';
      T(1:inter,:)=Tloc*T(1:inter,:);
      V=T*U';
      rankV=rankV+1;
   else
      if k>tol
	 rankV=rankV+1;
	 break;
      end
      if length(u)>=1
	 rowoffset=rowoffset+1;
      end
   end
end
T=T';
V=V';

%end .. @pol/private/hhr
