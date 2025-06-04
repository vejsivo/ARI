function[V,T,rankV] = hhl(U,tol)
%HHL
%
% This macro carries out a sequence of Householder transformations
% from the left to bring the constant input matrix in the form
% 
%           | * * * * |
%           | 0 * * * |
%           | 0 0 * * |
%           | 0 0 0 * |
%           | 0 0 0 0 |
%           | ....... |
%           | 0 0 0 0 |
%
% The command
%
%    [V,T,rankV] = HHL(U,tol)
%
% returns the output matrix V, the transformation matrix T with 
% V = T*U and the rank of the input matrix U.

% Author: Rens C.W. Strijbos, July 2, 1998.
% Copyright 1998 by Polyx, Ltd.


[ru,cu]=size(U);
   
T=eye(ru);
V=U;
rankV=0;

if nargin == 1
   tol = 1e5*eps;
end
rowoffset=0;
for i=1:cu
   u=V(i+rowoffset:ru,i);
   lengthu=length(u);
   k = norm(u);
   if k>tol & lengthu>1
      if u(1)>0
      	k = -k;
      end
      w = u;
      w(1) = w(1)-k;
      W = w/norm(w);
      Tloc = eye(lengthu)-2*W*W';
      T(i+rowoffset:ru,:)=Tloc*T(i+rowoffset:ru,:);
      V=T*U;
      rankV=rankV+1;
   else
      if k>tol
	 rankV=rankV+1;
	 break;
      end
      if length(u)>=1
	 rowoffset=rowoffset-1;
      end
   end
end

%end .. private/hhl

