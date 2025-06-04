function K = pkron(A,B)
%PKRON   Kronecker tensor product 
%        modified version of MATLAB function KRON to handle 3D arrays
%        which is usefull for the Polynomial Toolbox functionality
%        Necessary since 2010 when MATLAB KRON was limitted to 2D arrays 
%
%   For details, see the help of MATLAB KRON function
%   Works only with full arrays. With sparse arrays, use the original
%   MATLAB KRON

%   Copyright 2010 PolyX Ltd. 
%   $ M. Sebek $ $Date: 03-Dec-2011 $

[ma,na] = size(A);
[mb,nb] = size(B);

% if ~issparse(A) && ~issparse(B)  --- error message to be added

   [ia,ib] = meshgrid(1:ma,1:mb);
   [ja,jb] = meshgrid(1:na,1:nb);
   K = A(ia,ja).*B(ib,jb);
   
% end.. pkron
