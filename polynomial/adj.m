function [Aadj,dtA] = adj(A,arg2,arg3);
%ADJ   Adjoint of constant
%
% The command
%    AADJ = ADJ(A) 
% computes the adjoint of the square constant matrix A. 
% The adjoint is defined as AADJ(I,J) = DET(AIJ)*(-1)^(I+J), 
% where AIJ is A with its I-th column and J-th row omitted. 
% For any square matrix A, AADJ * A = A * AADJ = DET(A) * EYE(N), 
% where N is the size of A. If A is nonsingular then 
%    INV(A) = ADJ(A)/DET(A).
%
% The command
%    [AADJ,DETA] = ADJ(A) 
% computes both the adjoint and the determinant of A.
% If DETA is nonzero then AADJ/DETA is the inverse of A.
%
% This macro exists only for completeness.

%       Author: J. Jezek, 24-Feb-2003
%       Copyright(c) 2003 by Polyx, Ltd.

if nargin==3,
   if isempty(arg3) | isa(arg3,'char') | ...
         (isa(arg3,'double') & length(arg3==1)),
   else
      error('Invalid 3rd argument.');
   end;
else arg3 = [];
end;
if nargin>=2,
   if isempty(arg2) & isa(arg2,'char') | ...
         (isa(arg2,'double') & length(arg2)==1),
   else
      error('Invalid 2nd argument.');
   end;
else arg2 = [];
end;

Aadj = 0; dtA = 0;
if nargout<2,
   eval('Aadj = adj(pol(A),arg2,arg3);', ...
      'error(peel(lasterr));');
   Aadj = double(Aadj);
else
   eval('[Aadj,dtA] = adj(pol(A),arg2,arg3);', ...
      'error(peel(lasterr));');
   Aadj = double(Aadj); dtA = double(dtA);
end;

%end .. adj

