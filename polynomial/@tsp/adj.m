function [Aadj,dtA] = adj(A,arg2,arg3)
%ADJ     Adjoint of two-sided polynomial
%
% The command
%    AADJ = ADJ(A) 
% computes the adjoint of the square tsp matrix A. 
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
% The commands
%    AADJ = ADJ(A,'INT')
%    [AADJ,DETA] = ADJ(A,'INT') 
% use fast Fourier transformation and interpolation for computing 
% AADJ and DETA. This is the default method.
%
% The commands
%    AADJ = ADJ(A,'DEF')
%    [AADJ,DETA] = ADJ(A,'DEF') 
% proceed according to the definition of the adjoint matrix and use 
% fast Fourier transformation and interpolation for computing the
% subdeterminants.
%
% The commands
%    ADJ(A, TOL)
%    ADJ(A, METHOD, TOL) 
% work with zeroing specified by the input tolerance TOL. 
% METHOD is 'INT' or 'DEF'. The default value of TOL is the global 
% zeroing tolerance.
%
% See also TSP/INV, TSP/DET.

%     Author: J. Jezek  11-8-99
%     Copyright (c) 1999 by Polyx, Ltd.
%     $Revision: 3.0 $  $Date: 29-Sep-1999  13:00:00  $
%                       $Date: 22-May-2000  12:00:00  $
%                       $Date: 31-Oct-2000  13:00:00  $

ni = nargin; no = nargout;
if ni>=2,
   if ~isa(arg2,'double') & ~isa(arg2,'char'),
      error('Invalid 2nd argument.');
   end;
end;
if ni==3,
   if ~isa(arg3,'double') & ~isa(arg3,'char'),
      error('Invalid 3rd argument.');
   end;
end;

if no==1,
   PP = 0;
   if ni==1,
      eval('PP = adj(A.p);','error(peel(laserr));');
   elseif ni==2,
      eval('PP = adj(A.p,arg2);','error(peel(lasterr));');
   else
      eval('PP = adj(A.p,arg2,arg3);','error(peel(lasterr));');
   end;
   Aadj = tsp(PP); Aadj.h = A.h;
   Aadj = shift(Aadj,A.o*(A.s(1)-1));
else
   PP = 0; QQ = 0;
   if ni==1,
      eval('[PP,QQ] = adj(A.p);','error(peel(lasterr));');
   elseif ni==2,
      eval('[PP,QQ] = adj(A.p,arg2);','error(peel(lasterr));');
   else
      eval('[PP,QQ] = adj(A.p,arg2,arg3);','error(peel(lasterr));');
   end;
   Aadj = tsp(PP); Aadj.h = A.h;
   Aadj = shift(Aadj,A.o*(A.s(1)-1));
   dtA = tsp(QQ); dtA.h = A.h;
   dtA = shift(dtA,A.o*A.s(1));
end;

%end .. @tsp/adj
