function dtA = det(A,arg1,arg2)
%DET    Determinant of two-sided polynomial
%
% The command
%    DET(A) 
% computes the determinant of a square tsp matrix.
% The command
%    DET(A,METHOD)  
% allows the user to specify the method used for the
% computation of the determinant of the corresponding
% polynomial:
%   'fft': computation by interpolation and fast Fourier 
%          transformation. This is the default method.
%   'eig': computes the determinant as the characteristic
%          polynomial of the block companion matrix  
%          corresponding to the polynomial.  		
% The commands
%    DET(A,TOL) or DET(A,METHOD,TOL) 
% work with zeroing specified by the the input tolerance TOL.
%
% See also TSP/ADJ, TSP/RANK.

%      Author:  J.Jezek, 11-8-99
%      Copyright (c) 1998 by Polyx, Ltd.
%      $Revision: 3.0 $  $Date: 29-Sep-1999  13:00:00  $
%                        $Date: 01-Dec-1999  12:00:00  $
%                        $Date: 22-May-2000  12:00:00  $
%                        $Date: 31-Oct-2000  13:00:00  $

ni = nargin;
if ni>=2,
   if ~isa(arg1,'double') & ~isa(arg1,'char'),
      error('Invalid 2nd argument.');
   end;
end;
if ni==3,
   if ~isa(arg2,'double') & ~isa(arg2,'char'),
      error('Invalid 3rd argument.');
   end;
end;

PP = 0;
if ni==1,
   eval('PP = det(A.p);','error(peel(lasterr));');
elseif ni==2,
   eval('PP = det(A.p,arg1);','error(peel(lasterr));');
else
   eval('PP = det(A.p,arg1,arg2);','error(peel(lasterr));');
end;

dtA = tsp(PP); dtA.h = A.h;
dtA = shift(dtA,A.o*A.s(1));

%end .. @tsp/det