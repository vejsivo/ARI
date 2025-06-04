function dtA = det(A, arg2, arg3);
%DET   Determinant of polynomial
%
% The command
%    DET(A) 
% computes the determinant of a square polynomial matrix.
% The command
%    DET(A,METHOD)  
% allows the user to specify the used method:
%   'fft': computation by interpolation and fast Fourier 
%          transformation. This is the default method.
%   'eig': computes the determinant as the characteristic
%          polynomial of the block companion matrix  
%          corresponding to A.  		
% The commands
%    DET(A,TOL) or DET(A,METHOD,TOL) 
% work with zeroing specified by the the input tolerance TOL.
%
% See also POL/ADJ, POL/RANK, POL/ROOTS.

%	     Author(s): M. Hromcik, M. Sebek 20-5-98
%       Copyright (c) 1998 by ..., Ltd.
%       $Revision: 2.0 $  $Date: 20-May-1998 10:28:34   $
%       $Revision: 3.0 $  $Date: 11-Aug-1999   J.Jezek  $
%                         $Date: 08-Jul-2001   J.Jezek  $
%                         $Date: 15-Mar-2002   J.Jezek  cosmetics $

global PGLOBAL;

met = 'fft';
tol = PGLOBAL.ZEROING;

na = nargin;
switch na,
case 1,
case 2,
   if isempty(arg2), na = 1;
   end;
case 3,
   em2 = isempty(arg2); em3 = isempty(arg3);
   if em2 & em3, na =1;
   elseif ~em2 & em3, na = 2;
   elseif em2 & ~em3, na = 2; arg2 = arg3;
   end;    
end;

switch na,
case 1,
case 2,
   if isa(arg2,'double'), tol = arg2;
   elseif isa(arg2,'char'), met = arg2;
   else
      error('Invalid 2nd nonempty argument.');
   end;
case 3,
   if isa(arg2,'double') & isa(arg3,'char'),
      tol = arg2; met = arg3;
   elseif isa(arg2,'char') & isa(arg3,'double'),
      met = arg2; tol = arg3;
   else
      error('Invalid 2nd or 3rd argument.');
   end;
end;
if length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;

PI = pi;
Ac = A.c;
Asc = A.s;
As = Asc(1);
Ad = A.d;

if Asc(2)-As,
   error('Matrix is not square.');
end;

if As == 1, 
   dtA = A;
   return;
end;  

% expected degree of determinant:
coldegA = deg(A,'col'); rowdegA = deg(A,'row');
degofdet = min(sum(coldegA),sum(rowdegA));

% constant polynomial matrix:
if degofdet <= 0,
   dtA.d = degofdet;
   dtA.s = [1 1];
   Ac = cat(3, Ac, zeros([Asc 2]) );
   dtA.c = det(Ac(:,:,1));
   dtA.v = A.v;
   dtA.h = A.h;
   dtA.u = [];
   dtA.version = 3.0;
   dtA = class(dtA, 'pol');
   % no zeroing
   dtA = pclear(dtA);    %  added by J.Jezek  July 2001
   return
end;

dtA.d = degofdet;
dtA.s = [1 1];

switch met,
        	
case 'fft',		% MODIFIED INTERPOLATION with FFT & IFFT

	detA = zeros(1,degofdet+1);
	fq2 = floor(degofdet/2);
	
	% FFT of polynomial matrix A(s):
	Y = fft( Ac, degofdet+1, 3);
				
	% compute det(A(s)) at Fourier points:
	
  	if isreal(Ac),			% quicker for REAL POLYNOMIAL matrix

   % compute det( A(si) )
      for i = 1:fq2,			% lower half of Fourier points
         detAi = det(Y(:,:,1+i));		% = det( A(si) )
         detA(1+i) = detAi;
         detA(degofdet+2-i) = conj(detAi);
      end;

      % compute det(A(1)):
      detA(1) = det(Y(:,:,1)); 		% = det( A(1) )

      % even number of points - compute det(A(-1)):
      if rem(degofdet,2),		
         detA(fq2+2) = det(Y(:,:,fq2+2));	% = det( A(-1) )
      end;
	  
      % recover det(A(s)) (inverse FFT):
      detA = real( fft( conj(detA) ) ) ./ (degofdet+1);
	  	
   else				% COMPLEX POLYNOMIAL matrix

      % compute det( A(si) ):
      for i = 1:degofdet+1,			% all Fourier points
         detAi = det(Y(:,:,i));	  	% = det( A(si) )
         detA(i) = detAi;
      end;
		
      % recover det(A(s)) (inverse FFT):
      detA = conj( fft( conj(detA) ) ) ./ (degofdet+1);
   end;  
  	      	
case 'eig',		% GENERELIZED EIGENVALUES
   z = roots(A, 'eig', 'all');
   if sum( isnan( real(z) ) ),
      dtA = pol(0);
      return;
   end;
	
   kinf = find(~isinf(z));
   if isempty(kinf), 
      z = [];
   else,
      z = z(kinf);
   end;
	
   % first zeroing to cut "almost infinite" eigenvalues
   if tol > 0
      n = find(abs(z) > 1/tol);
   else
      n = [];
   end;

   if ~isempty(n),
      z(n) = [];
   end;
   if isempty(z)
      zet = 0;
   else
      % choice of zet which is not a zero
      f = 0;
      while ~isempty(f),
         zet = randn / norm(Ac(:,:));
         f = find( abs(z - zet) < eps );
      end;
   end;

   detA( degofdet-length(z)+1:degofdet+1 ) = poly(z);
     
   Az = polyval(A,zet);
   Azd = det(Az);
   if Azd == 0
      dtA = pol(0);
      return
   end
   dz = polyval(detA, zet);
   detA = detA * Azd/dz;
	
   if isreal(Ac),
      detA = real(detA);
   end;

   detA = fliplr(detA);

otherwise 
   error('Invalid command option.');   

end;		%switch method

dtA.c = permute(detA, [3,1,2]);
dtA.v = A.v;
dtA.h = A.h;
dtA.u = [];
dtA.version = 3.0;

% packaging:
dtA = class(dtA,'pol');

% zeroing
me = min(abs(nonzeros(Ac)));	% =?= min(abs(nonzeros(Ac))) ^ As;
if ~isempty(me),
  dtA = pzer(dtA,tol*me / As);
end;

%end .. @pol/det  	    
