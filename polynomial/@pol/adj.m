function [Aadj, dtA] = adj(A, arg2, arg3);
%ADJ    Adjoint of polynomial
%
% The command
%    AADJ = ADJ(A) 
% computes the adjoint of the square polynomial matrix A. 
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
%    AADJ = ADJ(A,'int')
%    [AADJ,DETA] = ADJ(A,'int') 
% use fast Fourier transformation and interpolation for computing 
% AADJ and DETA. This is the default method.
%
% The commands
%    AADJ = ADJ(A,'def')
%    [AADJ,DETA] = ADJ(A,'def') 
% proceed according to the definition of the adjoint matrix and use 
% fast Fourier transformation and interpolation for computing the
% subdeterminants.
%
% The commands
%    ADJ(A, TOL)
%    ADJ(A, METHOD, TOL) 
% work with zeroing specified by the input tolerance TOL. 
% METHOD is 'int' or 'def'. The default value of TOL is the global 
% zeroing tolerance.
%
% See also POL/INV, POL/DET.

%	     Author(s): M. Hromcik, M. Sebek 20-5-98
%	     Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 16-Sep-1998 10:28:34   $
%       $Revision: 3.0 $  $Date: 11-Aug-1999 12:00:00   J.Jezek  $
%                         $Date: 24-Jun-2001 J.Jezek  cosmetics  $ 
%                         $Date: 15-Mar-2002 J.Jezek  cosmetics  $
%                         $Date: 28-Mar-2003 J.Jezek  warning    $
%                         $Date: 09-Aug-2005 M.Hromcik modif for M-R14-SP2 $

global PGLOBAL;

tol = PGLOBAL.ZEROING;
met = 'int';

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

if isempty(Ac),
   Aadj = pol(zeros(Asc));
   if isempty(Ad), dtA  = pol(1);
   else dtA = pol(0);
   end;
   return;
end;  

% zeroing and singularity tolerance:
me = min(abs(nonzeros(Ac)));

k_sing = 1;
CondBnd = 1e14;

% expected degrees:
degofdet = As*Ad;
fq2 = floor(degofdet/2);

Aadjc = zeros([Asc, degofdet+1]);
detA = zeros(1, degofdet+1);

% FFT of polynomial matrix A(s):
Y = 0;
eval('Y = fft( Ac, degofdet+1, 3 );', 'Y = Ac(:,:,1);'); 

% =================================
%       Direct interpolation
% =================================
if strcmp(met, 'int'),

   saved_w = warning; warning off;
   anywarn = 0;

   % compute det(A(si)) and Aadj(si) = inv(A(si)) * det(A(si)) at Fourier points:
   if isreal(Ac),		% quicker for REAL POLYNOMIAL matrix

      % compute Aadj(si):
      for i = 1:fq2,		% lower half of Fourier points
         Ai = Y(:,:,i+1);
         detAi = det(Ai);		% = det( A(si)  ) 
         detA(1+i) = detAi;
         detA(degofdet+2-i) = conj(detAi);

         [Ui,Si,Vi] = svd(Ai);
         condA = Si(1,1)/Si(As,As);
    
         if abs(condA) < CondBnd,		% well conditioned - A(si)
            Aadj = (detAi*Vi/Si)*Ui';	% = inv(A(si)) * det(A(si)) = Aadj(si)
         else			% ill conditioned or singular
            anywarn = 1;
            tl = As * max(Si(:)) * eps;
            if sum(diag(Si)>tl) >= As-1,	% rank is close to As-1
               Aadj = Vi(:,end) * Ui(:,end)';  	% outter product of right and left null spaces    
               [YY,J] = max(abs(Aadj));
               [Amax,jj] = max(YY);
               ii = J(jj);
               k = ( det(Ai([1:jj-1,jj+1:As], [1:ii-1,ii+1:As]))*(-1)^(ii+jj) ) / Aadj(ii,jj);
               Aadj = k*Aadj;
            else				% rank < As-1
               Aadj = zeros(As);			% Aadj is zero
            end;
         end;  	
	     	  
         Aadjc(:,:,1+i) = Aadj;
         Aadjc(:,:,degofdet+2-i) = conj(Aadj); 
      end;	% for 

      % compute Aadj(1):
      A1 = Y(:,:,1);
      detA1 = det(A1);		% = det(A(1))
      detA(1) = detA1;

      [U1,S1,V1] = svd(A1);
      condA = S1(1,1)/S1(As,As);

      if abs(condA) < CondBnd,		% well conditioned - A(si)
         Aadj = (detA1*V1/S1)*U1';		% = inv(A(si)) * det(A(si)) = Aadj(si)
      else					% ill conditioned or singular
         anywarn = 1;
         tl = As * max(S1(:)) * eps; 
         if sum(diag(S1)>tl) >= As-1,	% rank is close to As-1
            Aadj = V1(:,end) * U1(:,end)';  	% outter product of right and left null spaces    
            [YY,J] = max(abs(Aadj));
            [Amax,jj] = max(YY);
            ii = J(jj);
            k = ( det(A1([1:jj-1,jj+1:As], [1:ii-1,ii+1:As]))*(-1)^(ii+jj) ) / Aadj(ii,jj);
            Aadj = k*Aadj;
         else				% rank < As-1
            Aadj = zeros(As);			% Aadj is zero
         end;
      end;  	
   
      Aadjc(:,:,1) = Aadj;

      % even number of points - compute det(Acmp(-1)):
      if rem(degofdet,2),
         i = fq2+2;		
         Ai = Y(:,:,i);
         detAi = det(Ai);
         detA(i) = detAi;		% = det(A(-1))

         [Ui,Si,Vi] = svd(Ai);
         condA = Si(1,1)/Si(As,As);
	
         if abs(condA) < CondBnd,		% well conditioned - A(si)
            Aadj = (detAi*Vi/Si)*Ui';		% = inv(A(si)) * det(A(si)) = Aadj(si)
         else				% ill conditioned or singular
            anywarn = 1;
            tl = As * max(Si(:)) * eps;
            if sum(diag(Si)>tl) >= As-1,	% rank is close to As-1
               Aadj = Vi(:,end) * Ui(:,end)';  	% outter product of right and left null spaces    
               [YY,J] = max(abs(Aadj));
               [Amax,jj] = max(YY);
               ii = J(jj);
               k = ( det(Ai([1:jj-1,jj+1:As], [1:ii-1,ii+1:As]))*(-1)^(ii+jj) ) / Aadj(ii,jj);
               Aadj = k*Aadj;
            else				% rank < As-1
               Aadj = zeros(As);			% Aadj is zero
            end;
         end;  	
    	    
         Aadjc(:,:,i) = Aadj;
      end;
	  
      warning(saved_w);

      % recover det(A(s)) and Aadj(s) (inverse FFT):
      if degofdet>=1,
         detA  = real( fft( conj(detA) ) ) ./ (degofdet+1);	
         Aadjc = real( fft( conj(Aadjc), [], 3 ) ) ./ (degofdet+1);	
      end;
 	  	
   else			% COMPLEX POLYNOMIAL matrix
      % compute det( A(si) ):
      for i = 1:degofdet+1,		% all Fourier points
         Ai = Y(:,:,i);
         detAi = det(Ai);		% = det( (A,si) )
         detA(i) = detAi;

         [Ui,Si,Vi] = svd(Ai);
         condA = Si(1,1)/Si(As,As);
	    
         if abs(condA) < CondBnd,		% well conditioned - A(si)
            Aadj = (detAi*Vi/Si)*Ui';	% = inv(A(si)) * det(A(si)) = Aadj(si)
         else			% ill conditioned or singular
            anywarn = 1;
            tl = As * max(Si(:)) * eps;
            if sum(diag(Si)>tl) >= As-1,	% rank is close to As-1
               Aadj = Vi(:,end) * Ui(:,end)';  	% outter product of right and left null spaces    
               [YY,J] = max(abs(Aadj));
               [Amax,jj] = max(YY);
               ii = J(jj);
               k = ( det(Ai([1:jj-1,jj+1:As], [1:ii-1,ii+1:As]))*(-1)^(ii+jj) ) / Aadj(ii,jj);
               Aadj = k*Aadj;
            else				% rank < As-1
               Aadj = zeros(As);			% Aadj is zero
            end;
         end;  	
    	    
         Aadjc(:,:,i) = Aadj;
      end;
      
      warning(saved_w);
      
      % recover det(A(s)) and Aadj(s) (inverse FFT):
      if degofdet>=1,
         detA  = conj( fft( conj(detA) ) ) ./ (degofdet+1);	
         Aadjc = conj( fft( conj(Aadjc), [], 3 ) ) ./ (degofdet+1);
      end;
   end;  	% if isreal ...
	
   if anywarn & all(As>1),
      warning(sprintf(['pol/adj\n',...
         '         Some interpolated matrices are ill conditioned.\n',...
         '	 Result may be inaccurate. Try the ''def'' method.']));
   end;	  		  
% ====================================
%    Definition-like computation
% ====================================      	
elseif strcmp(met, 'def'),
   if nargout == 2,

      % Compute DETERMINANT of A:
      if isreal(Ac),		% quicker for REAL POLYNOMIAL matrix
         % compute det( A(si) ):
         for i = 1:fq2,			% upper half of Fourier points
            detAi = det(Y(:,:,1+i));		
            detA(1+i) = detAi;
            detA(degofdet+2-i) = conj(detAi);
         end;

         % compute det(A(1)):
         detA(1) = det(Y(:,:,1)); 		% = det( pval(A,1) )

         % even number of points - compute det(A(-1)):
         if rem(degofdet,2),		
            detA(fq2+2) = det(Y(:,:,fq2+2));	% = det( (A, -1) )
         end;
   	    
         detA  = real( fft( conj(detA) ) ) ./ (degofdet+1);	
         
      else			% COMPLEX POLYNOMIAL matrix
         % compute det( A(si) ):
         for i = 1:degofdet+1,			% all Fourier points
            detAi = det(Y(:,:,i));	  	% = det( pval(A,si)*si^dif )
            detA(i) = detAi;
         end;

         % recover det(A(s)) (inverse FFT):
         if degofdet>=1,
            detA  = conj( fft( conj(detA) ) ) ./ (degofdet+1);	
         end;
      end;
   end; 	% if nargout..   

   % Compute ADJOINT af A:
   Aadjc = zeros([Asc, degofdet+1]);
   Ypom  = zeros([As-1 As-1 degofdet+1]);
   detAc = zeros(1, degofdet+1);

   if isreal(Ac),			% quicker for REAL POLYNOMIAL matrix	
      for ii = 1:As,
         for jj = 1:As,			
            Ypom(:,:,:) = Y([1:ii-1,ii+1:As],[1:jj-1,jj+1:As],:);	
				
            % compute det(Aij(s)) at Fourier points:
            for i = 1:fq2,			% upper half of Fourier points
               detAi = det(Ypom(:,:,1+i));		% = det( pval(A,si)*si^dif ), si = vec(i+1)
               detAc(1+i) = detAi;
               detAc(degofdet+2-i) = conj(detAi);
            end;

            % compute det(Aij(1)):
            detAc(1) = det(Ypom(:,:,1)); 		% = det( pval(A,1) )

            % even number of points - compute det(Acmp(-1)):
            if rem(degofdet,2),		
               detAc(fq2+2) = det(Ypom(:,:,fq2+2));	% = det( pval(A, -1)*(-1)^dif );
            end;

            Aadjc(jj,ii,:) = ((-1)^(ii+jj)) * detAc;

         end;
      end;	 
 	 
      % recover det(A(s)) and Aadj(s) (inverse FFT):
      if degofdet>=1,
         Aadjc = real( fft( conj(Aadjc), [], 3 ) ) ./ (degofdet+1); 	  	
      end;
	  
   else				% COMPLEX POLYNOMIAL matrix
      for ii=1:As,
         for jj = 1:As,			
            Ypom(:,:,:) = Y([1:ii-1,ii+1:As],[1:jj-1,jj+1:As],:);
            
            % compute det( A(si)*si^dif ):
            for i = 1:degofdet+1,			% all Fourier points
               detAi = det(Ypom(:,:,i));	  	% = det( pval(A,si)*si^dif )
               detAc(i) = detAi;
            end;
	
            Aadjc(jj,ii,:) = ((-1)^(ii+jj)) * detAc;
         end;
      end;

      % recover det(A(s)) and Aadj(s) (inverse FFT):
      if degofdet>=1,
         Aadjc = conj( fft( conj(Aadjc), [], 3 ) ) ./ (degofdet+1);
      end;
	  
   end; 	% if isreal

else error('Invalid command option.');
% =========================
end;
% =========================              	      	

% Aadj:
clear Aadj;
Aadj.d = degofdet;
Aadj.s = Asc;
Aadj.c = Aadjc;	
Aadj.v = A.v;
Aadj.h = A.h;
Aadj.u = [];
Aadj.version = 3.0;

% packaging:
Aadj = class(Aadj,'pol'); 

% det(A): 
degofdet = min( sum(deg(A, 'col')), sum(deg(A,'row')) );
dtA.d = degofdet;
degofdet(isinf(degofdet)) = 0;
dtA.s = [1 1];
dtAc = permute(detA, [3 1 2]);
dtA.c = dtAc(1:degofdet+1);
dtA.v = A.v;
dtA.h = A.h;
dtA.u = [];
dtA.version = 3.0;

% packaging:
dtA = class(dtA,'pol');  
 
% zeroing:
if ~isempty(me) & A.d > 0,
   Aadj = pzer(Aadj, tol * me/As);	% =?= pzer(Aadj, tol * me^(As-1));
   dtA  = pzer(dtA, tol * me/As);	% =?= pzer(dtA,  tol * me^As);
end;

%end .. @pol/adj
