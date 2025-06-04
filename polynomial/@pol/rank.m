function r = rank(A, arg2, arg3);
%RANK   Rank of polynomial
%
% The commands
%    RANK(A) 
%    RANK(A,'FFT') 
% provide an estimate of the number of linearly independent rows 
% or columns of the polynomial matrix A by computing the rank of A 
% evaluated at a suitable number of Fourier points.
%
% The number
%    RANK(A, TOL) or RANK(A, 'FFT', TOL) 
% is the minimal number of singular values of A evaluated at a suitable
% number of Fourier points that are larger than TOL. The commands RANK(A) 
% and RANK(A,'FFT') use the default tolerance TOL = MAX(SIZE(A))*NORM(A)*EPS).
%
% The command
%    RANK(A,'SYL') or RANK(A,'SYL',TOL) 
% estimates the rank of the polynomial matrix A by computing ranks of 
% appropriate Sylvester matrices. The algorithm is based on evaluation 
% of the rank of a constant matrices and therefore the optional input 
% argument TOL has the same meaning as for the standard MATLAB function 
% RANK.

%       Author(s): D. Henrion, M. Hromcik, M. Sebek, S. Pejchova 20-5-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 20-May-1998 10:28:34   $
%       $Revision: 3.0 $  $Date: 11-Aug-1999 12:00:00   J.Jezek  $
%                         $Date: 09-Nov-1999 10:30:00   J.Jezek  $
%                         $Date: 28-Feb-2003            J.Jezek  $

eval('A = pol(A);','error(peel(lasterr));');
met = ''; tol = [];
na = nargin;
if na>=2,
   if isa(arg2,'double') & (isempty(arg2) | ...
         (length(arg2)==1 & isreal(arg2) & arg2>=0)),
      tol = arg2;  % here, exceptionally, arg2<=1 is not tested.
                   % It caused troubles. J.Jezek
   elseif isa(arg2,'char') & ...
         (isempty(arg2) | strcmp(arg2,'fft') | strcmp(arg2,'syl')),
      met = arg2;
   else
      error('Invalid 2nd argument.');
   end;
end;
if na==3,
   if isa(arg3,'double') & (isempty(arg3) | ...
         (length(arg3)==1 & isreal(arg3) & arg3>=0)),
      tol = arg3;  % here, exceptionally, arg3<=1 is not tested.
                   % It caused troubles. J.Jezek
   elseif isa(arg3,'char') & ...
         (isempty(arg3) | strcmp(arg3,'fft') | strcmp(arg3,'syl')),
      met = arg3;
   else
      error('Invalid 3rd argument.');
   end;
end;

Default = isempty(tol);
if isempty(met), met = 'fft';
end;
Ac = A.c;
As = A.s;
minAs = min(As);
maxAs = max(As);
Ad = A.d;

if isempty(Ac),
  r = 0;
  return;
end;

if size(Ac, 3) == 1,
  Si = svd(Ac);
  if Default,
    tol = maxAs * max(Si(:)) * eps;
  end;
  r = sum(Si>tol);
  return;
end;

% ============ FFT interpolation approach ============== %
if strcmp(met, 'fft'),
	% expected #points:
	degofdet = min(As)*Ad;
	fq2 = floor(degofdet/2);

	% FFT of polynomial matrix A(s)- perform fast radix-2 algorithm:
		  
	Y = fft( Ac, 2^(ceil(log2(degofdet+1))), 3 ); 	  
	
	%eval('Y = fft( Ac, 2^(ceil(log2(degofdet+1))), 3 );', ...
	%     'Y = Ac; degofdet = 0; fq2 = -1;');
	
	if isreal(Ac),
	  R = zeros(1,fq2+2);
	  for i = 1:fq2+2,
	  	Si = svd(Y(:,:,i));
		if Default,
	    	 tol = maxAs * max(Si(:)) * eps;
		end;
	  	ri = sum(Si>tol);
		if ri == minAs,
		  r = minAs;
		  return;
		end;
	  	R(i) = ri;
	  end;
 	  	
	else				% COMPLEX POLYNOMIAL matrix
	  R = zeros(1,degofdet+1);
	  for i = 1:degofdet+1,
	  	Si = svd(Y(:,:,i));
		if Default,
	    	 tol = maxAs * max(Si(:)) * eps;
		end;
	  	ri = sum(Si>tol);
		if ri == minAs,
		  r = minAs;
		  return;
		end;
	  	R(i) = ri;
	  end;
	end;

	r = max(R);

% ============ Sylvester resultant approach ============== %
elseif strcmp(met, 'syl'),

 	% The rank of an arbitrary polynomial matrix is evaluated as the
 	% difference between ranks of two Sylvester matrices.

 	% order of Sylvester matrices, cf. [1]
 	cdeg = deg(A, 'col'); cdeg(cdeg < 0) = 0; cdeg = -sort(-cdeg);
 	rdeg = deg(A, 'row'); rdeg(rdeg < 0) = 0; rdeg = -sort(-rdeg);
 	order = min(sum(cdeg(1:(A.s(2)-1))), sum(rdeg(1:(A.s(1)-1))));

 	if Default, % built-in tolerance

  	   rrow = rank(sylv(A, order));
	   rcol = rank(sylv(A.', order));
	   if order > 0, 
	    rrow = rrow - rank(sylv(A, order-1));
	    rcol = rcol - rank(sylv(A.', order-1));
	   end;

 	else % user-supplied tolerance

	  rrow = rank(sylv(A, order), tol);
	  rcol = rank(sylv(A.', order), tol);
	  if order > 0,
	    rrow = rrow - rank(sylv(A, order-1), tol);
	    rcol = rcol - rank(sylv(A.', order-1), tol);
	  end;

	end;
	
	r = min(rrow, rcol);

end;

%end .. @pol/rank
