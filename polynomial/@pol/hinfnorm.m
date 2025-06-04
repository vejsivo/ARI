function [Hmax, w_opt] = hinfnorm(N,D,arg3,arg4);
%HINFNORM  H-infinity norm of a stable polynomial matrix fraction
% 
% The commands
%    HINFNORM(N,D)
%    HINFNORM(N,D,'l') 
% compute the H-infinity norm of the stable rational matrix G = D\N.
% The command
%    HINFNORM(N,D,'r') 
% computes the H-infinity norm of the stable rational matrix G = N/D.
% The function call
%    [PEAK, OMEGAMAX] = HINFNORM(N,D), 
% possibly combined with the 'l' or 'r' option, additionally returns 
% the frequency OMEGAMAX where the 2-norm of G(J*OMEGA) achieves its
% maximum.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is PGLOBAL.ZEROING^2. The tolerance is used to
% test whether the denominator has roots on the imaginary axis (for
% the 's' and 'p' variables) or on the unit circle (for the other 
% variables).
%
%See also POL/H2NORM, POL/NORM.

%       Author(s): M. Hromcik, M. Sebek 22-10-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 1.0 $  $Date: 22-Oct-1998 10:28:34   $
%		  $ Version 3 - 22-Nov-2000, Martin Hromcik $
%       $Revision $  $Date: 02-Jan-2002  J.Jezek  $

global PGLOBAL;
opt = 'l';		% D^-1 N
tol = PGLOBAL.ZEROING^2; % tolerance for |poles| == 0

pts = 50;		% # points for 's'

ex  = 10;		% 2^ex points for 'z','d'
stp = 4;		% 2^(ex-stp) points in the 1st run, 2^stp in the 2nd run
			% in case of MIMO discrete

isscalar = 0;

switch nargin,
  case 1,
    D = 1;   
  case 2,
  case 3,
    if ischar(arg3),
      opt = arg3;
    elseif isnumeric(arg3),
      tol = arg3;
    else error('Invalid command option.')       
    end;       
  case 4,
    if ischar(arg3),
      opt = arg3;
    elseif isnumeric(arg3),
      tol = arg3;
    else error('Invalid command option.')       
    end;       
    if ischar(arg4),
      opt = arg4;
    elseif isnumeric(arg4),
      tol = arg4;
    else error('Invalid command option.')       
    end;      
  otherwise error('Too many input arguments.');
end;	% switch

if ~all(size(tol)==1) | ~isreal(tol) | ...
      tol<0 | tol>1,
   error('Invalid tolerance.');
end;
 
eval('N = pol(N);', 'error(peel(lasterr));');
eval('D = pol(D);', 'error(peel(lasterr));');

Ns = N.s; Ds = D.s;
j = sqrt(-1);

[tv,var,N,D] = testvp(N,D);
if tv==2,
   [th,h,N,D] = testhp(N,D,var);
   if ~th,
      warning('Inconsistent sampling periods.');
   end;
   if strcmp(D.v,'z^-1'),
      dg = deg(D);
      N = shift(N,dg); D = rev(D,dg);
      if deg(D)>0, D.v = 'z';
      else D.v = '';
      end;
   else
      dg = deg(N);
      N = rev(N,dg); D = shift(D,dg);
      if deg(N)>0, N.v = 'z';
      else N.v = '';
      end;
   end;
   var = 'z';
elseif ~tv,
   warning('Inconsistent variables.');
end;

[th,h,N,D] = testhp(N,D,var);
if ~th, warning('Inconsistent sampling periods.');
end;

% Size consistency
if Ds(1)~=Ds(2),
   error('Denominator matrix is not square.');
end;

switch opt,
case 'l',	% D\N
   if Ns(1)~=Ds(2),
      if Ds(2)==1,
         D = diag(repmat(D,Ns(1)));
         Ds(1) = Ns(1); Ds(2) = Ns(1);
      else
         error('Matrices of inconsistent dimensions.'); 
      end;
   end;
case 'r',	% N/D
   if Ns(2)~=Ds(1),
      if Ds(1)==1,
         D = diag(repmat(D,Ns(2)));
         Ds(1) = Ns(2); Ds(2) = Ns(2);
      else
         error('Matrices of inconsistent dimensions.'); 
      end;   
   end;
otherwise error('Invalid command option.');
end;	

if Ns(1) == 0 | Ns(2) == 0,
   Hmax = 0; w_opt = []; return;
end;
if all([Ns Ds] == 1), isscalar = 1;
end;
Nc = N.c; Dc = D.c;
Nd = N.d; Dd = D.d;
Nv = N.v; Dv = D.v;

if isscalar,	
  if any(strcmp(var,{'s','p'})),
	% ..........................
	% SISO case continuous time:
	% ..........................
	
	if D == 0, 
	  warning('Denominator is zero to working precision.');
	  Hmax = Inf;
	  w_opt = 0;
	  return;
	end;
	  
	% Try w = Inf:
	if Nd > Dd,
	  Hmax = Inf;
	  w_opt = Inf;
	  warning('Fraction is nonproper. H-inf norm is infinite.');
 	  return;
   elseif Nd == Dd,
     Nc = N.c; 
     Dc = D.c;
 	  H_at_inf = abs(Nc(:,:,Nd+1)/Dc(:,:,Nd+1));   
	else
	  H_at_inf = 0;
   end;  
   	
	rootsD = roots(D);
	if any(abs(real(rootsD)) < tol), 
	  Hmax = Inf;
	  rootsD = imag(rootsD( abs(real(roots(D))) < tol));	% select imaginary roots
	  w_opt = min(abs(rootsD));
	  if isempty(w_opt), w_opt = NaN; end;
	  warning(sprintf(['Denominator has a purely imaginary root.\n',...
	           '	 H-inf norm is infinite.']));
	  return;
	end;  

	if any(real(rootsD)>0), 
	  warning('Denominator is unstable.');
	end;  

	w_max = 10*max(abs([rootsD; roots(N)]));
	if ~isempty(w_max),
	  delta = w_max/pts;
	  w = (0:delta:w_max)*j;
	else
	  delta = 1;
	  w = [0 1 2];
	end;

	pts = length(w);

   % Eval N, D and N/D at all omega's:
	if isempty(N.c), Hmax = 0; w_opt = 0; return; 
   else
       Nc = N.c;
       Nw = Nc(:,:,Nd+1);
	    for i=Nd:-1:1,
	      Nw = Nw.*w + Nc(:,:,i);
	    end;
	end;    
 	if isempty(D.c), H = Inf; return; 
	else
      Dc = D.c;
      Dw = Dc(:,:,Dd+1);
	    for i=Dd:-1:1,
	      Dw = Dw.*w + Dc(:,:,i);
	    end;
	end; 	     
	H = abs(Nw./Dw);
	[Hmax, ind] = max(H);
	w_opt = w(ind);
	  
	% Refine - search in <w(ind-1), w(ind+1)>
	delta2 = 2*delta/pts;
	if ind==1, 
	  w = j*(0:delta2:delta);
	elseif ind==pts,
	  w = j*(imag(w(ind-1)):delta2:imag(w(ind)));
	else 
	  w = ( imag(w(ind-1)):delta2:imag(w(ind+1)) )*j;
	end;  
	pts = length(w);
	Nw = Nc(Nd+1);
	for i=Nd-1:-1:0,
	    Nw = Nw.*w + Nc(:,:,i+1);
	end;
	Dw = Dc(Dd+1);
	for i=Dd-1:-1:0,
	    Dw = Dw.*w + Dc(:,:,i+1);
	end; 	     
	H = abs(Nw./Dw);
	[Hmax2, ind2] = max(H);
	w_opt2 = w(ind2); 	 
   	if Hmax2>Hmax,
   	   w_opt = w_opt2;
   	   Hmax = Hmax2;
   	end;
   	w_opt = imag(w_opt);
   	
   	if H_at_inf > Hmax, 
   	  Hmax = H_at_inf;
   	  w_opt = Inf;
   	end;  

  else	%discrete time variables   
	% ........................
	% SISO case discrete time:
	% ........................

	if D == 0, 
	  warning('Denominator is zero to working precision.');
	  Hmax = Inf;
	  w_opt = 0;
	  return;
	end;
	
	% Test if proper:
  	%if ~isproper(N,D,opt,tol)
	%  warning('Nonproper matrix fraction - the infinity-norm is infinite.');
	%  Hmax = Inf;
	%  w_opt = Inf;
	%  return;
	%end;   	  

	% One-step direct search using FFT
	rootsD = roots(D);
	if any(abs(abs(rootsD)-1) < tol),
	  warning(sprintf(['Denominator has a root on the unit circle.\n',... 
	  '	 H-inf norm is infinite.']));    	
	  Hmax = Inf;
	  rootsD = rootsD( abs(abs(rootsD)-1) < tol );
	  arg = angle(rootsD);
	  w_opt = abs(max(arg)); 
	else
	  if (any(strcmp(var,{'z','q'})) & any(abs(rootsD)>1)) ...
	     | (any(strcmp(var,{'z^-1','d'})) & any(abs(rootsD)<1)),
	    warning('Denominator is unstable.');
	  end;  

	  if N.d <= 0, Nw = repmat(Nc, [1 1 2^ex]);
	  else, Nw = fft(Nc, 2^ex, 3); end;
          if D.d <= 0, Dw = repmat(Dc, [1 1 2^ex]);
          else, Dw = fft(Dc, 2^ex, 3); end;
          Hw = abs( Nw ./ Dw );
          [Hmax, ind] = max(Hw);
          w_opt = 2*pi*(ind-1)/(2^ex);
        end;
  end;	% if strcmp        
	
else,
  if any(strcmp(var,{'s','p'})),
  	% ..........................
  	% MIMO case continuous time:
  	% ..........................

	if tol == PGLOBAL.ZEROING, 
	  rankD = rank(D);
	else,
	  rankD = rank(D, tol);
	end;    

	if rankD < Ds(1),
	  warning('Denominator matrix is singular to working precision.');
	  Hmax = Inf;
	  w_opt = 0;
	  return;
	end;
	  
	rootsD = roots(D);
	if any(abs(real(rootsD)) < tol), 
	  Hmax = Inf;
	  rootsD = imag(rootsD( abs(real(roots(D))) < tol));	% select imaginary roots
	  w_opt = min(abs(rootsD));
	  if isempty(w_opt), w_opt = NaN; end;
	  warning(sprintf(['Denominator matrix has a purely imaginary root.\n',...
	  '	 H-inf norm is infinite.']));
	  return;
	end;  

	if any(real(rootsD)>0), 
	  warning('Denominator matrix is unstable.');
	end;  

	w_max = 10*max(abs([rootsD; roots(N)]));
	if ~isempty(w_max),
	  delta = w_max/pts;
	  w = (0:delta:w_max)*j;
	else
	  delta = 1;
	  w = [0 1 2];
	end;

	pts = length(w);

     % Eval N, D and D\N (or N/D) at all omega's:
	  H_at_inf = -1;
	  switch opt,
	     case 'l',
	      % Try infinity:
         [Nr,Dr] = reverse(N,D,'l');
         Drc = Dr.c;
         Nrc = Nr.c;
         if rank(Drc(:,:,1)) < min(Dr.s),
  	 	     Hmax = Inf;
  	 	     w_opt = Inf;
   	  	  warning('Fraction is nonproper. H-inf norm is infinite.');
   		  return
   		else
   		  H_at_inf = max(svd(Drc(:,:,1)\Nrc(:,:,1)));
	    	end;  

   		Nw = polyval(N,w);
	   	Dw = polyval(D,w);
      	for i=1:pts,
   		   s = svd(Dw(:,:,i) \ Nw(:,:,i));
   		   H(i) = s(1);
   		end;

        case 'r',
     	  	% Try infinity:
         [Nr,Dr] = reverse(N,D,'r');
         Drc = Dr.c;
         Nrc = Nr.c;
         if rank(Drc(:,:,1)) < min(Dr.s),
   		  Hmax = Inf;
   		  w_opt = Inf;
   	  	  warning('Fraction is nonproper. H-inf norm is infinite.');
   		  return
   		else
   		  H_at_inf = max(svd(Nrc(:,:,1)/Drc(:,:,1)));
   		end;  
	
     		Nw = polyval(N,w);
		   Dw = polyval(D,w);
       	for i=1:pts,
   		   s = svd(Nw(:,:,i) / Dw(:,:,i));
   		   H(i) = s(1);
       	end;

   otherwise error('Invalid command option.');
  	end;	%switch ...     
     
   [Hmax, ind] = max(H);
  	w_opt = w(ind);

  	% Refine - search in < w(ind-1), w(ind+1) >
  	delta2 = 2*delta/pts;
  	if ind==1, 
  	  w = j*(0:delta2:delta);
 	elseif ind==pts,
  	  w = j*(imag(w(ind-1)):delta2:imag(w(ind)));
  	else 
  	  w = ( imag(w(ind-1)):delta2:imag(w(ind+1)) )*j;
  	end;  
  	pts = length(w);
  	Nw = polyval(N,w);
  	Dw = polyval(D,w);
  	switch opt,
   	  case 'l',
   	    	for i=1:pts,
   		   s = svd(Dw(:,:,i) \ Nw(:,:,i));
   		   H(i) = s(1);
   		end;
     	  case 'r',	
      	 	for i=1:pts,
   		   s = svd(Nw(:,:,i) / Dw(:,:,i));
   		   H(i) = s(1);
        	end;
    	  otherwise error('Invalid command option.');
  	end;	%switch ...     
  	[Hmax2, ind2] = max(H);
  	w_opt2 = w(ind2); 	 
  	if Hmax2>Hmax,
  	   w_opt = w_opt2;
  	   Hmax = Hmax2;
  	end;
  	w_opt = imag(w_opt);
  	if H_at_inf > Hmax,
    	  Hmax = H_at_inf;
    	  w_opt = Inf;
  	end;  

   else		 	
   	% ........................
  	% MIMO case discrete time:
  	% ........................
  	
  	% Test if proper:
  	%if ~isproper(N,D,opt,tol)
	%  warning('Nonproper matrix fraction - the infinity-norm is infinite.');
	%  Hmax = Inf;
	%  w_opt = Inf;
	%  return;
	%end;   

	% Test D rank:
  	if tol == PGLOBAL.ZEROING, rankD = rank(D);
  	else, rankD = rank(D,tol); end;
  	if rankD < Ds(1),
  	  warning('Denominator matrix is singular to working precision.');
  	  Hmax = Inf;
  	  w_opt = 0;
  	  return;
  	end;  
  	
  	rootsD = roots(D);
	if any(abs(abs(rootsD)-1) < tol),
	  warning(sprintf(['Denominator matrix has a root on the unit circle.\n',... 
	  '	 H-inf norm is infinite.']));    	
	  Hmax = Inf;
	  rootsD = rootsD( abs(abs(rootsD)-1) < tol );
	  arg = angle(rootsD);
	  w_opt = abs(max(arg));
	  return;
	end;

	if (any(strcmp(var,{'z','q'})) & any(abs(rootsD)>1)) ...
	     | (any(strcmp(var,{'z^-1','d'})) & any(abs(rootsD)<1)),
	   warning('Denominator matrix is unstable.');
	end;  
	  	
 	if N.d <= 0, Nw = repmat(Nc, [1 1 2^ex]);
   else, Nw = fft(Nc, 2^ex, 3);
   end;
  	if D.d <= 0, Dw = repmat(Dc, [1 1 2^ex]);
   else, Dw = fft(Dc, 2^ex, 3);
   end;

  	% 1st run - rough grid:	2^(ex-stp) points ( = 128)
  	sigmas = zeros(1, 2^(ex-stp));
  	ptr = 1;
	switch opt,
	  case 'l',
	    for i=1:2^stp:(2^ex)/2,	% lower half of Fourier points
	      Hi = Dw(:,:,i) \ Nw(:,:,i );
	      if any(isnan(Hi)) | any(isinf(Hi)), 
	        warning(sprintf('Denominator matrix has a root near the unit circle.')); 
 	        Hmax = Inf;
 	        w_opt = (2*pi/(2^ex))*(i-1);
 	        return;
 	      else
 	        sigmai = svd(Hi);
 	      end;    
 	      sigmas(ptr) = sigmai(1);
  	      ptr = ptr+1;
  	    end;
  	  case 'r',    
	    for i=1:2^stp:2^ex,
  	      Hi = Nw(:,:,i) / Dw(:,:,i);
	      if any(isnan(Hi)) | any(isinf(Hi)), 
	        warning(sprintf('Denominator matrix has a root near the unit circle.')); 
 	        Hmax = Inf;
 	        w_opt = (2*pi/(2^ex))*(i-1);
 	        return;
 	      else
 	        sigmai = svd(Hi);
 	      end;    
	      sigmas(ptr) = sigmai(1);
  	      ptr = ptr+1;
     	    end;
  	  otherwise error('Invalid command option.');
  	end;
  	
	[Hmax, ind1] = max(sigmas);
	
	% 2nd run 
	if ind1 == 1,
	  lo = 1; hi = 2^stp+1;
	elseif ind1 == 2^ex,
	  lo = 2^ex - 2^step; hi = 2^step;
	else  
	  lo = (ind1-2)*2^stp + 1; hi = ind1*2^stp + 1;
	end;
       
       	sigmas = zeros(1, 2*(2^stp-1) );
	switch opt,
	  case 'l',
	    for i=lo:1:hi,
  	      Hi = Dw(:,:,i) \ Nw(:,:,i);
   	      if any(isnan(Hi)) | any(isinf(Hi)), 
	        warning(sprintf('Denominator matrix has a root near the unit circle.')); 
 	        Hmax = Inf;
 	        w_opt = (2*pi/(2^ex))*(i-1);
 	        return;
 	      else
 	        sigmai = svd(Hi);
 	      end;    
	      sigmas(i) = sigmai(1);
  	    end;
  	  case 'r',    
	    for i=lo:1:hi,
  	      Hi = Nw(:,:,i) / Dw(:,:,i);
	      if any(isnan(Hi)) | any(isinf(Hi)), 
	        warning(sprintf('Denominator matrix has a root near the unit circle.')); 
 	        Hmax = Inf;
 	        w_opt = (2*pi/(2^ex))*(i-1);
 	        return;
 	      else
 	        sigmai = svd(Hi);
	      end;    
 	      sigmas(i) = sigmai(1);
  	    end;
  	  otherwise error('Invalid command option.');
  	end;
  	
	[Hmax2, ind2] = max(sigmas);	    
	w_opt = 2*pi*(ind2-1)/(2^ex);
	
    end;	% if any(str...    

end;	%if isscalar ...

%end .. @pol/hinfnorm
