function [P,S] = polfit(X,Y,N,tol);
%POLFIT Fit polynomial matrix to data
%
% The command
%    P = POLFIT(X,Y,N) 
% finds the polynomial matrix P whose entries have degrees 
% less than or equal to N(I,J) that fits the data 
%    P(X(K)) ~= Y(:,:,K), 
% in a least-squares	sense. If the third argument is missing
% then the overall degree N is set equal to N = LENGTH(X)-1.
%
% The command
%    [P,S] = POLFIT(X,Y,N) 
% returns the polynomial matrix P and a structure array S for use 
% with POL/POLYVAL to obtain error estimates on predictions. See 
% POLYFIT for details.
%
% The commands
%    POLFIT(X,Y,N,TOL)
%    POLFIT(X,Y,[],TOL) 
% work with zeroing specified by the input tolerance TOL.
%
% The length of X and the third size of Y must agree. The first 
% and second sizes of Y must equal the size of N, unless N is 
% a scalar. 
%
% See also POLYFIT, POLYVAL, ROOTS.

%       Author(s): M. Hromcik, M. Sebek 3-9-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 7-Sep-1998 10:28:34   $
%       Modified by J. Jezek, Aug 2001, arg checking

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;'); 

ni = nargin;
no = nargout;
if ni<2, error('Not enough input arguments.');
end

if ni == 2, 
 N = length(X)-1;
 tol = PGLOBAL.ZEROING;
elseif ni == 3,
 if isempty(N),
  N = length(X)-1;    
 end;
 tol = PGLOBAL.ZEROING;
elseif ni == 4,
 if isempty(N),
  N = length(X)-1;  
 end;
 if ~isa(tol,'double') | length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
    error('Invalid tolerance.');
 end;
end;  
  
mn = max(N(:));
sn = size(N);
sx = size(X);
sy = size(Y);

if ~isa(X, 'double'), error('1st argument must be a matrix.');
end;
if ~isa(Y, 'double'), error('2nd argument must be a matrix.');
end;
if ~isa(N, 'double') | any(any(floor(N)~=N)) | any(any(N<0)),
   error('The degree matrix must be a nonnegative integer array.');
end;
if any(sx==0), error('1st argument is empty.');
end;
if any(sy==0), error('2nd argument is empty.');
end;
if ~any(sx==1), error('1st argument must be a vector.');
end;

if ~isequal(max(sx),size(Y,3)), error('Inconsistent sizes of the data sets.');
end;

if any([size(Y,1) size(Y,2)] - sn),
   if all(sn == 1),
      N = repmat(N, sy);
      sn = sy;
   elseif sy(1) == 1 & sy(2) == 1,
      Y = repmat(Y, sn);  
      sy = sn;
   else,
      error('The matrix of entry degrees must match the data array.');
   end;
end;

P = zeros( [sy(1:2), mn+1] );
	
if no <= 1,
   for i = 1:sy(1),
    for j = 1:sy(2),
     Yi = Y(i,j,:);
	  Pi = polyfit( X(:), Yi(:), N(i,j) );
	  M = max(abs(Pi));
	  Pi( abs(Pi) < M*tol) = 0;
	  P( i,j,mn - N(i,j)+1:mn+1 ) = Pi;
	 end;  
	end;
	
elseif no == 2,
	for i = 1:sy(1),
	 for j = 1:sy(2),
	  Yi = Y(i,j,:);
	  [Pi , S(i,j)] = polyfit( X(:), Yi(:), N(i,j) );
	  M = max(abs(Pi));
	  Pi( abs(Pi) < M*tol) = 0;
	  P( i,j,mn - N(i,j)+1:mn+1 ) = Pi;
	 end; 
	end; 
		 
else error('Too many output arguments.');
end;

P = flipdim(P, 3);
P = pol(P(:,:), mn);
%P = pzer(P, 0);

%end .. polfit
