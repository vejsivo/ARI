function d = det2d(varargin);
%DET2D  Determinant of a 2-D polynomial matrix
%
% Given the square polynomial matrices  P0, P1, ... , PN  of the
% same size, the command
%    DC = DET2D(P0, P1, ... , PN)
% computes a column polynomial vector
%    DC = [D0(VAR1); D1(VAR1); ... ; DM(VAR1)],
% such that the 2-D scalar polynomial
%    D(VAR1,VAR2) = D0(VAR1) + ... + DM(VAR1) * VAR2^M
% is the determinant of the 2-D polynomial matrix
%    P(VAR1,VAR2) = P0(VAR1) + ... + PN(VAR1) * VAR2^N.
% VAR1 is the variable of the input polynomial matrices P0, ... ,PN. 
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also POL/DET.

%       Author(s): M. Hromcik, M. Sebek 22-10-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 1.0 $  $Date: 29-Oct-1998 10:28:34   $
%       $Revision: 2.0 $  $Date: 10-Oct-2000 10:30:00 - Version 3.0, M. Hromcik $

global PGLOBAL;
eval('tol = PGLOBAL.ZEROING;', 'pinit;');

% ...............
% Process inputs:
% ...............

ni = nargin;

if ni == 0,
  error('Not enough input arguments.');

elseif ni == 1
  eval('d = det( pol(varargin{ni}) );', 'error(peel(lasterr));');
  return;

else,	% ni >= 2
  endin = varargin{ni};
  if isnumeric(endin) & all(size(endin)==1),
    tol = endin;
    %eval 'PP = [varargin{1:ni-1}];', ...
  	% 'error(''Inconsistent dimensions of the input matrices'')');
    try
        PP = [varargin{1:ni-1}];
    catch
        error('Inconsistent dimensions of the input matrices.');
    end;
    dr = ni-2;			% dr == degree of P in the 2nd variable VAR2.
  else
    dr = ni-1;
    %eval( 'PP = [varargin{:}];', ...
  	% 'error(''Inconsistent dimensions of the input matrices'')');
    try
        PP = [varargin{:}];
    catch
        error('Inconsistent dimensions of the input matrices.');
    end;
  end; 
end;  

if isempty(PP),
  d = pol(1);
  return;
end;
  
P = pol(PP);    
Pc = P.c;
Ps = P.s;
if Ps(1) ~= Ps(2)/(dr+1), error('All matrices must be square matrices of the same size.'); end;

% ........................................
% Evaluate P(r,s) at Fourier points ri,sj:
% ........................................

siz = size(P,1);
ds = P.d;			% ds == degree of P in the 1st variable VAR1.
degdet_s = siz*ds;		% upper bound for the detrminant degree in VAR1
if isinf(degdet_s), 
  d = pol(0);
  return;
end;
  
degdet_r = siz*dr;		% upper bound for the detrminant degree in VAR2
  
X = fft(Pc, degdet_s+1, 3);	% degree of det(P) in VAR1 <= size(P)*deg(P,VAR1)

X = reshape(X, [siz, siz, (degdet_s+1)*(dr+1)]);
X2 = zeros([siz, siz*(degdet_s+1), dr+1]);
for i=1:(degdet_s+1),
   X2(:, ((i-1)*siz+1): i*siz, :) = X(:,:, ((i-1)*(dr+1)+1):i*(dr+1));
end;  

Y = fft(X2,degdet_r+1,3);	% degree of det(P) in VAR2 <= size(P)*deg(P,VAR2)	

% .............................................
% Interpolate through determinants of P(ri,sj):
% .............................................

dt = zeros((degdet_s+1)*(degdet_r+1),1);
Y = reshape(Y, [siz,siz,(degdet_s+1)*(degdet_r+1)]);
for i=1:(degdet_s+1)*(degdet_r+1),
  dt(i) = det(Y(:,:,i));
end;
dt = reshape(dt, [degdet_s+1, degdet_r+1]);

% ...................................
% Backward transform of determinants:
% ...................................

if isreal(Pc),
  d = real(fftn(conj(dt))) / ((degdet_s+1)*(degdet_r+1));
else
  d = conj(fftn(conj(dt))) / ((degdet_s+1)*(degdet_r+1));
end;    

% Pre-zeroing:
d = d.';
me = min(abs(nonzeros(Pc)));
d( abs(d) < tol*me ) = 0;
[ds1,ds2] = size(d);
d = d( 1:max( [find(any(d,2))] ) , 1:max( [find(any(d,1))] ) );

% POL object:
if isempty(d), d = pol(0);
else d = pol(d, size(d,2)-1); end;

%end .. det2d
