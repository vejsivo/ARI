function [Ts,m] = shift(T,n)
%SHIFT     Shift two-sided polynomial
%                
% TS = SHIFT(T,N) with tsp matrix T and integer scalar N computes
%     TS(z) = T(z)* z^N .
%
% N can be also an integer matrix of the same dimensions as T
% unless T is scalar. The shift is performed elementwise:
% every T(i,j) is shifted by N(i,j). Scalar T can be shifted
% by any integer matrix N, the dimensions of the result being
% the same as those of N.
%
% PS = SHIFT(T) computes
%     PS(z) = T(z)* z^M
% with least possible scalar integer M such that the result PS
% is polynomial. M may be positive, zero or negative.
%
% [PS,M] = SHIFT(T)
% acts as above but returns also M.
%
% See also TSP/TIMES, TSP/MTIMES, TSP/POWER, TSP/MPOWER.

%        Author:  J. Jezek  11-8-99
%        Copyright (c) by Polyx, Ltd.
%        $Revision: 3.0 $ S. Pejchova  $Date: 17-Sep-1999 18:30:34   $
%                       $ J. Jezek     $Date: 13-Oct-1999 12:00:00   $
%                                      $Date: 22-May-2000 14:00:00   $

if nargin==2,
   if ~isa(n,'double') | ndims(n)>2 | ~isreal(n) | ...
         (~isempty(n) & any(any(floor(n)~=n)) ),
      error('Invalid shift; must be integer.');
   end;
   Tr = T.r;
   if length(n)==1,            % n is scalar
      Ts = T;
      if ~isempty(Tr) & isfinite(Tr),
         Ts.o = Ts.o + n;
         Ts.d = Ts.d + n;
         Ts.t = Ts.t + n;
      end;
      if (Tr>0)|Ts.o, Ts.v = 'z';
      else, Ts.v = '';
      end;
      
   else                        % n is matrix
      
      nsiz = size(n); Tsiz = T.s;
      if all(Tsiz==1),                 % T is scalar
         T = ones(nsiz)*T;
         Tsiz = nsiz;
      else                             % T is matrix
         if any(Tsiz~=nsiz),
            error('Matrices not of the same dimensions.');
         end;
      end;
      
      if ~isempty(Tr) & isfinite(Tr) & ~isempty(n),
         nmin = min(min(n)); nmax = max(max(n));
         nmaxmin = nmax - nmin;
         Tp = T.p; Tpc = Tp.c; Tpd = Tp.d;
         Tspc = zeros(Tsiz(1),Tsiz(2),Tpd+1+nmaxmin);
         for i = 1:Tsiz(1),
            for j = 1:Tsiz(2),
               nijn = n(i,j)-nmin+1;
               Tspc(i,j,nijn:nijn+Tpd) = Tpc(i,j,1:1+Tpd);
            end;
         end;
         
         Tspd = Tpd+nmaxmin;
         Tsp = pol(Tspc(:,:),Tspd,'z');
         
         Ts.p = Tsp;
         Ts.s = Tsiz;
         Ts.d = Tsp.d-T.o-nmin;
         Ts.t = T.o+nmin;
         Ts.o = T.o+nmin;
         Ts.r = Tsp.d;
         Ts.v = 'z';
         Ts.h = T.h;
         Ts.u = [];
         Ts.version = 3.0;
         Ts = class(Ts,'tsp');
      else
         Ts = T;
      end;
   end;
   Ts = tclear(Ts);   % end  nargin==2

else                  % nargin==1
   Ts = T.p;
   m = -T.o;
end;

%end .. @tsp/shift
 