function T = tsp(varargin)
%TSP    Create two-sided polynomial
%
% The commmand
%    T =  TSP(N,P) 
% creates a two-sided polynomial matrix object 
%    T(z) = NK*z^-k + ...N1*z^-1 + (N0+P0) + P1*z + P2*z^2 + .. + PM*z^m
% from the polynomial matrix N (of degree k) and polynomial matrix P (of degree m).
% The part consisting of negative degrees of 'z' is a 'negative part' of TSP object
% and similarly the 'positive part' of TSP contains the positive degrees of 'z'.
% The variable string 'z' ('q') or 'z^-1' ('d') of the input polynomial matrix determines
% if it is transformed to a 'positive' or 'negative' part of TSP object respectively.
% The input order determines transformation to a 'negative' or 'positive' part
% from polynomial matrices with variable strings 's' and 'p'.
% 
% The matrix coefficients of T are ordered according to ascending powers of 'z'.
%
% The commmand
%    T =  TSP(A) 
% returns T = A if A is already a two-sided polynomial matrix object. If A is a 
% (standard Matlab) constant matrix then it returns the corresponding 
% zero degree two sided polynomial matrix object. If A is a polynomial matrix
% with variable string 's' or 'p', T is a two-sided polynomial matrix containing
% only a 'positive' part.
%
% The variable string of TSP is always  'z', it is independent of variable strings
% of input polynomial objects.

% The internal structure of the two-sided polynomial matrix object T is as follows:
%
% T.p is a POL object created from inputs N and P and shifted to positive powers of 'z'.
%
% T.s is the size of T(z), a two-element row vector [nrow, ncol] containing
% the number of rows and columns;
%
% T.d is the degree of T(z), T.d = m . For zero tsp, T.d = -Inf .
%
% T.t is the trailing degree of T(z), T.t = -k . For zero tsp, T.t = +Inf .
%
% T.o is a offset of T.p or the degree of 'z' to obtain from T.p the original,
% if both input polynomial objects N and P are cleared then T.o=-k (degree of N).
%
% T.r is the range of T(z), it is the degree of the POL object T.p, if both
% input polynomial objects N and P are cleared then T.r = m+k.
%
% T.v is the variable string - always 'z'; for empty and constant TSP objects is T.v=''.
%
% T.h is the sampling period;
%
% T.u is a user data field;
%
% T.version is the version number, always 3.0;

%       Author(s):  S. Pejchova, M. Sebek 01-7-99
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 3.0 $Date: 11-Aug-1999 12:00:00  J.Jezek  $
%                      $Date: 27-Mar-2000 10:52:11   $
%                      $Date: 22-May-2000 10:30:00  J.Jezek  $
%                      $Date: 11-Jul-2000 13:40:00  S.Pejchova  $
%                      $Date: 24-Jul-2000  J.Jezek, S.Pejchova  $
%                      $Date: 31-Jul-2000  J.Jezek  $
%                      $Date: 15-Jul-2001  J.Jezek  $
%                      $Date: 24-Feb-2002  J.Jezek  $
%                      $Date: 28-Feb-2003  J.Jezek  $

global PGLOBAL;

na = nargin;
if na>=1 & isa(varargin{1},'tsp'),

   % Quick exit for TSP(T) with T of class TSP
   T = varargin{1};
   if na > 1,
      props(T,varargin{2:na});
%      error('Use PROPS to modify the properties of the TSP object.');
   end;
   return;
end;

if na>=3 & (isa(varargin{1},'double') & ... %state space arguments
      isa(varargin{2},'double') & isa(varargin{3},'double')),
   A = varargin{1}; B = varargin{2}; C = varargin{3};
   n = size(A,1); m = size(B,2); p = size(C,1);
   
   tol = [];
   lastarg = varargin{na};
   if isa(lastarg,'char'),
      if ~isempty(lastarg),
         switch lastarg,
         case {'z','zi','z^-1'},
         otherwise
            tol = str2double(lastarg);
            if isnan(tol),
               error('Invalid last argument.');
            end;
         end;
      end;
      na = na-1;
      lastarg = varargin{na};
      if isa(lastarg,'char'),
         if ~isempty(lastarg),
            switch lastarg,
            case {'z','zi','z^-1'},
            otherwise
               tol = str2double(lastarg);
               if isnan(tol),
                  error('Invalid last but one argument.');
               end;
            end;
         end;
         na = na-1;
      end;
   end;
   if isempty(tol), tol = PGLOBAL.ZEROING;
   end;
   
   savedvar = PGLOBAL.VARIABLE; PGLOBAL.VARIABLE = 'z';
   if na<=4,
      if na==3,
         D = zeros(p,m);
      elseif na==4,
         D = varargin{4};
         if ~isa(D,'double') & ~isa(D,'pol'),
            error('Invalid 4th argument.');
         end;
      end;
      N = 0; P = 0;
      if p>=m,
         eval('[N,D] = ss2rmf(A,B,C,D,tol);','error(peel(lasterr));');
         P = rdf(N,D);
      else
         eval('[N,D] = ss2lmf(A,B,C,D,tol);','error(peel(lasterr));');
         P = ldf(D,N);
      end;
   
   elseif na==5,
      D = varargin{4}; E = varargin{5};
      if ~isa(D,'double'),
         error('Invalid 4th argument.');
      end;
      if ~isa(E,'double'),
         error('Invalid 5th argument.');
      end;
      N = 0; P = 0;
      if p>=m,
         eval('[N,D] = dss2rmf(A,B,C,D,E,tol);','error(peel(lasterr));');
         P = rdf(N,D);
      else
         eval('[N,D] = dss2lmf(A,B,C,D,E,tol);','error(peel(lasterr));');
         P = ldf(D,N);
      end;
   else
      error('Too many input arguments.');
   end;
   PGLOBAL.VARIABLE = savedvar;
   
   P.v = 'z';
   eval('T = tsp(P);', 'error(peel(lasterr));');
   return;
end;
   
if na>=1 & isa(varargin{1},'sym'),    %symbolic argument
   lv = length(varargin);
   eval('T = tsp(sdf(varargin{1:lv}));', 'error(peel(lasterr));');
   return;
end;

superiorto('double','pol');

eval('PGLOBAL.VARIABLE;','painit;');   
  
N = pol([]); P = pol([]); nrow = 0; ncol = 0;
PropValStart = 0;   

if na==1 | ((na>1) & isa(varargin{2},'char')),
   eval('P = pol(varargin{1});', ...
        'error(''Argument is not convertible to tsp.'');');
   [nrow, ncol] = size(P);
   N=pol(zeros(nrow,ncol));
   if na>1, PropValStart = 2;
   end;   
elseif na>1,
   eval('N = pol(varargin{1}); P = pol(varargin{2});', ...
      'error(''Arguments are not convertible to tsp.'');');
   
   %Dimensions checking
   [td,P,N] = testdp(P,N);
   if td==0,
      error('Matrices not of the same dimensions.');
   end;
   [rP, cP] = size(P); [rN, cN] = size(N);
   nrow = rP; ncol = cP;
   if na>2, PropValStart = 3;
   end      
end;

Zmx = zeros(nrow,ncol);
Q = pol(Zmx);

% Default property values
T.p = Q;
T.s = [nrow, ncol];
T.d = [];
T.t = [];
T.o = 0;
T.r = [];
T.v = '';
T.h = [];
T.u = [];
T.version = 3.0;

wr1 = 0; 
switch  N.v,
case {'d','s','p'},
   props(N,'z^-1'); wr1=1;
case 'q',
   props(N,'z');  wr1=1;
end;

switch  P.v,
case {'q','s','p'},
   props(P,'z');  wr1=1;
case 'd',
   props(P,'z^-1'); wr1=1;
end;
switch length([N.v,P.v]),
case 5,
   if strcmp(P.v,'z^-1'), Q=P; P=N; N=Q;
   end;
   Nc=N.c; Pc=P.c; 
   P0=Pc(:,:,1)+Nc(:,:,1);
   Nc=Nc(:,:,end:-1:2); Pc=Pc(:,:,2:end);
   Nd=size(Nc,3);
   Qd=size(Pc,3)+Nd;
   Qc=[Nc(:,:),P0,Pc(:,:)];
   Q=pol(Qc,Qd,'z');
   [th,h,N,P] = testhp(N,P,'z');
   T.o=-Nd;
case {0,1,2},
   [th,h,N,P] = testhp(N,P,P.v); Q = P+N;
case {4,8},
   [th,h,N,P] = testhp(N,P,P.v); Q = P+N;
   if ~isempty(Q.v),
      Qc=Q.c; Qd=Q.d;
      Qc=Qc(:,:,end:-1:1);
      Q=pol(Qc(:,:),Qd,'z');
      T.o=-Qd;   
   end;
end;

if ~isempty(Zmx),
   T.p = Q;
   T.r = Q.d;
   if isinf(Q.d), 
      T.d = -inf; T.t = +inf;
   else,
      T.d = T.r+T.o; T.t = T.o;
   end;
   if (Q.d>0)|T.o,
      T.v = 'z';
   end;
end;

T = class(T,'tsp');
T = tclear(T);

% Set Property Values
if (PropValStart>0) & (PropValStart <= na),
   for ii = PropValStart:na,
      Pv = varargin{ii};
      if isa(Pv,'pol'),
         [vs1,vs2,vd]=size(Pv);
         if all([vs1,vs2,vd]==1)&(~any(Pv.c(:,:)-[0,1])),
            Pv = Pv.v;
         end; 
      end;      
      if isa(Pv,'double') & (isempty(Pv) | ...
            (length(Pv)==1 & isreal(Pv) & isfinite(Pv) & Pv>=0)),
         h = Pv; th = 1;
      elseif ~isa(Pv,'char')|(~strcmp(Pv,'z')&(~strcmp(Pv,'z^-1'))&(~strcmp(Pv,'zi'))),
         T.u = Pv; 
      else, 
         wr1 = 0;
      end; 
   end;
end;

if ~isempty(T.v)&wr1,
   warning('Variable symbol is changed to ''z''.');
end;

if th==0,
   warning('Inconsistent sampling periods.');
end;
T.h = h; T.p.h = h;

%end .. @tsp/tsp
