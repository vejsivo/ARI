function M = mdf(varargin)
%MDF     Create matrix-denominator fraction
%
% The command  M = MDF(N,D) creates matrix-denominator fraction M from
% given polynomial matrices N and D. The (i,j)-th entry of M is  
% N(i,j)/D(i,j). All entries of D must be nonzero. Polynomials N,D
% must be finite (Inf, NaN not allowed).
%
% This syntax may be followed by PropertyValue arguments.
% See PROPS for details on assignable properties.
%
% Matrices N,D must have the same dimensions, unless one of them is
% scalar. In such a case, this scalar polynomial is divided by (or into)
% every entry of the second matrix.
%
% The command  M = MDF(N)  returns M = N if N is already a matrix-den
% fraction. Otherwise, it means  M = MDF(N,1).
%
% The variable symbols of polynomials N,D should be the same. When they
% are not, a warning is issued and the symbols are changed to
% PGLOBAL.VARIABLE. However, if one symbol is 'z' and the other 'z^-1'
% then the symbols play a role, no warning being issued, the resulting
% symbol is taken from N.
%
% A matrix-denominator fraction from polynomials N,D can be also
% created by operations  N./D  or  D.\N  .
%
% The command
%    M = MDF(A,B,C)  or  M = MDF(A,B,C,D)  or  M = MDF(A,B,C,D,E)
% creates matrix-den fraction M from state space matrices A,B,C and
% possibly D,E. In this case, an optional input argument TOL may
% specify the zeroing tolerance to be used instead the standard
% one; this argument must have a form of character string, e.g.
% MDF(A,B,C,'10^-6'). Also in this case, an optional input 
% argument VAR may specify the variable symbol to be used
% instead the standard one. It must have a form of character
% string; possible values are 's', 'p', 'z', 'zi', 'z^-1',
% 'q' or 'd'.
%
% The command
%    M = MDF(SYS)
% creates matrix-den fraction M from Control System Toolbox object
% SYS. Also in this case, the optional arguments TOL or VAR
% may be used.
%
% The command
%    M = MDF(SYMB)
% creates matrix-den fraction M from Symbolic Toolbox object SYMB.
% Optional argument VAR may be used.
%
% See also POL/LDIVIDE, POL/RDIVIDE.

% For the structure of object "matrix-denominator fraction",
% see FRAC/FRAC.

%       Author:  J. Jezek,  24-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $Revision: $  $Date: 21-Apr-2000 $
%                     $Date: 15-Jun-2000 $
%                     $Date: 31-Jul-2000  D.Henrion $
%                     $Date: 03-Aug-2000  J.Jezek   $
%                     $Date: 06-Nov-2000  J.Jezek   $
%                     $Date: 24-Jan-2002  J.Jezek   $
%                     $Date: 30-Sep-2002  J.Jezek   $
%                     $Date: 28-Feb-2003  J.Jezek   $
%                     $Date: 23-Nov-2009  M.Sebek   $ 
%                           fixed sampling period for "z" 
%                           for state-space input                        

global PGLOBAL;

ni = nargin;
if ni>=1 & isa(varargin{1},'mdf'),      % quick exit
   M = varargin{1};
   if ni==1 | (ni==2 & isa(varargin{2},'double')), %TOL ignored
      if strcmp(PGLOBAL.COPRIME,'cop'), M = coprime(M);
      end;
      if strcmp(PGLOBAL.REDUCE,'red'), M = reduce(M);
      else M = smreduce(M);
      end;
   else
      props(M,varargin{2:ni});
   end;
   return;
end;

eval('PGLOBAL.VARIABLE;','painit;');

if ni>=3 & (isa(varargin{1},'double') & ... %state space arguments
      isa(varargin{2},'double') & isa(varargin{3},'double')),
   A = varargin{1}; B = varargin{2}; C = varargin{3};
   n = size(A,1); m = size(B,2); p = size(C,1);

   var = PGLOBAL.VARIABLE; tol = []; proper = logical(0);
   lastarg = varargin{ni};
   if isa(lastarg,'char'),
      if ~isempty(lastarg),
         switch lastarg,
         case {'s','p','z','q','d','zi','z^-1'}
            var = lastarg;
         otherwise
            tol = str2double(lastarg);
            if isnan(tol),
               error('Invalid last argument.');
            end;
         end;
      end;
      ni = ni-1;
      lastarg = varargin{ni};
      if isa(lastarg,'char'),
         if ~isempty(lastarg),
            switch lastarg,
            case {'s','p','z','q','d','zi','z^-1'}
               var = lastarg;
            otherwise
               tol = str2double(lastarg);
               if isnan(tol),
                  error('Invalid last but one argument.');
               end;
            end;
         end;
         ni = ni-1;
      end;
   end;
   if isempty(tol), tol = PGLOBAL.ZEROING;
   end;
   
   savedvar = PGLOBAL.VARIABLE; PGLOBAL.VARIABLE = 'z';
   if ni<=4,
      proper = logical(1);
      if ni==3,
         D = zeros(p,m);
      elseif ni==4,
         D = varargin{4};
         if ~isa(D,'double') & ~isa(D,'pol'),
            error('Invalid 4th argument.');
         end;
      end;
      N = 0;
      if p>=m,
         eval('[N,D] = ss2rmf(A,B,C,D,tol);','error(peel(lasterr));');
         [AdjD,D] = adj(D,tol);
         N = mtimes(N,AdjD,tol);
      else
         eval('[N,D] = ss2lmf(A,B,C,D,tol);','error(peel(lasterr));');
         [AdjD,D] = adj(D,tol);
         N = mtimes(AdjD,N,tol);
      end;
      if ~all(N.s==D.s),
         D = repmat(D,N.s);
      end;
      
   elseif ni==5,
      D = varargin{4}; E = varargin{5};
      if ~isa(D,'double'),
         error('Invalid 4th argument.');
      end;
      if ~isa(E,'double'),
         error('Invalid 5th argument.');
      end;
      N = 0;
      if p>=m,
         eval('[N,D] = dss2rmf(A,B,C,D,E,tol);','error(peel(lasterr));');
         [AdjD,D] = adj(D,tol);
         N = mtimes(N,AdjD,tol);
      else
         eval('[N,D] = dss2lmf(A,B,C,D,E,tol);','error(peel(lasterr));');
         [AdjD,D] = adj(D,tol);
         N = mtimes(AdjD,N,tol);
      end;
      if ~all(N.s==D.s),
         D = repmat(D,N.s);
      end;
   else
      error('Too many input arguments.');
   end;
   PGLOBAL.VARIABLE = savedvar;

   F = frac(N,D);
   M.title = 'mdf';
   superiorto('double','pol','tsp');
   M = class(M,'mdf',F);

   switch var,
%   case {'s','p','q','z'}  % M.Sebek, Nov 2009
%      M.frac.v = var;
%      M.frac.h = 0;
   case {'s','p'}
   M.frac.v = var;
   M.frac.h = 0;
   case {'q','z'} % M.Sebek, Nov 2009: sampling period correctly set to default
   M.frac.v = var;
   case {'d'}
      M.frac.v = 'd'; M = reverse(M);
   case {'zi','z^-1'}
      M.frac.v = 'z'; M = reverse(M);
   end;
   
   if proper, props(M,'prop',tol);
   end;
   
   if strcmp(PGLOBAL.COPRIME,'cop'), M = coprime(M,tol);
   end;
   if strcmp(PGLOBAL.REDUCE,'red'), M = reduce(M,tol);
   else M = smreduce(M);
   end;
     
elseif ni>=1 & isa(varargin{1},'lti'),    % LTI argument
   sys = varargin{1};
   [p,m] = size(sys);
   var = ''; tol = [];
   if ni>=2,
      for i = 2:ni,
         arg = varargin{i};
         if isa(arg,'char'),
            switch arg,
            case {'s','p','z','q','d','zi','z^-1'}
               var = arg;
            otherwise
               tol = str2num(arg);
               if length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
                  error(['Invalid ',nth(i),' argument.']);
               end;               
            end;
         elseif isa(arg,'double'),
            tol = arg;
            if length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
               error(['Invalid ',nth(i),' argument.']);
            end;
         else
            error(['Invalid ',nth(i),' argument.']);            
         end;
      end;
   end;
   if isempty(tol), tol = PGLOBAL.ZEROING;
   end;
   
   if isa(sys,'ss'),      % SS argument
      if m==0 | p==0,
         N = pol(zeros(p,m)); D = N;
      else
         N = 0;
         if p>=m,
            eval('[N,D] = lti2rmf(sys,tol);','error(peel(lasterr));');
            [AdjD,D] = adj(D,tol);
            N = times(N,AdjD,tol);
         else
            eval('[N,D] = lti2lmf(sys,tol);','error(peel(lasterr));');
            [AdjD,D] = adj(D,tol);
            N = times(AdjD,N,tol);
         end;
      end;
            
   else                   % TF of ZPK argument
      if isa(sys,'zpk'), sys = tf(sys);
      end;
      num = sys.num; den = sys.den;
      N = pol(zeros(p,m)); D = N;
      if p>0 & m>0,
         for i = 1:p,
            for j = 1:m,
               numij = num{i,j}; denij = den{i,j};
               d = size(numij,2)-1;
               N(i,j) = lop(numij,d); D(i,j) = lop(denij,d);
            end;
         end;
         N.v = sys.var; D.v = sys.var;
         N.h = sys.ts; D.h = sys.ts;
      end;
   end;
   
   if ~all(N.s==D.s),
      if all(N.s==1),
         N = repmat(N,D.s);
      elseif all(D.s==1),
         D = repmat(D,N.s);
      end;
   end;

   F = frac(N,D);
   M.title = 'mdf';
   superiorto('double','pol','tsp');
   M = class(M,'mdf',F);

   switch var,
   case {'s','p','q','z'}
      M.frac.v = var;
   case {'d'}
      M.frac.v = 'd'; M = reverse(M);
   case {'zi','z^-1'}
      M.frac.v = 'z'; M = reverse(M);
   end;
   
   if strcmp(PGLOBAL.COPRIME,'cop'), M = coprime(M,tol);
   end;
   if strcmp(PGLOBAL.REDUCE,'red'), M = reduce(M,tol);
   else M = smreduce(M);
   end;
   
elseif ni==0,      % no arguments
   F = frac([],1);
   M.title = 'mdf';
   superiorto('double','pol','tsp');
   M = class(M,'mdf',F);
   isproper(M);

elseif ni>=1 & isa(varargin{1},'sym'),    % symbolic argument

   % Conversion of a symbolic fraction from the Symbolic Toolbox
   % to a matrix-den fraction
   % Author: D. Henrion, July 31, 2000

   symb = varargin{1};
   var = '';

   % parse input arguments: variable symbol
   if ni>2,
    error('Too many input arguments.');
   elseif ni==2,
    arg = varargin{2};
    if ~isa(arg,'char'),
       eval('arg = pol(arg);',...
          'error(''Invalid 2nd argument; must be a valid variable symbol.'');');
       [vs1,vs2,vd] = size(arg);
       if all([vs1,vs2,vd]==1) & all(arg.c(:,:)==[0,1]),
          arg = arg.v;
       else
          error('Invalid 2nd argument; must be a valid variable symbol.');
       end;
    end;
    var = arg;
   end;
      
   % unspecified variable symbol: use that of symbolic expression
   if isempty(var),
    var = findsym(symb);
   end;
 
   % extract numerator and denominator with Symbolic Toolbox macro NUMDEN
   [num, den] = numden(symb);

   % conversion from symbolic to polynomial
   polynum = 0; polyden = 0;
   eval('polynum = pol(num,var);','error(peel(lasterr));');
   eval('polyden = pol(den,var);','error(peel(lasterr));');

   % creation of MDF object
   M = 0;
   eval('M = mdf(polynum,polyden);','error(peel(lasterr));');
   switch var,
   case {'s','p','q','z'}
      M.frac.v = var;
   case {'d'}
      M.frac.v = 'd'; M = reverse(M);
   case {'zi','z^-1'}
      M.frac.v = 'z'; M = reverse(M);
   end;

else               % num,den arguments
   N = 0;
   eval('N = pol(varargin{1});', ...
      'error(''Invalid 1st argument.'');');
   if ni==1,
      D = pol(1); PropValStart = 2;
   else
      D = 0;
      eval('D = pol(varargin{2});', ...
         'error(''Invalid 2nd argument.'');');
      PropValStart = 3;
   end;
   
   [td,N,D] = testdp(N,D);
   if ~td, error('Matrices not of the same dimensions.');
   end;
   
   if ~all(all(isfinite(N))) | ~all(all(isfinite(D))),
      error('Polynomial is not finite.');
   end;
   if any(any(D==0)),
      error('Denominator is zero.');
   end;
   
   F = frac(N,D);
   M.title = 'mdf';
   superiorto('double','pol','tsp');
   M = class(M,'mdf',F);
   
   if PropValStart<=ni,
      props(M,varargin{PropValStart:ni});
   end;
      
end;

%end .. @mdf/mdf
