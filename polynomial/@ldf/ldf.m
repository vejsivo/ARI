function M = ldf(varargin)
%LDF     Create left-denominator fraction
%
% The commmand
%    M = LDF(D,N)
% creates a left-denominator matrix fraction
%    M(v) = (D(v))^-1 * N(v)
% from the square nonsingular polynomial matrix D (denominator)
% and the polynomial matrix N (numerator). The sizes of matrices
% N,D must be compatible unless one of them is scalar. Polynomials
% N,D must be finite (Inf, NaN not allowed).
%
% This syntax may be followed by PropertyValue arguments.
% See PROPS for details on assignable properties.
%
% The commmand
%    M = LDF(N) 
% returns M = N if N is already a left-den fraction. Otherwise,
% it means  M = LDF(1,N).
%
% Both polynomials D,N should have the same variable symbol. When not,
% a warning is issued and 'v' copies the current value of 
% PGLOBAL.VARIABLE. However, if one of the symbols is 'z' and the
% other 'z^-1' then they play a role, no warning, 'v' is taken from N.
%
% A left-den fraction from polynomials D,N can be also created by
% operation  D\N .
%
% The command
%    M = LDF(A,B,C)  or  M = LDF(A,B,C,D)  or  M = LDF(A,B,C,D,E)
% creates left-den fraction M from state space matrices A,B,C and
% possibly D,E. In this case, an optional input argument TOL may
% specify the zeroing tolerance to be used instead the standard
% one; this argument must have a form of character string, e.g.
% LDF(A,B,C,'10^-6'). Also in this case, an optional input 
% argument VAR may specify the variable symbol to be used
% instead the standard one. It must have a form of character
% string; possible values are 's', 'p', 'z', 'zi', 'z^-1', 'q'
% or 'd'.
%
% The command
%    M = LDF(SYS)
% creates left-den fraction M from Control System Toolbox object
% SYS. Also in this case, the optional arguments TOL or VAR
% may be used.
%
% The command
%    M = LDF(SYMB)
% creates left-den fraction M from Symbolic Toolbox object SYMB.
% Optional argument VAR may be used.
%
% See also POL/MLDIVIDE, POL/LDF, RDF/LDF.

% For the internal structure of object "left-denominator
% fraction", see FRAC/FRAC.
%
%       Author(s):  S. Pejchova, M. Sebek 14-7-99
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 15-Jul-1999 16:20:11   $
%                         $Date: 15-Jun-2000  J. Jezek  $
%                         $Date: 31-Jul-2000  D.Henrion $
%                         $Date: 03-Aug-2000  J.Jezek   $
%                         $Date: 02-Nov-2000  J.Jezek   $
%                         $Date; 24-Jan-2002  J.Jezek   $
%                         $Date: 28-Feb-2003  J.Jezek   $
%                         $Date: 24-Nov-2009  M.Sebek   $ 
%                           fixed sampling period for "z" 
%                           for state-space input                        


global PGLOBAL;

ni = nargin;
if ni>=1 & isa(varargin{1},'ldf'),      % quick exit
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
      eval('[N,D] = ss2lmf(A,B,C,D,tol);','error(peel(lasterr));');
      
   elseif ni==5,
      D = varargin{4}; E = varargin{5};
      if ~isa(D,'double'),
         error('Invalid 4th argument.');
      end;
      if ~isa(E,'double'),
         error('Invalid 5th argument.');
      end;
      N = 0;
      eval('[N,D] = dss2lmf(A,B,C,D,E,tol);','error(peel(lasterr));');
      
   else
      error('Too many input arguments.');
   end;
   PGLOBAL.VARIABLE = savedvar;
   
   F = frac(N,D);
   M.title = 'ldf';
   superiorto('double','pol','tsp');
   M = class(M,'ldf',F);

   switch var,
%   case {'s','p','q','z'} % M.Sebek, Nov 2009
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
   
   N = 0; D = 0;
   eval('[N,D] = lti2lmf(sys,tol);','error(peel(lasterr));');
   F = frac(N,D);
   M.title = 'ldf';
   superiorto('double','pol','tsp');
   M = class(M,'ldf',F);

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
   M.title = 'ldf';
   superiorto('double','pol','tsp');
   M = class(M,'ldf',F);
   isproper(M);

elseif ni>=1 & isa(varargin{1},'sym'),    % symbolic argument

   % Conversion of a symbolic fraction from the Symbolic Toolbox
   % to a left-den fraction
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

   % creation of LDF object
   M = 0;
   eval('M = ldf(mdf(polynum,polyden));','error(peel(lasterr));');
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
      D = N;
      eval('N = pol(varargin{2});', ...
         'error(''Invalid 2nd argument.'');');
      PropValStart = 3;
   end;
   
   if ~all(all(isfinite(N))) | ~all(all(isfinite(D))),
      error('Polynomial is not finite.');
   end;
   if D.s(1)~= D.s(2),
      error('Denominator matrix is not square.');
   end;
   if issingular(D),
      error('Denominator matrix is singular.');
   end;

   if D.s(2)~=N.s(1),
      if all(D.s==1),
         D = D * eye(N.s(1));
      elseif  all(N.s==1),
         N = N * eye(D.s);
      else
         error('Matrices of inconsistent dimensions.');
      end;
   end;
   
   F = frac(N,D);
   M.title = 'ldf';
   superiorto('double','pol','tsp');
   M = class(M,'ldf',F);
   
   if PropValStart<=ni,
      props(M,varargin{PropValStart:ni});
   end;
      
end;

%end .. @ldf/ldf
