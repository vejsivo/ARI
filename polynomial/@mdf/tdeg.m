function [D,L] = tdeg(A,varargin)
%TDEG    Trailing degree of matrix-den fraction
%
% Fraction A(v) is treated as a function of complex
% variable 'v' and its Laurent series is taken.
% The trailing degree TDEG(A) shows what is the
% lowest power of 'v' in the series. When 'v' is 'z','q',
% 's' or 'p', the series is taken in point Inf,
% when 'z^-1' or 'd', in point 0.
%
% In the command  D = TDEG(A,STRING) , the second input
% argument STRING means:
%   D = TDEG(A)        default, the same as TDEG(A,'mat').
%   D = TDEG(A,'mat')  the trailing degree of matrix A.
%   D = TDEG(A,'ent')  the matrix of trailing degrees of the A entries.
%   D = TDEG(A,'row')  the column vector of row trailing degrees of A.
%   D = TDEG(A,'col')  the row vector of column trailing degrees of A.
%
% For fractions in 'z' or 'z^-1', it is possible to specify
% (by an optional input argument VAR) that the triling degree is to be
% understood by the highest power of 'z' or 'z^-1'.
%
% In the command  [D,L] = TDEG(A,STRING) , the second
% output argument L returns the correspoding matrix
% of leading coefficients, the same as first output
% argument in function TCOEF.
%
% See also POL/TDEG, POL/TCOEF, MDF/TCOEF.

%       Author:  J. Jezek, 28-Oct-2002
%       Copyright(c) 2002 by Polyx, Ltd.
%       $ Revision $  $ Date 07-Nov-2002 $
%                     $ Date 17-Nov-2002  z,z^-1  $

if ~isa(A,'mdf'),
   error('Some argument but not 1st is invalidly mdf.');
end;
Av = A.frac.v; var = Av; 
if strcmp(Av,'z^-1') | strcmp(Av,'d'),
   A = reverse(A);
end;

string = 'mat';
li = length(varargin);
if li>0,
   for i = 1:li,
      arg = varargin{i};
      if isa(arg,'pol'),
         [vs1,vs2,vd] = size(arg);
         if all([vs1,vs2,vd]==1) & all(arg.c(:,:)==[0 1]),
            arg = arg.v;
         else
            error(['Invalid ',nth(i+1),' argument.']);
         end;
      end;
      if ischar(arg),
         if strcmp(arg,'zi'), arg = 'z^-1';
         end;
         I = strmatch(arg,{'s';'p';'z^-1';'d';'z';'q'},'exact');
         if ~isempty(I),
            var = arg;
            if ~isempty(Av) & ~strcmp(Av,var),
               if (strcmp(Av,'z') & strcmp(var,'z^-1')) | ...
                  (strcmp(Av,'z^-1') & strcmp(var,'z')),
               else
                  error('Invalid variable symbol.');
               end;
            end; 
         else
            string = arg;
         end;
      else
         error(['Invalid ',nth(i+1),' argument.']);
      end;
   end;
end;        
   
no = nargout;
if strcmp(var,'z^-1'),
   if no<2,
      eval('D = deg(A,string);','error(peel(lasterr));');
   else
      eval('[D,L] = deg(A,string);','error(peel(lasterr));');
   end;
   D = -D;
   return;
end;

Nd = deg(A.frac.num);
As1 = A.frac.s(1); As2 = A.frac.s(2);
L = zeros(As1,As2);
A.frac.v = 'z';

if strcmp(string,'ent'),
   if As1==1 & As2>0, string = 'col';
   elseif As2==1 & As1>0, string = 'row';
   end;
end;

switch string,
case 'mat',
   if isempty(Nd),
      D = [];
      return;
   end;
   eval('A = tsp(A); ok = logical(1);','ok = logical(0);');
   if ok,
      if no~=2,
         D = tdeg(A);
      else
         [D,L] = tdeg(A);
      end;
   else
      D = -Inf; L = NaN*L;
   end;
   
case 'ent',
   D = zeros(As1,As2);
   if isempty(Nd), return;
   end;
   for i = 1:As1,
      for j = 1:As2,
         num = A.frac.num(i,j); den = A.frac.den(i,j);
         Aij = mdf(num,den);
         eval('Aij = tsp(Aij); ok = logical(1);','ok - logical(0);');
         if ok,
            if no~=2,
               D(i,j) = tdeg(Aij);
            else
               [D(i,j),L(i,j)] = tdeg(Aij);
            end;
         else
            D(i,j) = -Inf; L(i,j) = NaN;
         end;
      end;
   end;
   
case 'col',
   D = zeros(1,As2);
   if isempty(Nd), return;
   end;
   for j = 1:As2,
      Aj = mdf(A.frac.num(:,j),A.frac.den(:,j));
      eval('Aj = tsp(Aj); ok = logical(1);','ok = logical(0);');
      if ok,
         if no~=2,
            D(j) = tdeg(Aj);
         else
            [D(j),L(:,j)] = tdeg(Aj);
         end;
      else
         D(j) = -Inf; L(:,j) = NaN;
      end;
   end;
   
case 'row',
   D = zeros(As1,1);
   if isempty(Nd), return;
   end;
   for i = 1:As1,
      Ai = mdf(A.frac.num(i,:),A.frac.den(i,:));
      eval('Ai = tsp(Ai); ok = logical(1);','ok = logical(0);');
      if ok,
         if no~=2,
            D(i) = tdeg(Ai);
         else
            [D(i),L(i,:)] = tdeg(Ai);
         end;
      else
         D(i) = -Inf; L(i,:) = NaN;
      end;
   end;
   
otherwise,
   error('Invalid command option.');
end;  %switch string

%end .. @mdf/tdeg
