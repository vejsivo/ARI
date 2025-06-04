function [L,D] = tcoef(A,varargin)
%TCOEF   Trailing coefficient of scalar-den fraction
%
% Fraction A(v) is treated as a function of complex
% variable 'v' and its Laurent series is taken.
% The trailing coefficient TCOEF(A) is the coefficient
% with the lowest power of 'v' in the series.
% When 'v' is 'z','q','s' or 'p', the series is taken
% in point Inf, when 'z^-1' or 'd', in point 0.
%
% In the command  D = TCOEF(A,STRING) , the second input
% argument STRING means:
%   D = TCOEF(A)        default, the same as TCOEF(A,'mat').
%   D = TCOEF(A,'mat')  the trailing coefficient matrix of A.
%   D = TCOEF(A,'ent')  the matrix of scalar trailing coefficients
%                         of A entries.
%   D = TCOEF(A,'row')  the row trailing coefficient matrix of A.
%   D = TCOEF(A,'col')  the column trailing coefficient matrix of A.
%
% For fractions in 'z' or 'z^-1', it is possible to specify
% (by an optional input argument VAR) that the trailing coefficient
% is to be understood by the lowest power of 'z' or 'z^-1'.

% In the command  [L,D] - TCOEF(A,STRING) , the second
% output argument D returns the corresponding matrix or vector of
% trailing degrees, the same as the first output argument in
% function TDEG.
%
% See also POL/TCOEF, POL/TDEG, SDF/TDEG.

%       Author:  J. Jezek, 28-Oct-2002
%       Copyright(c) 2002 by Polyx, Ltd.
%       % Revision $  $ Date 07-Nov-2002 $
%                     $ Date 17-Nov-2002  z,z^-1  $

if ~isa(A,'sdf'),
   error('Some argument but not 1st is invalidly sdf.');
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
      eval('L = lcoef(A,string);','error(peel(lasterr));');
   else
      eval('[L,D] = lcoef(A,string);','error(peel(lasterr));');
      D = -D;
   end;
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
         L = tcoef(A);
      else
         [L,D] = tcoef(A);
      end;
   else
      L = NaN*L; D = -Inf;
   end;
   
case 'ent',
   D = zeros(As1,As2);
   if isempty(Nd), return;
   end;
   den = A.frac.den;
   for i = 1:As1,
      for j = 1:As2,
         num = A.frac.num(i,j);
         Aij = sdf(num,den);
         eval('Aij = tsp(Aij); ok = logical(1);','ok = logical(0);');
         if ok,
            if no~=2,
               L(i,j) = tcoef(Aij);
            else
               [L(i,j),D(i,j)] = tcoef(Aij);
            end;
         else
            L(i,j) = NaN; D(i,j) = -Inf;
         end;
      end;
   end;
   
case 'col',
   D = zeros(1,As2);
   if isempty(Nd), return;
   end;
   den = A.frac.den;
   for j = 1:As2,
      Aj = sdf(A.frac.num(:,j),den);
      eval('Aj = tsp(Aj); ok = logical(1);','ok = logical(0);');
      if ok,
         if no~=2,
            L(:,j) = tcoef(Aj);
         else
            [L(:,j),D(j)] = tcoef(Aj);
         end;
      else
         L(:,j) = NaN; D(j) = -Inf;
      end;
   end;
   
case 'row',
   D = zeros(As1,1);
   if isempty(Nd), return;
   end;
   den = A.frac.den;
   for i = 1:As1,
      Ai = sdf(A.frac.num(i,:),den);
      eval('Ai = tsp(Ai); ok = logical(1);','ok = logical(0);');
      if ok,
         if no~=2,
            L(i,:) = tcoef(Ai);
         else
            [L(i,:),D(i)] = tcoef(Ai);
         end;
      else
         L(i,:) = NaN; D(i) = -Inf;
      end;
   end;
   
otherwise,
   error('Invalid command option.');
end;  %switch string

%end .. @sdf/tcoef
