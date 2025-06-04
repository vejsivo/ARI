function [D,L] = deg(A,varargin)
%DEG     Degree of matrix-den fraction
%
% Fraction A(v) is treated as a function of complex
% variable 'v' and its Laurent series is taken.
% The degree  DEG(A) shows what is the highest power
% of 'v' in the series. When 'v' is 'z','q','s' or 'p',
% the series is taken in point Inf, when 'z^-1' or 'd',
% in point 0. 
%
% In the command  D = DEG(A,STRING) , the second input
% argument STRING means:
%   D = DEG(A)        default, the same as DEG(A,'mat').
%   D = DEG(A,'mat')  the degree of matrix A.
%   D = DEG(A,'ent')  the matrix of degrees of the A entries.
%   D = DEG(A,'row')  the column vector of row degrees of A.
%   D = DEG(A,'col')  the row vector of column degrees of A.
%
% For fractions in 'z' or 'z^-1', it is possible to specify
% (by an optional input argument VAR) that the degree is to be
% understood by the highest power of 'z' or 'z^-1'.
%
% In the command  [D,L] = DEG(A,STRING) , the second
% output argument L returns the correspoding matrix
% of leading coefficients, the same as first output
% argument in function LCOEF. 
%
% See also POL/DEG, POL/LCOEF, MDF/LCOEF.

%       Author:  J. Jezek, 28-Oct-2002
%       Copyright(c) 2002 by Polyx, Ltd.
%       $ Revision $  $ Date 17-Nov-2002  z,z^-1  $

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
      eval('D = tdeg(A,string);','error(peel(lasterr));');
   else
      eval('[D,L] = tdeg(A,string);','error(peel(lasterr));');
   end;
   D = -D;
   return;
end;

Nd = deg(A.frac.num);
As1 = A.frac.s(1); As2 = A.frac.s(2);
L = zeros(As1,As2);

if no~=2,
   Dn = deg(A.frac.num,'ent');
   Dd = deg(A.frac.den,'ent');
else
   [Dn,Ln] = deg(A.frac.num,'ent');
   [Dd,Ld] = deg(A.frac.den,'ent');
end;
dif = Dn-Dd;

if strcmp(string,'ent'),
   if As1==1 & As2>0, string = 'col';
   elseif As2==1 & As1>0, string = 'row';
   end;
end;

switch string,
case 'mat',
   if isempty(Nd), D = []; return;
   end;
   D = max(max(dif));
   if no==2,
      for i = 1:As1,
         for j = 1:As2,
            if dif(i,j)==D,
               L(i,j) = Ln(i,j)/Ld(i,j);
            end;
         end;
      end;
   end;

case 'ent',
   if isempty(Nd), D = zeros(As1,As2); return;
   end;
   D = dif;
   if no==2,
      L = Ln./Ld;
   end;
   
case 'row',
   if isempty(Nd), D = zeros(As1,0); return;
   end;
   D = max(dif,[],2);
   if no==2,
      for i = 1:As1,
         Di = D(i);
         for j = 1:As2,
            if dif(i,j)==Di,
               L(i,j) = Ln(i,j)/Ld(i,j);
            end;
         end;
      end;
   end;
   
case 'col',
   if isempty(Nd), D = zeros(0,As2); return;
   end;
   D = max(dif,[],1);
   if no==2,
      for j = 1:As2,
         Dj = D(j);
         for i = 1:As1,
            if dif(i,j)==Dj,
               L(i,j) = Ln(i,j)/Ld(i,j);
            end;
         end;
      end;
   end;
      
otherwise,
   error('Invalid command option.');
end;  %switch string

%end .. @mdf/deg
   