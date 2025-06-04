function [D,L] = deg(A,varargin)
%DEG     Degree of scalar-den fraction
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
% See also POL/DEG, POL/LCOEF, SDF/LCOEF.

%       Author:  J. Jezek, 28-Oct-2002
%       Copyright(c) 2002 by Polyx, Ltd.
%       $ Revision $  $ Date 17-Nov-2002  z,z^-1  $

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
      eval('D = tdeg(A,string);','error(peel(lasterr));');
   else
      eval('[D,L] = tdeg(A,string);','error(peel(lasterr));');
   end;
   D = -D;
   return;
end;

Nd = deg(A.frac.num); Dd = deg(A.frac.den);
Anc = A.frac.num.c; Adc = A.frac.den.c;
Adcl = Adc(1,1,Dd+1);
As1 = A.frac.s(1); As2 = A.frac.s(2);
L = zeros(As1,As2);

if strcmp(string,'ent'),
   if As1==1 & As2>0, string = 'col';
   elseif As2==1 & As1>0, string = 'row';
   end;
end;

switch string,
case 'mat',
   if isempty(Nd), D = []; return;
   end;
   D = Nd-Dd;
   if (no==2) & isfinite(Nd),
      L = Anc(:,:,Nd+1)/Adcl;
   end;
case 'ent',
   if isempty(Nd), D = zeros(As1,As2); return;
   end;
   if no~=2,
      D = deg(A.frac.num,'ent') - Dd;
   else
      [D,L] = deg(A.frac.num,'ent');
      D = D-Dd; L = L/Adcl;
   end;
case 'row',
   if isempty(Nd), D = zeros(As1,0); return;
   end;
   if no~=2,
      D = deg(A.frac.num,'row') - Dd;
   else
      [D,L] = deg(A.frac.num,'row');
      D = D-Dd; L = L/Adcl;
   end;
case 'col',
   if isempty(Nd), D = zeros(0,As2); return;
   end;
   if no~=2,
      D = deg(A.frac.num,'col') - Dd;
   else
      [D,L] = deg(A.frac.num,'col');
      D = D-Dd; L = L/Adcl;
   end;
otherwise,
   error('Invalid command option.');
end;  %switch string

%end .. @sdf/deg
