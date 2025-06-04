function [L,D] = lcoef(A,varargin)
%LCOEF   Leading coefficient of scalar-den fration
%
% Fraction A(v) is treated as a function of complex
% variable 'v' and its Laurent series is taken.
% The leading coefficient LCOEF(A) is the coefficient
% with the highest power of 'v' in the series.
% When 'v' is 'z','q','s' or 'p', the series is taken
% in point Inf, when 'z^-1' or 'd', in point 0.
%
% In the command  L = LCOEF(A,STRING) , the second
% input argument STRING means:
%   L = LCOEF(A)        default, the same as LCOEF(A,'mat').
%   L = LCOEF(A,'mat')  the leading coefficient matrix of A.
%   L = LCOEF(A,'ent')  the matrix of scalar leading coefficients
%                         of A entries.
%   L = LCOEF(A,'row')  the row leading coefficient matrix of A.
%   L = LCOEF(A,'col')  the column leading coefficient matrix of A.
%
% For fractions in 'z' or 'z^-1', it is possible to specify
% (by an optional input argument VAR) that the leading coefficient
% is to be understood by the highest power of 'z' or 'z^-1'.

% In the command  [L,D] - LCOEF(A,STRING) , the second
% output argument D returns the corresponding matrix or vector of
% degrees, the same as the first output argument in function DEG.
%
% See also POL/LCOEF, POL/DEG, SDF/DEG.

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
      eval('L = tcoef(A,string);','error(peel(lasterr));');
   else
      eval('[L,D] = tcoef(A,string);','error(peel(lasterr));');
      D = -D;
   end;
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
   L = Anc(:,:,Nd+1)/Adcl; D = Nd-Dd;
case 'ent',
   if isempty(Nd), D = zeros(As1,As2); return;
   end;
   if no~=2,
      L = lcoef(A.frac.num,'ent')/Adcl;
   else
      [L,D] = lcoef(A.frac.num,'ent');
      L = L/Adcl; D = D-Dd;
   end;
case 'row',
   if isempty(Nd), D = zeros(As1,0); return;
   end;
   if no~=2,
      L = lcoef(A.frac.num,'row')/Adcl;
   else
      [L,D] = lcoef(A.frac.num,'row');
      L = L/Adcl; D = D-Dd;
   end;
case 'col',
   if isempty(Nd), D = zeros(0,As2); return;
   end;
   if no~=2,
      L = lcoef(A.frac.num,'col')/Adcl;
   else
      [L,D] = lcoef(A.frac.num,'col');
      L = L/Adcl; D = D-Dd;
   end;
otherwise,
   error('Invalid command option.');
end;  %switch string

%end .. @sdf/lcoef

