function G = samph(F,T,varargin)
%SAMPH    Sample and hold continuous time polynomial
%
% For polynomial with variable symbol 's' or 'p', 
% the command  G = SAMPH(F,T,TAU,L)  returns
% polynomial or fraction G in variable 'z','q','z^-1'
% or 'd', the result of holding with interval L and
% sampling with period T and phase TAU.
%
% Argument TAU can also be, instead of scalar, a vector.
% In such a case, the results form a cell vector.
%
% When some of the sampling points TAU, T+TAU, T+2*TAU...
% coincides with 0 or L, the degree of F must not be
% greater than 0 (otherwise, in time domain, sampling
% of a Dirac impulse would be attempted).
%
% The arguments T,TAU,L are optional, the defaults being
% T=1, TAU=0, L=T. An optional argument TOL may specify the
% zeroing tolerance (for the above mentioned coincidence 
% of points). This argument must have a form of character
% string, e.g. '10^-4'. An optional argument VAR may
% specify the required variable symbol of the discrete
% time result.
%
% For more details, sse FRAC/SAMPH, POL/SAMP.

%        Author: J.Jezek, 03-Jul-2000
%        Copyright(c) 2000 by Polyx, Ltd.
%        $ Revision $  $ Date 04-Oct-2000 $
%                      $ Date 01-Dec-2000 $
%                      $ Date 16-Jan-2001 $

global PGLOBAL;

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
if ~isa(F,'pol'),
   error('Some argument but not 1st is invalidly polynomial.');
end;

Fv = F.v;
if ~isempty(Fv) & ~strcmp(Fv,'s') & ~strcmp(Fv,'p'),
   error('Invalid variable symbol; must be continuous time.');
end;

if ni>=2 & ~isempty(T),
   if ~isa(T,'double') | length(T)~=1 | ~isreal(T) | T<=0,
      error('Invalid sampling period.');
   end;
else T = 1;
end;

tau = 0; L = T;
tol = PGLOBAL.ZEROING; Gv = PGLOBAL.DISCRVAR;
lv = length(varargin); start = 1;
if lv>=1,
   arg1 = varargin{1};
   if isa(arg1,'double'),
      if ~isempty(arg1),
         if ndims(arg1)==2 & any(size(arg1)==1) & ...
               isreal(arg1) & all(arg1>=0),
            tau = sort(arg1);
         else
            error('Invalid sampling phase.');
         end;
      end;
      
      start = 2;
      if lv>=2,
         arg2 = varargin{2};
         if isa(arg2,'double'),
            if ~isempty(arg2),
               if length(arg2)==1 & isreal(arg2) & arg2>0,
                  L = arg2;
               else
                  error('Invalid holding interval.');
               end;
            end;
            start = 3;
         end;
      end;
   end;
end;

for i = start:lv,
   arg = varargin{i};
   if isa(arg,'pol'),
      [vs1,vs2,vd] = size(arg);
      if all ([vs1,vs2,vd]==1) & all(arg.c(:,:)==[0,1]),
         arg = arg.v;
      else
         error(['Invalid ',nth(i+start-1),' argument.']);
      end;
   end;
   
   if isa(arg,'char'),
      num = str2num(arg);
      if isa(num,'double') & ~isempty(num),
         if length(num)==1 & isreal(num) & num>=0 & num<1,
            tol = num;
         else
            error('Invalid tolerance.');
         end;
      elseif ~isempty(arg),
         if strcmp(arg,'zi'), arg = 'z^-1';
         end;
         I = strmatch(arg,{'z';'q';'z^-1';'d'},'exact');
         if ~isempty(I), Gv = arg;
         else error(['Invalid ',nth(i+start-1),' argument.']);
         end;
      end;
   else
      error(['Invalid ',nth(i+start-1),' argument.']);
   end;
end;

ltau = length(tau); stau = size(tau);
GG = pol(zeros(size(F)));
if ltau==1,
   G = GG;
else
   G = cell(stau);
   for k = 1:ltau,
      G{k} = GG;
   end;
end;
if isempty(F) | all(all(F==0)),
   return;
end;

dF = deg(F); Fc0 = F.c(:,:,1);
for k = 1:ltau,
   tauk = tau(k);
   D = tauk/T; D1 = (L-tauk)/T;
   if dF>0 & (D<=tol | D1-floor(D1)<=tol | floor(D1)+1-D1<tol),
      error('In time domain, sampling of Dirac impulse.');
   end;

   Q = floor((L-tol-tauk)/T);
   if Q>=0,
      Gc = repmat(Fc0,1,Q+1);
      switch Gv,
      case {'z^-1','d'},
         GG = pol(Gc,Q,Gv);
      case 'z',
         GG = pol(Gc,Q,'z^-1'); GG = reverse(sdf(GG));
         GG = defract(GG);
      case 'q',
         GG = pol(Gc,Q,'z^-1'); GG = reverse(sdf(GG)); GG.v = 'q';
         GG = defract(GG);
      end;
      GG.h = T;
      if ltau==1, G = GG;
      else G{k} = GG;
      end;
   end;
end;

%end .. @pol/samph

    

