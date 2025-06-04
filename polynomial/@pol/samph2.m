function G = samph2(F,T,varargin)
%SAMPH2     Sample and second order hold continuous time polynomial
%
% For polynomial F with variable symbol 's' or 'p', the command
%  G = SAMPH2(F,T,TAU,L)  returns polynomial or fraction G,
% the result of second order holding with interval L and
% sampling with period T and phase TAU.
%
% When some of the sampling points TAU, T+TAU, T+2*TAU...
% coincides with 0, L,2*L or 3*L, the degree of F must not
% be greater than 2 (otherwise, in time domain, sampling
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
% For more details, see FRAC/SAMPH2, POL/SAMPH, POL/SAMP.
%
%    Author: J.Jezek, 13-Jun-2000
%    Copyright(c) 2000 by Polyx,Ltd.
%    $ Revision $  $ Date 04-Oct-2000 $
%                  $ Date 01-Dec-2000 $
%                  $ Date 16-Jan-2001 $

global PGLOBAL;

ni = nargin;
if ni<2,
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
if dF>=1, Fc1 = F.c(:,:,2);
end;
if dF>=2, Fc2 = F.c(:,:,3);
end;

for k = 1:ltau,
   tauk = tau(k);
   D = tauk/T; D1 = (L-tauk)/T; D2 = (2*L-tauk)/T; D3 = (3*L-tauk)/T;
   if dF>2 & (D<=tol | D1-floor(D1)<=tol | floor(D1)+1-D1<tol | ...
         D2-floor(D2)<=tol | floor(D2)+1-D2<=tol | ...
         D3-floor(D3)<=tol | floor(D3)+1-D3<=tol),
      error('In time domain, sampling of Dirac impulse.');
   end;

   Q1 = floor((L-tol-tauk)/T); Q2 = floor((2*L-tol-tauk)/T);
   Q3 = floor((3*L-tol-tauk)/T);
   if Q3<0, return;
   end;
   p1 = 0:Q1; p2 = Q1+1:Q2; p3 = Q2+1:Q3;
   t1 = tauk+p1*T; t2 = tauk+p2*T; t3 = tauk+p3*T;
   L2 = L^2;
   v1 = t1.^2/2; v2 = (-2*t2.^2 + 6*L*t2 - 3*L2)/2;
   v3 = (t3.^2 - 6*L*t3 + 9*L2)/2;
   v = [v1,v2,v3];
   Gc = kron(v,Fc0);
   if dF>=1,
      v1 = t1; v2 = -2*t2 + 3*L; v3 = t3-3*L;
      v = [v1,v2,v3];
      Gc1 = kron(v,Fc1);
      Gc = Gc+Gc1;
      if dF>=2,
         v1(:) = 1; v2(:) = -2; v3(:) = 1;
         v = [v1,v2,v3];
         Gc2 = kron(v,Fc2);
         Gc = Gc+Gc2;
      end;
   end;

   switch Gv,
   case {'z^-1','d'},
      GG = pol(Gc,Q3,Gv);
   case 'z',
      GG = pol(Gc,Q3,'z^-1'); GG = reverse(sdf(GG));
      GG = defract(GG);
   case 'q',
      GG = pol(Gc,Q3,'z^-1'); GG = reverse(sdf(GG)); GG.v = 'q';
      GG = defract(GG);
   end;
   GG.h = T;
   if ltau==1, G = GG;
   else G{k} = GG;
   end;
end;

%end .. @pol/samph2

    

