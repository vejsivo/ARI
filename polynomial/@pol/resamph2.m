function G = resamph2(F,varargin)
%RESAMPH2     Resample and second order hold discrete time polynomial
%
% For polynomial F in variable 'z','q','z^-1 or 'd',
% the command  G = RESAMPH2(F,K,H,L)  returns polynomial
% or fraction G, the result of second order holding with
% interval L and resampling with period K and phase H.
%
% For more details, see FRAC/RESAMPH2.
%
%    Author: J.Jezek, 30-Jun-2000
%    Copyright(c) 2000 by Polyx,Ltd.
%    $ Revision $  $ Date 04-Oct-2000 $
%                  $ Date 02-Feb-2001 $

global PGLOBAL;

if nargin<1,
   error('Not enough input arguments.');
end;
if ~isa(F,'pol'),
   error('Some argument but not 1st is invalidly polynomial.');
end;

Fv = F.v;
if isempty(Fv), Fv = PGLOBAL.DISCRVAR;
end;
lv = length(varargin);

switch Fv,
case {'z^-1','d'},
   k = 2; h = 0; l = k;
   if lv>=1,
      arg1 = varargin{1};
      if isa(arg1,'double'),
         if ~isempty(arg1),
            if length(arg1)==1 & isreal(arg1) & ...
                  arg1>0 & floor(arg1)==arg1,
               k = arg1; l = k;
            else
               error('Invalid resampling period.');
            end;
         end;
         if lv>=2,
            arg2 = varargin{2};
            if isa(arg2,'double'),
               if ~isempty(arg2),
                  if ndims(arg2)==2 & any(size(arg2)==1) & ...
                        isreal(arg2) & all(floor(arg2)==arg2) & ...
                        all(arg2>=0),
                     h = arg2;
                  else
                     error('Invalid resampling phase.');
                  end;
               end;
               if lv>=3,
                  arg3 = varargin{3};
                  if isa(arg3,'double'),
                     if ~isempty(arg3),
                        if length(arg3)==1 & isreal(arg3) & ...
                              floor(arg3)==arg3 & arg3>=1,
                           l = arg3;
                        else
                           error('Invalid holding interval.');
                        end;
                     end;
                  end;
               end;
            end;
         end;
      end;   
   end;
   
   Fh = F.h;
   if isempty(Fh) | ~isfinite(Fh), Fh = 1;
   end;
   switch Fv,
   case 'z^-1', zizi = zi(Fh);
      F = F*rdiv(zizi^2*(1-zizi^l)^3,(1-zizi)^3);
   case 'd', dd = d(Fh);
      F = F*rdiv(dd^2*(1-dd^l)^3,(1-dd)^3);
   end;
   F.h = Fh;
   G = resamp(F,k,h);
   
case {'z','q'},
   F = sdf(F);
   eval('G = resamph2(F,varargin{1:lv});',...
      'error(peel(lasterr));');
   
otherwise
   error('Invalid variable symbol; must be discrete time.');
end;

%end .. @pol/resamph2
