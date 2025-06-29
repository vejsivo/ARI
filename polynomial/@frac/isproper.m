function t = isproper(F,arg2,arg3)
%ISPROPER    Test if fraction is proper
%
% For scalar-den fraction or matrix-den fraction F,
% the expression  ISPROPER(F)  returns 1 if F is proper,
% otherwise 0.
%
% For fraction in 's','p','z','q', the properness means the
% behaviour in point Infinity. For fraction in 'z^-1' or 'd',
% in point 0.
%
% An optional argument 'strictly' requires testing whether
% the fraction is strictly proper. It may be shortened down to 's'.
%  
% An optional argument TOL may specify the zeroing tolerance
% to be used instead of the standard one.

%       Author:  J. Jezek  22-Mar-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 09-Jul-2001 $
%                     $ Date 25-Jul-2002 $
%                     $ Date 06-Jul-2002 $
%                     $ Date 14-Oct-2002 $

global PGLOBAL;

strict = 0; tol = PGLOBAL.ZEROING;
ni = nargin;

if ni>=2,
   arg2err = 0;
   if isa(arg2,'char'),
      if isempty(arg2), arg2 = [];
      else
         I = strmatch(arg2,{'strictly';'strictl';'strict';'stric'; ...
            'stri';'str';'st';'s'},'exact');
         if ~isempty(I), strict = 1; arg2 = [];
         else arg2 = str2num(arg2);
            if isempty(arg2), arg2err = 1;
            end;
         end;
      end;
   end;
   if isa(arg2,'double'),
      if ~isempty(arg2),
         if length(arg2)==1 & isreal(arg2) & arg2>=0 & arg2<=1,
            tol = arg2;
         else
            arg2err = 1;
         end;
      end;
   else
      arg2err = 1;
   end;
   if arg2err,
      error('Invalid 2nd argument.');
   end;
end;

if ni>=3,
   arg3err = 0;
   if isa(arg3,'char'),
      if isempty(arg3), arg3 = [];
      else
         I = strmatch(arg3,{'strictly';'strictl';'strict';'stric'; ...
            'stri';'str';'st';'s'},'exact');
         if ~isempty(I), strict = 1; arg3 = [];
         else arg3 = str2num(arg3);
            if isempty(arg3), arg3err = 1;
            end;
         end;
      end;
   end;
   if isa(arg3,'double'),
      if ~isempty(arg3),
         if length(arg3)==1 & isreal(arg3) & arg3>=0 & arg3<=1,
            tol = arg3;
         else
            arg3err = 1;
         end;
      end;
   else
      arg3err = 1;
   end;
   if arg3err,
      error('Invalid 3rd argument.');
   end;
end;

if strcmp(F.p,'nprop') & F.tp==tol,      % quick exit
   t = logical(0); return;
elseif ~strict & strcmp(F.p,'prop') & F.tp==tol,
   t = logical(1); return;
end;

Fcl = class(F);
if strcmp(Fcl,'frac'),
   error('Invalid 1st argument.');
end;

if any(F.s==0),
   t = logical(1);
   
else
   if strcmp(F.v,'z^-1') | strcmp(F.v,'d'),
      FF = reverse(F);
   else FF = F;
   end;

   Nd = deg(FF.num,'ent'); Dd = deg(FF.den,'ent');
   t = (all(all(Nd<=Dd)));
end;

name = inputname(1);
if ~isempty(name),
   if t, props(F,'prop',tol);
   else  props(F,'nprop',tol);
   end;
   assignin('caller',inputname(1),F);
end;

if strict & all(F.s>0),
   t = all(Nd<Dd);
end;

%end .. @frac/isproper
