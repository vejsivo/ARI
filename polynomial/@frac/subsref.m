function P = subsref(R,S)
%SUBSREF    Subscripted reference for fraction
%
% For fraction R,
%   R.num  returns the numerator 
%   R.den          the denominator
%   R.s            the size
%   R.v            the variable symbol
%   R.h            the sampling period
%   R.u            the user data
%   R.c            the coprime flag
%   R.r            the reduced flag
%   R.p            the proper  flag
%   R.tc           the coprime tolerance
%   R.tp           the proper  tolerance
%   R.version      the version, always 3.0
%
% For the numerator and denominator, the subscripting in reference
% may be further elaborated, as for polynomials, e.g.  R.num(i,j)
% or  R.den{k} 

%     Author:   J.Jezek  18-Feb-2000
%     Copyright(c) 2000 by Polyx, Ltd.
%     $ Revision $  $ Date 06-Oct-2002 $
%                   $ Date 14-Oct-2002 $

St1 = S(1).type;
if strcmp(St1,'.'),
   flstr = S(1).subs;
   if strcmp(flstr,'version'), P = R.version;
   elseif strcmp(flstr,'den'), P = R.den;
   else   
      switch lower(flstr(1)),
      case 'n', P = R.num;
      case 's', P = R.s;
      case 'v', P = R.v;
      case 'h', P = R.h;
      case 'u', P = R.u;
      case 'c', P = R.c;
      case 'r', P = R.r;
      case 'p', P = R.p;
      case 't',
         if strcmp(flstr,'tc'), P = R.tc;
         elseif strcmp(flstr,'tp'), P = R.tp;
         else error('Invalid field in subscripted reference.');
         end;
      otherwise
         error('Invalid field in subscripted reference.');
      end;
   end;
else
   error('Not implemented subscripted reference.');
end;
if length(S)>1,
   eval('P = subsref(P,S(2:end));','error(peel(lasterr));');
end;
      
%end .. @frac/subsref
