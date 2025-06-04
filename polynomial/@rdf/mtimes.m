function H = mtimes(F,G,arg3,arg4);
%MTIMES (*)   Matrix multiply right-den fractions
%                   H = F*G
%
% The command  H = F*G  or  H = MTIMES(F,G)  is the matrix product
% of right-denominator fractions F and G. Any scalar may
% multiply anything. Otherwise, the number of columns in the
% numerator of F must equal the number of rows in the numerator of G.
% The routine runs using the global zeroing tolerance. The resulting
% fraction is right (it corresponds to the first or to the only
% fractional input argument).
%
% The variable symbols of F and G should be the same. When not,
% a warning is issued and both symbols are changed to the standard one.
% However, if one symbol is 'z' and the other 'z^-1' then the symbols
% play a role, no warning, the resulting symbol is taken from F.
%
% An optional numerical input argument may specify the tolerance to
% be used. An optional 'ldf' or 'rdf' input argument may specify
% the required class of the resulting fraction.
%
% See also RDF/TIMES.

%         Author:  J. Jezek  18-Nov-1999
%         Copyright(c) 1999 by Polyx, Ltd.
%         $ Revision $  $ Date 25-Apr-2000 $
%                       $ Date 30-May-2000 $
%                       $ Date 01-Nov-2000 $
%                       $ Date 30-Sep-2002 $
%                       $ Date 14-Oct-2002 $

global PGLOBAL;

tol = PGLOBAL.ZEROING; hand = 'rdf';
na = nargin;
if na<2,
   error('Not enough input arguments.');
elseif na>2 & ~isempty(arg3),
   if isa(arg3,'double'), tol = arg3;
   elseif isa(arg3,'char'), hand = arg3;
   else error('Invalid 3rd argument.');
   end;
end;
if na>3 & ~isempty(arg4),
   if isa(arg4,'double'), tol = arg4;
   elseif isa(arg4,'char'), hand = arg4;
   else error('Invalid 4th argument.');
   end;
end;

eval('F = rdf(F);','error(peel(lasterr));');
Fs = size(F);

if strcmp(hand,'rdf'),
   if isa(G,'ldf'),
      [tv,Hv,F,G] = testvf(F,G);
      if ~tv, warning('Inconsistent variables.');
      end;
      [th,Hh,F,G] = testhf(F,G,Hv);
      if ~th, warning('Inconsistent sampling periods.');
      end;
      Fn = F.frac.num; Fd = F.frac.den;
      Gn = G.num; Gd = G.den; GdFd = 0;
      eval('GdFd = mtimes(Gd,Fd,tol);','error(peel(lasterr));');
      Gs = size(G);
      if Fs(2)~=Gs(1),
         if all(Fs==1),
            H = rdf(ldf(GdFd,mtimes(Fn,Gn,tol)));
         elseif all(Gs==1),
            H = rdf(mtimes(Fn,Gn,tol),GdFd);
         else
            error('Matrices of inconsistent dimensions.');
         end;
      else
         [X Y] = axby0(GdFd,-Gn,tol);
         Hn = mtimes(Fn,X,tol);
         H = rdf(Hn,Y);
      end;
      if strcmp(F.frac.p,'prop') & strcmp(G.p,'prop') & ...
            F.frac.tp==G.tp,
         props(H,'prop',F.frac.tp);
      end;
   else
      eval('G = rdf(G);','error(peel(lasterr));');
      [tv,Hv,F,G] = testvf(F,G);
      if ~tv, warning('Inconsistent variables.');
      end;
      [th,Hh,F,G] = testhf(F,G,Hv);
      if ~th, warning('Inconsistent sampling periods.');
      end;
      Fn = F.frac.num; Fd = F.frac.den;
      Gn = G.frac.num; Gd = G.frac.den;
      Gs = size(G);
      if Fs(2)~=Gs(1),
         if all(Fs==1) | all(Gs==1),
            FnGn = 0;
            eval('FnGn = mtimes(Fn,Gn,tol);','error(peel(lasterr));');
            H = rdf(FnGn,mtimes(Fd,Gd,tol));
         else          
            error('Matrices of inconsistent dimensions.');
         end;
      else
         X = 0; Y = 0;
         eval('[X Y] = axby0(Fd,-Gn,tol);','error(peel(lasterr));');
         Hn = mtimes(Fn,X,tol); Hd = mtimes(Gd,Y,tol);
         H = rdf(Hn,Hd);
      end;
      if strcmp(F.frac.p,'prop') & strcmp(G.frac.p,'prop') & ...
            F.frac.tp==G.frac.tp,
         props(H,'prop',F.frac.tp);
      end;
   end;
elseif strcmp(hand,'ldf'),
   eval('G = ldf(G);','error(peel(lasterr));');
   [tv,Hv,F,G] = testvf(F,G);
   if ~tv, warning('Inconsistent variables.');
   end;
   [th,Hh,F,G] = testhf(F,G,Hv);
   if ~th, warning('Inconsistent sampling periods.');
   end;
   Fn = F.frac.num; Fd = F.frac.den;
   Gn = G.num; Gd = G.den; GdFd = 0;
   eval('GdFd = mtimes(Gd,Fd,tol);','error(peel(lasterr));');
   Gs = size(G);
   if Fs(2)~=Gs(1),
      if all(Fs==1),
         H = ldf(GdFd,mtimes(Fn,Gn,tol));
      elseif all(Gs==1),
         H = ldf(rdf(mtimes(Fn,Gn,tol),GdFd));
         return;
      else          
         error('Matrices of inconsistent dimensions.');
      end;
   else
      [X Y] = xayb0(GdFd,-Fn,tol);
      Hn = mtimes(X,Gn,tol);
      H = ldf(Y,Hn);
   end;
   if strcmp(F.frac.p,'prop') & strcmp(G.p,'prop') & ...
         F.frac.tp==G.tp,
      props(H,'prop',F.frac.tp);
   end;
else
   error('Invalid command option.');
end;

if strcmp(PGLOBAL.COPRIME,'cop'), H = coprime(H,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'), H = reduce(H,tol);
else H = smreduce(H);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), H = defract(H);
end;

%end .. @rdf/mtimes   
