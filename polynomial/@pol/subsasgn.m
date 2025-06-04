function PP = subsasgn(PP,S,R)
%SUBSASGN  Subscripted assignment for polynomial
%
% If P is the polynomial matrix 
%    P0 + P1*s + P2*s^2 + ... + Pd*s^d
% then
%    P{0:d} = R is [P0 P1 P2 ... Pd] = R,
%    P{k} = R is the  k-term of P, Pk = R.
%
% P(ri,ci) = R sets the polynomial submatrix of P of row indices ri
% and column indices ci equal to R. 
%
% P.var  = R  sets v, the indeterminate string variable of P.
% P.h    = R  sets h, the sampling period of P.
% P.user = R  sets u, the user data of P.
%
% The above subscripted assignments may be combined, i.e.,
% P{di}(ri,ci) = R  sets the constant coefficient matrix
% corresponding to a submatrix P equal to R.
%
% See also POL/SUBSREF.

%       Author(s):  S. Pejchova, M. Sebek 5-3-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 14-Apr-2000 10:47:34   $
%       $Revision: 3.0 $  $Date: 20-Jun-2001 10:45:00  J.Jezek  $
%                         $Date: 02-Feb-2002           J.Jezek  $
%                         $Date: 28-Feb-2002           J.Jezek  $
%                         $Date: 28-Jun-2002           J.Jezek  $
%                         $Date: 28-Feb-2003           J.Jezek  $
%                         $Date: 27-Apr-2003           J.Jezek  $

global PGLOBAL;

if nargin<3,
   error('Not enough input arguments.');
end;
if ~isstruct(S),
   error('Invalid 2nd argument.');
end;

S_lg = length(S); S1 = S(1);
St1 = S1.type; Ss1 = S1.subs;
Ss1_lg = length(Ss1);
r_ind = NaN; d_ind = NaN; one_ind = 1;
S_par.type = '()';
S_par.subs{1} = ':';
S_par.subs{2} = ':';
S_par.subs{3} = ':';

switch St1,
case '.',
   flstr = lower(Ss1); flstr1 = flstr(1);
   switch flstr1,
   case 'f',
      if S_lg==1,
         error('Invalid field in subscripted assignment.');
      end;
      S(1) = []; St1 = S(1).type; Ss1 = S(1).subs;
      if ~strcmp(St1,'.'),
         error('Invalid combination in subscripted assignment.');
      end;
      flstr = lower(Ss1); flstr1 = flstr(1);
      if ~strcmp(flstr1,'n'),
         error('Invalid field in subscripted assignment.');
      end;
      S(1) = [];
      eval(['PP',expand(S),'=R;'],'error(peel(lasterr));');
   case 'n',
      S(1) = [];
      eval(['PP',expand(S),'=R;'],'error(peel(lasterr));');      
   case 'v',
      if S_lg>1,
         error('Invalid combination in subscripted assignment.');
      end;
      if isa(R,'pol'),
         [vs1,vs2,vd] = size(R);
         if all([vs1,vs2,vd]==1)&(~any(R.c(:,:)-[0,1])),
            R = R.v;
         else
            error('Invalid variable symbol.');
         end;
      end;
      if isstr(R),
         if isempty(R),
            if (isempty(PP.d)|(PP.d<=0)), PP.v = '';
            else error('Invalid variable symbol.');
            end;
         else
            if strcmp(R,'zi'), R = 'z^-1';
            end;
            I = strmatch(R,{'s';'p';'z^-1';'d';'z';'q'},'exact');
            if ~isempty(I) | ((isempty(PP.d)|(PP.d<=0)) & isempty(R)),
               Oldv = PP.v; PP.v = R;
               J = strmatch(R,{'z^-1';'d';'z';'q'},'exact');
               if ~isempty(J) & (strcmp(Oldv,'s') | strcmp(Oldv,'p')),
                  PP.h = 1;
               end;
            else,
               error('Invalid variable symbol.');
            end;
            if (isempty(PP.d)|(PP.d<=0)), PP.v = '';
            end;
         end;
         if isempty(PP.v), PP.h = [];
         elseif strcmp(PP.v,'s') | strcmp(PP.v,'p'), PP.h = 0;
         end;
      else, 
         error('Invalid variable symbol.');
      end;
      
   case 'h',
      if S_lg>1,
         error('Invalid combination in subscripted assignment.');
      end;
      if isempty(R),
         I = strmatch(PP.v,{'z^-1';'z';'d';'q'},'exact');
         if ~isempty(I), PP.h = [];
         end;
      elseif ~isa(R,'double') | length(R)~=1 | ...
             ~isreal(R) | R<0 | isinf(R),
         error('Invalid sampling period.');
      elseif strcmp(PP.v,'s') | strcmp(PP.v,'p'),
         if R>0,
            error('Invalid object for nonzero sampling period.');
         end;
         PP.h = R;
      elseif ~isempty(PP.v),
         if R==0,
            error('Invalid object for zero sampling period.');
         end;
         PP.h = R;
      end;
      
   case 'u',
      S(1) = [];
      if isempty(S),
         PP.u = R;
      else
         PPu = PP.u;
         eval('PPu = subsasgn(PPu,S,R);', ...
            ['if isa(PPu,''pol'') | isa(PPu,''tsp'') | isa(PPu,''frac''),', ...
               'error(peel(lasterr)); else error(peel(lasterr)); end;']);
         PP.u = PPu;
      end;
      
   otherwise
      error('Invalid field in subscripted assignment.');
   end;
   
case '()',
   if S_lg>1,
      i = 2; Sti = S(2).type; Ssi = S(2).subs;
      if ~strcmp(Sti,'.'),
         error('Invalid combination in subscripted assignment.');
      end;
      if strcmp(Ssi,'f'),
         if S_lg==2,
            error('Invalid field in subscripted assignment.');
         end;
         i = 3; Sti = S(3).type; Ssi = S(3).subs;
      end;
      if ~strcmp(Ssi,'n'),
         error('Invalid field in subscripted assignment.');
      end;
      while i<S_lg,
         i = i+1; Sti = S(i).type; Ssi = S(i).subs;
         if strcmp(Sti,'.'),
            if ~strcmp(Ssi,'n'),
               error('Invalid field in subscripted assignment.');
            end;
         elseif strcmp(Sti,'{}'),
            eval('R = double(R);', 'error(peel(lasterr));');
            if length(Ssi)>1,
               error('More than one subscript in {} .');
            elseif ~isempty(PP.d) & isempty(Ssi{1}),
               error('Invalid reference to degree.');
            else
               d_ind = Ssi{1};
            end;
         else
            error('Invalid combination in subscripted assignment.');
         end;
         
      end;
   end;
   
   r_ind = 1;
   S_par.subs{1} = Ss1{1};
   if Ss1_lg==2,
      one_ind = 0; S_par.subs{2} = Ss1{2};
   end;  
   if Ss1_lg>2,
      error('More than two subscripts.');
   end;
   
case '{}',
   eval('R = double(R);', 'error(peel(lasterr));');
   if length(Ss1)>1,
      error('More than one subscript in {} .');
   elseif ~isempty(PP.d) & isempty(Ss1{1}),
      error('Invalid reference to degree.');
   else
      d_ind = Ss1{1};
   end;
   if S_lg==2,
      St2 = S(2).type;
      if strcmp(St2,'()'),
         Ss2 = S(2).subs; Ss2_lg = length(Ss2);
         r_ind = 1;
         S_par.subs{1} = Ss2{1};
         if Ss2_lg==2,
            one_ind = 0; S_par.subs{2} = Ss2{2};
         end;
         if Ss2_lg>2,
            error('More than two subscripts.');
         end;
      else
         error('Invalid combination in subscripted assignment.');
      end;
   elseif S_lg>2,
      error('Invalid combination in subscripted assignment.');
   end;
end;

if isnan(r_ind) & isnan(d_ind),
   return;
end;

var1 = PP.v; var2 = '';
if isa(R,'pol'),
   var2 = R.v;
   if isempty(R.d) | isinf(R.d),
      Rc = zeros(R.s);
   else
      Rc = R.c;
   end;
   if (strcmp(var1,'z') & strcmp(var2,'z^-1')) | ...
         (strcmp(var1,'z^-1') & strcmp(var2,'z')),
      R = tsp(R);
   end;
end;
if isa(R,'tsp'),
   PP = tsp(PP);
   eval(['PP',expand(S),'=R;'],'error(peel(lasterr));');
   return;
elseif isa(R,'double') & ndims(R)<4,
   Rc = R;
elseif ~isa(R,'pol'),
   error('Argument cannot be assigned to polynomial.');
end;

Pd = PP.d; ps1 = PP.s(1); ps2 = PP.s(2);
if isempty(Pd) | isinf(Pd),
   Pd = 0; PP.c = zeros(ps1,ps2);
end;
Pc = PP.c; sn1 = ps1; sn2 = ps2;
if isnan(r_ind), one_ind = 0;
end;
if one_ind & ~(isempty(S_par.subs{1})),
   sn1=ps1*ps2; sn2=1;
   Pc = reshape(PP.c,[sn1,sn2,Pd+1]);
   S_par.subs{2}=[1];
end;
[rc1,rc2,rc3]=size(Rc);
max_d=1;
if ~isnan(d_ind),
   S_par.subs{3}=d_ind;
   if isa(d_ind,'double'),
      if isfinite(d_ind), 
         max_d=max(d_ind+1); 
      elseif isinf(PP.d),
         fd_inf=find(isinf(d_ind)&d_ind<0);
         d_ind(fd_inf)=0;
      else
         error('Subscript is not finite.');
      end;
      if d_ind<0,
         if isempty(var1) | strcmp(var1,'z') | strcmp(var1,'z^-1'),
            PP = tsp(PP);
            eval('PP = subsasgn(PP,S,R);','error(peel(lasterr));');
            return;
         else
            error('Subscript in {} is negative.');
         end;
      end;
      S_par.subs{3}=d_ind+1;   
      if ndims(Rc)==2 & rc2,
         ld=length(d_ind);
         rn2=rc2/length(d_ind);
         if rn2==round(rn2),
            Rc=reshape(Rc,rc1,rn2,ld); rc2=rn2; rc3=ld;
         else, error('Subscripted assignment dimension mismatch.');
         end;
      end;         
   end;
elseif ~(isempty(S_par.subs{1})) & ~(isempty(S_par.subs{2})),
   Rc=cat(3,Rc,zeros(rc1,rc2,max([rc3,size(Pc,3)])-rc3)); 
end;
s3=max([rc3,size(Pc,3),max_d]);
Pc=cat(3,Pc,zeros(sn1,sn2,s3-size(Pc,3)));
i1 = S_par.subs{1}; i2 = S_par.subs{2};
if (isnumeric(i1) & ~isempty(i1) & any(i1<=0)) | ...
   (isnumeric(i2) & ~isempty(i2) & any(i2<=0))
      error('Index into matrix is negative or zero.');
end;
if any(~isfinite(i1)) | any(~isfinite(i2)),
   error('Subscript is not finite.');
end;      
if rc1 | rc2
   if ~isempty(Pc),
      eval('Pc = subsasgn(Pc,S_par,Rc);','error(peelf(lasterr));');
   end;
else,
   [sn1,sn2,sn3] = size(Pc);
   if isnumeric(i1) & isnumeric(i2) & length(i2)==1,
      if i2==1 & sn2==1,  
         if sn1==0,
            error('Improper null index in null matrix.');
         end;         
         Pci = zeros(sn1,sn2);
         nPc = zeros(sn2,sn1-length(i1),sn3);  
         if sn3>0,
            for i = 1:sn3,
               Pci = Pc(:,:,i);
               eval('Pci(i1) = [];', 'error(peelf(lasterr));');
               nPc(:,:,i) = Pci;
            end;
         end;
         Pc = nPc;
      else
         error('Indexed empty matrix assignment is not allowed.');
      end;
   else
      if one_ind & ischar(i1) & strcmp(i1,':'),
         Pc = zeros(0,0,0);
         sn1 = 0; sn2 = 0; ps1 = 0; ps2 = 0;
      else
         if ~isempty(Pc),
            eval('Pc = subsasgn(Pc,S_par,[]);','error(peelf(lasterr));');
         end;
      end;
   end;
end;
[sn1,sn2,sn3]=size(Pc);
if ~isempty(Pc),
   if one_ind &(prod([sn1,sn2,sn3])==prod([ps1,ps2,s3])),
      Pc=reshape(Pc,ps1,ps2,s3); sn1=ps1; sn2=ps2;
   end
   PP.c=Pc; PP.s(1)=sn1; PP.s(2)=sn2; PP.d=sn3-1;
   if ((isempty(var1))|(isempty(var2))),
      PP.v = [var1,var2];
      if (sn3 > 1)&(isempty(PP.v)), PP.v = PGLOBAL.VARIABLE; end;
   else, PP.v = var1;
     if ~strcmp(var1,var2), warning('Inconsistent variables.'); end;
   end;
   if isa(R,'pol'),
      [th,h,PP,R] = testhp(PP,R,PP.v);
      if ~th, warning('Inconsistent sampling periods.');
      end;
      PP.h = h;
   end;
   PP = pclear(PP);
else,
   if sn1==0 | sn2==0,
      PP.c=zeros(sn1,sn2,0);
      PP.s(1)=sn1; PP.s(2)=sn2;
   else,
      PP.c=[]; PP.s=[0 0];
   end;
   PP.d=[]; PP.v=''; PP.h=[];
end;

%end .. @pol/subsasgn