function M = mrand(varargin)
%MRAND  Create a matrix-denominator fraction with random real coefficients
%
% The command
%     M = MRAND(DEGM,I,J)  
% generates a "random"  I-by-J matrix-denominator fraction M with normally
% distributed coefficients. The degrees of numerators and denominators are DEGM.
% If J is missing then a square I-by-I numerator and denominator are created. 
% If also I is missing then the scalar matrix-denominator fraction is created.
%
%     M = MRAND(DEGM,I,J,OPTIONS)  
% To generate the matrix-denominator fraction M with some specified
% attributes, include OPTIONS among the input arguments.
% OPTIONS is one of the strings 'STA', 'BISTA', 'PROP'  or 'PS'. 
%
%     M = MRAND(DEGM,I,J,'STA')  
% The OPTIONS = 'STA' generates M a 'stable' matrix-denominator fraction 
% (see the macro ISSTABLE for a more information).
%
%     M = MRAND(DEGM,I,J,'BISTA')  
% The OPTIONS = 'BISTA' generates M a 'bistable' matrix-denominator fraction 
% (both - the numerator and the denominator are stable).
%
%     M = MRAND(DEGM,I,J,'PROP')  
% The OPTIONS = 'PROP' generates M a 'proper' matrix-denominator fraction 
% (see the macro ISPROPER for a more information).
%
%     M = MRAND(DEGM,I,J,'PS')  or  M = MRAND(DEGM,I,J,'PROP','STA')  
% The OPTIONS = 'PS' generates M a 'proper' and 'stable' matrix-denominator 
% fraction.
%
% Any of these syntaxes may be combined with the string 'INT' to produce
% "small integer" coefficients and a VARIABLE. VARIABLE is one of the 
% strings 'S', 'P', 'D', 'Q', 'Z', 'Z^-1'  or 'ZI'. 


%       Author(s):  S. Pejchova  31-08-00
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 19-Oct-2000  S. Pejchova   $

global PGLOBAL;
eval('PGLOBAL.VARIABLE;', 'painit;');
na = nargin;
if na > 7, error('Too many input arguments.'); end;

i_tr = '';  vvx = '';    staM = '';  proM='';  se_db = 0; M=[];

degM = max([1,round(abs(2*randn(1)))]); str_deg=num2str(degM);
rM = max([1,round(4*rand(1))]); 
cM = max([1,round(4*rand(1))]); 
vvx=PGLOBAL.VARIABLE;
if na>0,
 for kk=1:na
   argm=varargin{kk};
   if isa(argm,'double')&(ndims(argm)==2)&(~any(size(argm)>1)),
      argcol=argm(:); argcol((argcol<0)&(isinf(argcol)))=0;
      if any(~isfinite(argcol)),
         error('Nan and +Inf not allowed for degree and size.');
      elseif (se_db<2)&((norm(argcol-round(argcol))>PGLOBAL.ZEROING)|(any(argcol<0))),
         error('The degree and size must be nonnegative integers.');
      end;
      if ~se_db,
         if isempty(argm), 
            [rM,cM]=size(argm); degM=[]; str_deg='[]'; 
         else, 
            degM=argm; rM=1; cM=1; str_deg=num2str(degM); 
         end;
         se_db=1;
      elseif (isempty(argm))|(any(isinf(argm))),
            error('The size must be a nonempty and finite number.'); 
      elseif se_db==1,
         rM=argm; cM=rM; se_db=2;
      elseif se_db==2,
         cM=argm; se_db=3;
      end;
   elseif ischar(argm),
      switch argm,
       case {'sta','bista'},
         staM=argm;
       case 'prop',
         proM=argm;
       case 'ps',
         proM='prop'; staM='sta';
       case {'s','d','p','q','z','z^-1','zi'},
          vvx=argm;
       case 'int', 
          i_tr=[',''int'''];   
       otherwise,
         error('Invalid command option.');
      end;
   else,
      error(['Invalid ',nth(kk),' argument.']);
   end;
 end;  
end;
vvx=[',''',vvx,''''];
Nm=pol(zeros(rM,cM)); Dn=pol(ones(rM,cM));

if ~isempty(staM)&~isempty(degM)&~isinf(degM),
   for mm=1:rM, for nn=1:cM,
         eval(['Dn(mm,nn)=prand(',str_deg,',''sta'',1',i_tr,vvx,');'],... 
            'error(peel(lasterr));'); 
   end;  end;
elseif ~isinf(degM), 
   eval(['Dn=prand(',str_deg,',',num2str(rM),',',num2str(cM),i_tr,vvx,');'],... 
      'error(peel(lasterr));'); 
end;
Dnc=Dn.c;
if size(Dnc,3), 
   Dn = Dn + ~sum(abs(Dnc),3); 
end;

if isempty(degM)|(~isempty(degM) & ~isinf(degM)),
   if ~isempty(proM)&~isempty(degM),
      d1=deg(Dn,'ent'); 
      eval(['Nm=prand(d1',i_tr,vvx,');'],'error(peel(lasterr));'); 
   elseif strcmp(staM,'bista')&~isempty(degM), 
      for mm=1:rM, for nn=1:cM, 
            eval(['Nm(mm,nn)=prand(',str_deg,',''sta'',1',i_tr,vvx,');'],... 
               'error(peel(lasterr));'); 
      end;  end;
   else, 
      eval(['Nm=prand(',str_deg,',',num2str(rM),',',num2str(cM),i_tr,vvx,');'],...
         'error(peel(lasterr));'); 
   end; 
end;
M=mdf(Nm,Dn);  wr=0;

if ~isempty(proM)&~isproper(M), wr=1; 
elseif  ~isempty(staM)&~isstable(M), wr=1;, 
end;

if strcmp(staM,'bista'),
   stab1=zeros(rM,cM);
   for mm=1:rM, for nn=1:cM, 
         stab1(mm,nn) = isstable(Nm(mm,nn)); %POL/ISSTABLE on numerator
   end;  end;
   wr = wr|(~all(all(stab1)));
end;   
if wr, warning(' The random fraction has not expected properties. Try once more.');  end;

%end .. mrand
