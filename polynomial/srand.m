function S = srand(varargin)
%SRAND  Create a scalar-denominator fraction with random real coefficients
%
% The command
%     S = SRAND(DEGS,I,J)  
% generates a "random"  I-by-J scalar-denominator fraction S with normally
% distributed coefficients. The degrees of numerators and denominator are DEGS.
% If J is missing then a square I-by-I numerator is created. 
% If also I is missing then the scalar fraction is created.
%
%     S = SRAND(DEGS,I,J,OPTIONS)  
% To generate the scalar-denominator fraction S with some specified
% attributes, include OPTIONS among the input arguments.
% OPTIONS is one of the strings 'STA', 'BISTA', 'PROP'  or 'PS'. 
%
%     S = SRAND(DEGS,I,J,'STA')  
% The OPTIONS = 'STA' generates S a 'stable' scalar-denominator fraction 
% (see the macro ISSTABLE for a more information).
%
%     S = SRAND(DEGS,I,J,'BISTA')  
% The OPTIONS = 'BISTA' generates S a 'bistable' scalar-denominator fraction 
% (both - the numerator and the denominator are stable).
%
%     S = SRAND(DEGS,I,J,'PROP')  
% The OPTIONS = 'PROP' generates S a 'proper' scalar-denominator fraction 
% (see the macro ISPROPER for a more information).
%
%     S = SRAND(DEGS,I,J,'PS')  or  S = SRAND(DEGS,I,J,'PROP','STA')  
% The OPTIONS = 'PS' generates S a 'proper' and 'stable' scalar-denominator 
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

i_tr = '';  staS = '';  proS='';  se_db = 0; S=[]; str_sta='';

degS = max([1,round(abs(2*randn(1)))]); str_deg=num2str(degS);
rS = max([1,round(4*rand(1))]); 
cS = max([1,round(4*rand(1))]); 
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
            [rS,cS]=size(argm); degS=[]; str_deg='[]'; 
         else, 
            degS=argm; rS=1; cS=1; str_deg=num2str(degS); 
         end;
         se_db=1;
      elseif (isempty(argm))|(any(isinf(argm))),
            error('The size must be a nonempty and finite number.'); 
      elseif se_db==1,
         rS=argm; cS=rS; se_db=2;
      elseif se_db==2,
         cS=argm; se_db=3;
      end;
   elseif ischar(argm),
      switch argm,
       case {'sta','bista'},
         staS=argm; str_sta=[',''sta'''];
       case 'prop',
         proS=argm;
       case 'ps',
         proS='prop'; staS='sta';
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
Nm=pol(zeros(rS,cS)) ; Dn=pol(1);
if isempty(degS)|(~isempty(degS)&~isinf(degS)),
   for ii=1:10, 
      eval(['Dn=prand(',str_deg,str_sta,',1',i_tr,vvx,');'],... 
         'error(peel(lasterr));');
      if isempty(degS), break; end;
      if deg(Dn)==degS, break; end; 
   end; 
end; 
Dnc=Dn.c;
if size(Dnc,3), 
   Dn = Dn + ~sum(abs(Dnc),3); 
end;

if (isempty(degS)|(rS~=cS)|(~isempty(degS)&isfinite(degS)))&~strcmp(staS,'bista'), 
   for jj=1:10, 
      eval(['Nm=prand(',str_deg,',',num2str(rS),',',num2str(cS),i_tr,vvx,');'],... 
         'error(peel(lasterr));');  
      if isempty(degS)| isempty(proS)| all(all(deg(Nm,'ent')==degS)), break; end; 
   end;
elseif ~isinf(degS), 
   for jj=1:10, 
      eval(['Nm=prand(',num2str(degS*rS),',''sta'',',num2str(rS),i_tr,vvx,');'],... 
         'error(peel(lasterr));'); 
      if deg(Nm)<=degS, break; end; 
   end; 
end; 

S=sdf(Nm,Dn);  wr=0;
if ~isempty(proS)&~isproper(S), wr=1; 
elseif  ~isempty(staS)&~isstable(S), wr=1;, 
end;
if strcmp(staS,'bista')&(rS==cS)&~isstable(S.n), wr=1; end;

if wr, warning(' The random fraction has not expected properties. Try once more.');  end;

%end .. srand
