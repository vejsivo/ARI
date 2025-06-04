function R = rrand(varargin)
%RRAND Create a right-denominator fraction with random real coefficients
%
% The command
%     R = RRAND(DEGR,I,J)  
% generates a "random"  I-by-J right-denominator fraction R with normally
% distributed coefficients. The degrees of numerator and denominator are DEGR.
% If J is missing then a square I-by-I numerator and denominator are created. 
% If also I is missing then the scalar right-denominator fraction is created.
%
%     R = RRAND(DEGR,I,J,OPTIONS)  
% To generate the right-denominator fraction R with some specified
% attributes, include OPTIONS among the input arguments.
% OPTIONS is one of the strings 'STA', 'BISTA', 'PROP'  or 'PS'. 
%
%     R = RRAND(DEGR,I,J,'STA')  
% The OPTIONS = 'STA' generates R a 'stable' right-denominator fraction 
% (see the macro ISSTABLE for a more information).
%
%     R = RRAND(DEGR,I,J,'BISTA')  
% The OPTIONS = 'BISTA' generates R a 'bistable' right-denominator fraction 
% (both - the numerator and the denominator are stable).
%
%     R = RRAND(DEGR,I,J,'PROP')  
% The OPTIONS = 'PROP' generates R a 'proper' right-denominator fraction 
% (see the macro ISPROPER for a more information).
%
%     R = RRAND(DEGR,I,J,'PS')  or  R = RRAND(DEGR,I,J,'PROP','STA')  
% The OPTIONS = 'PS' generates R a 'proper' and 'stable' right-denominator 
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

i_tr = '';  staR = '';  proR='';  se_db = 0; R=[];

degR = max([1,round(abs(2*randn(1)))]); str_deg=num2str(degR);
rR = max([1,round(4*rand(1))]); 
cR = max([1,round(4*rand(1))]); 
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
            [rR,cR]=size(argm); degR=[]; str_deg='[]'; 
         else, 
            degR=argm; rR=1; cR=1; str_deg=num2str(degR); 
         end;
         se_db=1;
      elseif (isempty(argm))|(any(isinf(argm))),
            error('The size must be a nonempty and finite number.'); 
      elseif se_db==1,
         rR=argm; cR=rR; se_db=2;
      elseif se_db==2,
         cR=argm; se_db=3;
      end;
   elseif ischar(argm),
      switch argm,
       case {'sta','bista'},
         staR=argm;
       case 'prop',
         proR=argm;
       case 'ps',
         proR='prop'; staR='sta';
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
Nm=[]; Dn=[];
for ii=1:10, 
   if ~isempty(staR)&~isempty(degR)&~isinf(degR),
      eval(['Dn=prand(',num2str(degR*cR),',''sta'',',num2str(cR),i_tr,vvx,');'],...
         'error(peel(lasterr));'); 
   elseif isinf(degR), Dn=eye(cR); 
   else,
      eval(['Dn=prand(',str_deg,',',num2str(cR),i_tr,vvx,');'],...
         'error(peel(lasterr));'); 
   end; 
   if ~isempty(proR)&~isempty(degR)&~isinf(degR),  Dn=colred(Dn); end;
   if isempty(degR)|isinf(degR), break; end;
   if (deg(Dn)==degR)&(~issingular(Dn)), break; end;
end;

if ~isempty(proR)&~isempty(degR)&~isinf(degR),
   d1=deg(Dn,'col');
   n1=repmat(d1,[rR,1]);
   eval(['Nm=prand(n1',i_tr,vvx,');'],'error(peel(lasterr));'); 
elseif strcmp(staR,'bista')&~isempty(degR)&~isinf(degR)&(rR==cR),
   for jj=1:10,
      eval(['Nm=prand(',num2str(degR*min([rR,cR])),',''sta'',',num2str(cR),i_tr,vvx,');'],...
         'error(peel(lasterr));');
      if deg(Nm)<=degR, break; end; 
   end;
else,
   eval(['Nm=prand(',str_deg,',',num2str(rR),',',num2str(cR),i_tr,vvx,');'],...
      'error(peel(lasterr));'); 
end; 
R=rdf(Nm,Dn);  wr=0;
if ~isempty(proR)&~isproper(R), wr=1; 
elseif  ~isempty(staR)&~isstable(R), wr=1;, 
end;
if strcmp(staR,'bista')&(rR==cR)&~isstable(R.n), wr=1; end;

if wr, warning(' The random fraction has not expected properties. Try once more.');  end;

%end .. rrand
