function L = lrand(varargin)
%LRAND Create a left-denominator fraction with random real coefficients
%
% The command
%     L = LRAND(DEGL,I,J)  
% generates a "random"  I-by-J left-denominator fraction L with normally
% distributed coefficients. The degrees of numerator and denominator are DEGL.
% If J is missing then a square I-by-I numerator and denominator are created. 
% If also I is missing then the scalar left-denominator fraction is created.
%
%     L = LRAND(DEGL,I,J,OPTIONS)  
% To generate the left-denominator fraction L with some specified
% attributes, include OPTIONS among the input arguments.
% OPTIONS is one of the strings 'sta', 'bista', 'prop'  or 'ps'. 
%
%     L = LRAND(DEGL,I,J,'sta')  
% The OPTIONS = 'sta' generates L a 'stable' left-denominator fraction 
% (see the macro ISSTABLE for a more information).
%
%     L = LRAND(DEGL,'bista')  
% The OPTIONS = 'bista' generates L a 'bistable' left-denominator fraction 
% (both - the numerator and the denominator are stable).
%
%     L = LRAND(DEGL,I,J,'prop')  
% The OPTIONS = 'prop' generates L a 'proper' left-denominator fraction 
% (see the macro ISPROPER for a more information).
%
%     L = LRAND(DEGL,I,J,'sprop')  
% The OPTIONS = 'sprop' generates L a 'strictly proper' left-denominator fraction 
% (see the macro ISPROPER for a more information).
%
%     L = LRAND(DEGL,I,J,'ps')  or  L = LRAND(DEGL,I,J,'prop','sta') 
% The OPTIONS = 'ps' generates L a 'proper' and 'stable' left-denominator 
% fraction.
%
%     L = LRAND(DEGL,I,J,'sps')  or  L = LRAND(DEGL,I,J,'sprop','sta') 
% The OPTIONS = 'sps' generates L a 'strictly proper' and 'stable' 
% left-denominator fraction.
%
% Any of these syntaxes may be combined with the string 'int' to produce
% "small integer" coefficients and a VARIABLE. VARIABLE is one of the 
% strings 's', 'p', 'd', 'q', 'z', 'z^-1'  or 'zi'. 


%       Author(s):  S. Pejchova  31-08-00
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 19-Oct-2000  S. Pejchova   $
%                         $Date: 19-Nov-2002  S. Pejchova   $  

global PGLOBAL;
eval('PGLOBAL.VARIABLE;', 'painit;');

na = nargin;
if na > 7, error('Too many input arguments.'); end;

i_tr = '';  vx = '';    staL = '';  proL='';  se_db = 0; L=[];

degL = max([1,round(abs(2*randn(1)))]); str_deg=num2str(degL);
rL = max([1,round(4*rand(1))]); 
cL = max([1,round(4*rand(1))]); 
vx=PGLOBAL.VARIABLE;
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
            [rL,cL]=size(argm); degL=[]; str_deg='[]'; 
         else, 
            degL=argm; rL=1; cL=1; str_deg=num2str(degL); 
         end;
         se_db=1;
      elseif (isempty(argm))|(any(isinf(argm))),
            error('The size must be a nonempty and finite number.'); 
      elseif se_db==1,
         rL=argm; cL=rL; se_db=2;
      elseif se_db==2,
         cL=argm; se_db=3;
      end;
   elseif ischar(argm),
      argm=lower(argm);
      switch argm,
       case {'sta','bista'},
         staL=argm;
       case {'prop','sprop'},
         proL=argm;
       case 'ps',
         proL='prop'; staL='sta';
       case 'sps',
         proL='sprop'; staL='sta';
       case {'s','d','p','q','z','z^-1','zi'},
          vx=argm;
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
vvx=[',''',vx,''''];
Nm=[]; Dn=[];
for ii=1:10, 
   if ~isempty(staL)&~isempty(degL)&~isinf(degL),
      eval(['Dn=prand(',num2str(degL*rL),',''sta'',',num2str(rL),i_tr,vvx,');'],...
         'error(peel(lasterr));'); 
   elseif isinf(degL), Dn=eye(rL); 
   else,
      eval(['Dn=prand(',str_deg,',',num2str(rL),i_tr,vvx,');'],...
         'error(peel(lasterr));'); 
   end; 
   Dold=Dn;
   if ~isempty(proL)&~isempty(degL)&~isinf(degL),  Dn=rowred(Dn); end;
   if isempty(degL)|isinf(degL), break; end;
   if (deg(Dn)==degL)&(~issingular(Dn))&(abs(norm(Dn-Dold))<=PGLOBAL.ZEROING), break; end;
end;

if ~isempty(proL)&~isempty(degL)&~isinf(degL),
   d1=deg(Dn,'row');
   if strcmp(proL,'sprop'), 
      d1=d1-1; d1(d1<0)=-inf; 
   end;
   n1=repmat(d1,[1,cL]);
   eval(['Nm=prand(n1',i_tr,vvx,');'],'error(peel(lasterr));'); 
   if strcmp(proL,'sprop'), 
      if strcmp(vx,'d'),   Nm=d*Nm; 
      elseif strcmp(vx,'z^-1')|strcmp(vx,'zi'),  Nm=(z^-1)*Nm;  
      end; 
   end;
elseif strcmp(staL,'bista')&~isempty(degL)&~isinf(degL)&(rL==cL),
   for jj=1:10,
      eval(['Nm=prand(',num2str(degL*min([rL,cL])),',''sta'',',num2str(rL),i_tr,vvx,');'],...
         'error(peel(lasterr));');
      if deg(Nm)<=degL, break; end; 
   end;
else,
   eval(['Nm=prand(',str_deg,',',num2str(rL),',',num2str(cL),i_tr,vvx,');'],...
      'error(peel(lasterr));'); 
end; 
L=ldf(Dn,Nm);  wr=0;
if strcmp(proL,'prop')&~isproper(L), wr=1; 
elseif strcmp(proL,'sprop')&(isproper(L)~=2), wr=1;
elseif  ~isempty(staL)&~isstable(L), wr=1;, 
end;
if strcmp(staL,'bista')&(rL==cL)&~isstable(L.n), wr=1; end;

if wr, warning(' The random fraction has not expected properties. Try once more.');  end;

%end .. lrand
