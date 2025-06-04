function T = trand(varargin)
%TRAND Create a two-sided polynomial matrix with random real coefficients
%
% The command
%     T = TRAND(DEGT,I,J)  
% generates a "random"  I-by-J two-sided polynomial matrix T of trailing
% degree 0 leading degree DEGT with normally distributed coefficients. 
% If J is missing then a square I-by-I matrix is created. If also I 
% is missing then the scalar two-sided polynomial is created.
%
% If DEGT is an integer two-element vector  [DEGT1 DEGT2]  then
% DEGT1 and DEGT2 specify the trailing and the leading degrees of T. 
%
%     T = TRAND(DEGT,I,J,'ENT')  
% To generate the degrees of the entries of T randomly (within {DEGT(1),DEGT(2)}), 
% include string 'ENT' among the input arguments.
%
% Any of these syntaxes may be combined with the string 'INT' to produce
% "small integer" coefficients.

%       Author(s):  S. Pejchova  31-08-00
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 26-Sep-2000  S. Pejchova   $

global PGLOBAL;
na = nargin;
if na > 6, error('Too many input arguments.'); end;

i_tr = '';  typeT = '';  se_db = 0; T=[];

str_deg = num2str(max([1,round(abs(10*rand(1)))]));
o_s = round(2*randn(1));
rT = max([1,round(4*rand(1))]); 
cT = max([1,round(4*rand(1))]); 
eval('PGLOBAL.VARIABLE;', 'painit;');
if na>0,
 for kk=1:na
   argm=varargin{kk};
   if isa(argm,'double')&(ndims(argm)==2)&(~all(size(argm)>1))&length(argm)<3,
      if any(isnan(argm(:))),
         error('NaN is not allowed for degree or size.');
      end;
      if ~se_db,
         o_s =0;
         if isempty(argm), [rT,cT]=size(argm); str_deg='[]';
         else,
            argcol=argm; argcol(isinf(argcol))=0;
            if norm(argcol-round(argcol))>PGLOBAL.ZEROING,
               error('The degree must be an integer.');
            end;
            if all(isinf(argm)), degT=-inf; 
            elseif length(argm)==1, degT=[argm,0]; 
            else,
               degT=argm;
               degT(find(isinf(degT)))=0;
            end;
            degT=sort(degT);
            if ~isinf(degT(1)),  o_s=degT(1);  degT=degT(2)-degT(1);  end;
            rT=1; cT=1;  str_deg=num2str(degT);
         end; 
         se_db=1; 
      elseif any(size(argm)~=1)|isinf(argm)| (norm(argm-round(argm))>PGLOBAL.ZEROING), 
         error('The size must be a nonempty and finite integer.');
      elseif se_db==1,
         rT=argm; cT=rT; se_db=2;
      elseif se_db==2,
         cT=argm; se_db=3;
      end;
   elseif ischar(argm),
      switch argm,
       case 'ent',
         typeT=[',''ent'''];
       case {'z','z^-1','zi'},
       case 'int', i_tr=[',''int'''];;
       otherwise,
         error('Invalid command option.');
      end;
   else,
      error(['Invalid ',nth(kk),' argument.']);
   end;
 end;
end;
eval(['T=prand(',str_deg,',',num2str(rT),',',num2str(cT),typeT,i_tr,',''z'');'],...
   'error(peel(lasterr));');
T=tsp(T);
if o_s,  T=shift(T,o_s);  end;

%end .. trand
