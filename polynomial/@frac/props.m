function Po = props(PP,varargin)
%PROPS  Display/modify properties of fraction object
%
% The command
%   PROPS(P,VALUE)  
% sets the property of P corresponding to VALUE equal to the value VALUE.
% The command
%   PROPS(P,Value1,Value2,...)  
% sets multiple property values of P equal to the values Value1, Value2,...
% 
% For admissible VALUEs, call  PROPS(P)  and see the display.
% When VALUE is 'cop','ncop','cop?','prop','nprop','prop?',
% it must be followed by the corresponding tolerance,
% e.g.  PROPS(P,'cop',10^-4)  or  PROPS(P,'prop?',[]) .
% Tolerance "[]" means "standard" or "don't care".
%
% The commmand
%   PROPS(P)  
% from the kayboard displays all the properites of P and
% their admissible values.
%
% See also: POL/PROPS, TSP/PROPS.

%       Author(s):  S. Pejchova  16-03-2000
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date:20-Mar-2000 17:00:50  J.Jezek $
%                         $Date: 05-Jun-2000 09:45:00  J.Jezek $
%       $Revision: 3.0 $  $Date: 10-Jul-2000 16:46:11  S.Pejchova $
%                         $Date: 04-Oct-2000   S.Pejchova $
%                         $Date: 30-Jun-2002   J.Jezek    $
%                         $Date: 14-Oct-2002   J.Jezek    $
%                         $Date: 28-Feb-2003   J.Jezek    $

global PGLOBAL;

ni = nargin; no = nargout;
name = inputname(1);
if isempty(name)&ni>1,
   error('1st argument is not a named variable.');
else,
   switch class(PP),
   case 'frac', colxx = ['POLYNOMIAL FRACTION'];
   case 'rdf',colxx = ['RIGHT-DENOMINATOR FRACTION'];
   case 'ldf',colxx = ['LEFT-DENOMINATOR FRACTION'];
   case 'mdf',colxx = ['MATRIX-DENOMINATOR FRACTION'];
   case 'sdf',colxx = ['SCALAR-DENOMINATOR FRACTION'];
   otherwise, error('1st argument is not a polynomial fraction.'); 
   end;
end;

% Now we have set(P,V1,V2, ...)
ni1 = ni-1; i = 1;
while i<=ni1,
   toltoberead = '';
   Value = varargin{i}; 
   if isa(Value,'pol'),
      [vs1,vs2,vd] = size(Value);
      if all([vs1,vs2,vd]==1)&(~any(Value.c(:,:)-[0,1])),
         Value = Value.v;
      end; 
   end;
   if isstr(Value),
      if strcmp(Value,'zi'), Value = 'z^-1'; end;
      I = strmatch(Value,{'s';'p';'z^-1';'d';'z';'q'},'exact');
      if ~isempty(I),
         Oldvalue = PP.v;
         PP.v = Value; PP.num.v = Value; PP.den.v = Value;
         J = strmatch(Value,{'z^-1';'d';'z';'q'},'exact');
         K = strmatch(Oldvalue,{'s';'p'},'exact');
         if ~isempty(J) & ~isempty(K),
            PP.h = 1; PP.num.h = 1; PP.den.h = 1;
         end;
         J = strmatch(Value,{'s';'p'},'exact');
         if ~isempty(J),
            PP.h = 0; PP.num.h = 0; PP.den.h = 0;
         end;
         PP.c = 'cop?'; PP.r = 'red?'; PP.p = 'prop?';
         PP.tc = []; PP.tp = [];
      elseif strcmp(Value,'red')|strcmp(Value,'reduce')|strcmp(Value,'reduced'),
         PP.r = 'red';
      elseif strcmp(Value,'nred')|strcmp(Value,'nreduce')|...
            strcmp(Value,'nreduced'),
         PP.r = 'nred';
      elseif strcmp(Value,'red?')|strcmp(Value,'reduce?')|...
            strcmp(Value,'reduced?'),
         PP.r = 'red?';
      elseif strcmp(Value,'cop')|strcmp(Value,'coprime'),
         PP.c = 'cop'; toltoberead = 'c';
      elseif strcmp(Value,'ncop')|strcmp(Value,'ncoprime'),
         PP.c = 'ncop'; toltoberead = 'c';
      elseif strcmp(Value,'cop?')|strcmp(Value,'coprime?'),
         PP.c = 'cop?'; toltoberead = 'ignore'; PP.tc = [];
      elseif strcmp(Value,'prop')|strcmp(Value,'proper'),
         PP.p = 'prop'; toltoberead = 'p';
      elseif strcmp(Value,'nprop')|strcmp(Value,'nproper'),
         PP.p = 'nprop'; toltoberead = 'p';
      elseif strcmp(Value,'prop?')|strcmp(Value,'proper?'),
         PP.p = 'prop?'; toltoberead = 'ignore'; PP.tp = [];
      else,
         PP.u = Value;
      end;
      if ~isempty(toltoberead);
         i = i+1;
         if i>ni1,
            error('Tolerance for ''cop'' or ''prop'' is missing.');
         end;
         tol = varargin{i};
         if isempty(tol),
            tol = PGLOBAL.ZEROING;
         elseif ~isa(tol,'double') | length(tol)~=1 | ...
               ~isreal(tol) | tol<0 | tol>1,
            error('Tolerance for ''cop'' or ''prop'' is invalid or missing,');
         end;
         if strcmp(toltoberead,'c'),
            PP.tc = tol;
         elseif strcmp(toltoberead,'p'),
            PP.tp = tol;
         end;
      end;
   elseif isa(Value,'double') & (isempty(Value) | ...
         (length(Value)==1 & isreal(Value) & isfinite(Value) & Value>=0)),
      var = PP.v;
      if isempty(var), Value = [];
      elseif strcmp(var,'s') | strcmp(var,'p'), Value = 0;
      end;
      PP.h = Value; PP.num.h = Value; PP.den.h = Value;
   else, PP.u = Value;
   end; % isstr(Va...
   i = i+1;
end % while
if ni>1,
   if isempty(PP.num.v)&(isempty(PP.den.v)), 
      PP.v = ''; 
   end; 
   % Assign PP in caller's workspace
   assignin('caller',name,PP);
end;

% Return Props and set nothing if nargin==1
colA = char('size','numerator degree','denominator degree    ');
Pnd = PP.num.d; Pdd = PP.den.d;
dgB1 = int2str(Pnd);
if isempty(Pnd), dgB1 = '[]'; 
elseif isinf(Pnd), dgB1 = '-Inf';
end;
dgB2 = int2str(Pdd);
if isempty(Pdd), dgB1 = '[]'; 
elseif isinf(Pdd), dgB1 = '-Inf';
end;
crpB3 = PP.c; crpB4 = PP.r; crpB5 = PP.p;
crtc = num2str(PP.tc); crtp = num2str(PP.tp);
if isempty(crtc), crtc = '[]';
end;
if isempty(crtp), crtp = '[]';
end;
colB = char([int2str(PP.s(1)),'-by-',int2str(PP.s(2))],dgB1,dgB2);

P_out=char(' ',colxx,[colA,colB],' ');
h1 = cl2str(PP.h);
usdt = cl2str(PP.u);
col1 = char('PROPERTY NAME:','variable symbol   ','sampling period',...
   'coprime flag','reduce flag','proper flag','coprime tolerance',...
   'proper tolerance','user data');
col2 = char('CURRENT VALUE:   ',PP.v,h1,crpB3,crpB4,crpB5,crtc,crtp,usdt);
col3 = char('  AVAILABLE VALUES:',[' ''s'',''p'',''z^-1'',''d'',''z'',''q'''],...
      '  0 for c-time, nonneg or [] for d-time',...
      ' ''cop'', ''ncop'', ''cop?'' ',...
      ' ''red'', ''nred'', ''red?'' ',...
      ' ''prop'', ''nprop'', ''prop?'' ',...
      '  any number from 0 to 1, or []',...
      '  any number from 0 to 1, or []',...      
      '  arbitrary');
P_out=char(P_out,[col1,col2,col3],' ');
if no==0 & length(dbstack)==1,
   disp(P_out);
else,
   Po=P_out; 
end;

%end .. @frac/props
