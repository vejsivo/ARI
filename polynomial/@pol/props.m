function Po = props(P,varargin)
%PROPS  Display/modify properties of polynomial object
%
% The commmand
%   PROPS(P,VALUE)  
% sets the property of P corresponding to VALUE equal to the value 
% VALUE. The commmand
%   PROPS(P,Value1,Value2,...)  
% sets multiple property values of P equal to the values Value1, 
% Value2,... The command
%   PROPS(P)  
% displays all the properties of P and their admissible values.
%
% See also: TSP/PROPS, FRAC/PROPS.

%       Author(s):  S. Pejchova, M. Sebek 5-3-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 16-Mar-2000 13:57:50   $
%       $Revision: 3.0 $  $Date: 06-Jun-2000 10:15:00  J.Jezek  $        
%       $Revision: 3.0 $  $Date: 11-Jul-2000 10:59:11  S.Pejchova $
%                         $Date: 28-Feb-2003           J.Jezek  $

ni = nargin;

name = inputname(1);
if isempty(name)&ni>1,
   error('1st argument is not a named variable.');
elseif ~isa(P,'pol'),
   error('1st argument is a nonpolynomial object.');
end;
% Now we have set(P,V1,V2, ...)
for ii=1:ni-1,
   Value = varargin{ii}; 
   if isa(Value,'pol'),
      [vs1,vs2,vd]=size(Value);
      if all([vs1,vs2,vd]==1)&(~any(Value.c(:,:)-[0,1])), Value=Value.v; end; 
   end;
   if isstr(Value),
      if strcmp(Value,'zi'), Value='z^-1'; end;
      I=strmatch(Value,{'s';'p';'z^-1';'d';'z';'q'},'exact');
      if ~isempty(I),
         Oldvalue = P.v; P.v = Value;
         J = strmatch(Value,{'z^-1';'d';'z';'q'},'exact');
         if ~isempty(J) & (strcmp(Oldvalue,'s') | strcmp(Oldvalue,'p')),
            P.h = 1;
         end;
      else, P.u = Value;
      end;
   elseif isa(Value,'double') & (isempty(Value) | ...
         (length(Value)==1 & isreal(Value) & isfinite(Value) & Value>=0)),
      P.h = Value;
   else, P.u = Value;
   end;
end % for
if ni>1,
   if isempty(P.d)|(P.d<=0), P.v='';  end; 
   if isempty(P.v), P.h = []; 
   elseif strcmp(P.v,'s') | strcmp(P.v,'p'), P.h = 0;   
   end;
   % Assign P in caller's workspace
   assignin('caller',name,P); 
end;

% Return Props and set nothing if nargin==1
[sP1 ,sP2,dP]= size(P);
colA='degree';
colB=[' ',int2str(dP)];
if isempty(dP),
   colxx=['EMPTY POLYNOMIAL MATRIX']; colA=' ';
elseif isinf(dP),
   colxx=['ZERO POLYNOMIAL MATRIX']; colB='-Inf';
elseif dP==0,
   colxx=['CONSTANT POLYNOMIAL MATRIX'];
   colA=' '; colB=' ';
else,
   colxx=['POLYNOMIAL MATRIX'];
end;
P_out=[char('size       ',colA),...
      char([' ',int2str(sP1),'-by-',int2str(sP2)],colB)];
P_out=char(' ',colxx,P_out);
if isinf(dP)|dP>0, P_out=char(P_out,' '); end;

col1=char('PROPERTY NAME:','variable symbol   ','sampling period','user data');
h1=cl2str(P.h); usdt=cl2str(P.u);
col2=char('CURRENT VALUE:   ',P.v,h1,usdt);
col3=char('  AVAILABLE VALUES:',[' ''s'',''p'',''z^-1'',''d'',''z'',''q'''],...
     '  0 for c-time, nonneg or [] for d-time','  arbitrary');
P_out=char(P_out,[col1,col2,col3],' ');
if ni==1,
   disp(P_out);
else,
   Po=P_out; 
end;
   

%end .. @pol/props
