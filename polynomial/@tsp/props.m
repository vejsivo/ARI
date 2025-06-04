function T_o=props(T,varargin)
%PROPS  Display/modify properties of two-sided polynomial object
%
% The commmand
%   PROPS(T,VALUE)  
% sets the property 'user data' of T  equal to the value  VALUE.
% The commmand
%   PROPS(T)  
% from the keyboard displays all the properties of T and
% their admissible values.
%
% See also POL/PROPS, FRAC/PROPS.

%       Author(s):  S. Pejchova  21-07-99
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 16-Mar-2000 14:25:50  $
%                         $Date: 22-May-2000  J.Jezek  $
%                         $Date: 11-Jul-2000 11:53:11  S.Pejchova $
%                         $Date: 30-Jun-2002 J.Jezek   $
%                         $Date: 28-Feb-2002 J.Jezek   $

ni = nargin; no = nargout;
name = inputname(1);
if isempty(name)&ni>1,
   error('Invalid 1st argument; not a named variable.');
elseif ~isa(T,'tsp'),
   error('Invalid 1st argument.');
end;
% Now we have set(T,V1,V2, ...)
for ii=1:ni-1,
   Value = varargin{ii}; 
   if isa(Value,'pol'),
      [vs1,vs2,vd]=size(Value);
      if all([vs1,vs2,vd]==1)&(~any(Value.c(:,:)-[0,1])), 
         Value = Value.v;
      end; 
   end;
   if isstr(Value),
      if strcmp(Value,'z') | strcmp(Value,'zi') | strcmp(Value,'z^-1'), 
      else, T.u = Value;
      end;
   elseif isa(Value,'double') & (isempty(Value) | ...
         (length(Value)==1 & isreal(Value) & isfinite(Value) & Value>=0)),
      T.h = Value;
   else, T.u = Value;
   end;
end % for
if ni>1,
   if (isempty(T.d)|(~T.o & T.r<=0)),  
      T.v = ''; 
   end;
   if isempty(T.v), T.h = []; end;
   % Assign T in caller's workspace
   assignin('caller',name,T); 
end;

% Return Props and set nothing if nargin==1
colA=char('degree','trailing degree    ');
sgB1=' ';  sgB2=' ';
[sT1 ,sT2] = size(T); Tr = T.r; To=T.o; degT=T.d;  tdegT=T.t;
if degT>0, sgB1='+'; elseif degT<0, sgB1='-'; end;
if tdegT>0, sgB2='+'; elseif tdegT<0, sgB2='-'; end;
colB=char([sgB1,int2str(abs(degT))],[sgB2,int2str(abs(tdegT))]);
if isempty(Tr),
   colxx=['EMPTY TWO-SIDED POLYNOMIAL MATRIX']; colA=' '; colB=' ';
elseif isinf(Tr),
   colxx=['ZERO TWO-SIDED POLYNOMIAL MATRIX']; colB=char('-Inf','+Inf');
elseif (~To & Tr<=0),
   colxx=['CONSTANT TWO-SIDED POLYNOMIAL MATRIX'];
   colA=' '; colB=' ';
else,
   colxx=['TWO-SIDED POLYNOMIAL MATRIX'];
end;

T_out=[char('size       ',colA),...
      char([' ',int2str(sT1),'-by-',int2str(sT2)],colB)];
T_out=char(' ',colxx,T_out);
if isinf(Tr)|Tr>0,, T_out=char(T_out,' '); end;
col1=char('PROPERTY NAME:','variable symbol   ','sampling period','user data');
h1=cl2str(T.h); usdt=cl2str(T.u);
col2=char('CURRENT VALUE:   ',T.v,h1,usdt);
col3=char('  AVAILABLE VALUES:',['  ''z'''],'  nonnegative or []','  arbitrary');
T_out=char(T_out,[col1,col2,col3],' ');
if no==0 & length(dbstack)==1,
   disp(T_out);
else,
   T_o=T_out; 
end 

%end .. @tsp/props
