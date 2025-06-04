function pprop(P,varargin)
%PPROP  Sets the properties of the polynomial object P
%
% The commmand
%   PPROP(P,VALUE)  
% sets the property of P corresponding to VALUE equal to the value 
% VALUE. The commmand
%   PPROP(P,Value1,Value2,...)  
% sets multiple property values of P equal to the values Value1, 
% Value2,... The commmand
%   PPROP(P)  
% displays all the properies of P and their admissible values.

%       Author(s):  S. Pejchova, M. Sebek 5-3-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 01-Feb-1999 16:44:50   $

ni = nargin; P_cLAss=0;
if ni < 1,
   error('Not enough input arguments.');
elseif isa(P,'pol'),
   P_cLAss=1;  
end;

[sP1 ,sP2,dP]= size(P);
if ni==1,
  % Return Props and set nothing if nargin==1
  disp(' ');
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
  if P_cLAss,
     disp(colxx);
     disp([char('size       ',colA),...
           char([' ',int2str(sP1),'-by-',int2str(sP2)],colB)]);
     if isinf(dP)|dP>0, disp(' '); end;
     col1=char('PROPERTY NAME:','variable symbol   ','user data');
     usdt=P.u;
  else, usdt=P;
  end;
  if ~isstr(usdt),
    if isa(usdt,'double') & ndims(usdt)==2 & (isempty(usdt) | ... 
          (size(usdt,1)<=1 & size(usdt,2)<25)),
        if isempty(usdt) & ~isequal(size(usdt),[0 0]),
         usdt = sprintf('[%dx%d] double',size(usdt,1),size(usdt,2));
        else
         usdt = mat2str(usdt,3);
        end
    elseif isa(usdt,'cell') & isempty(usdt),
      if isequal(size(usdt),[0 0]),
         usdt = '{}';
      else
         usdt = sprintf('{%dx%d} cell',size(usdt,1),size(usdt,2));
      end
    else
      % Too big to be displayed
      us_str = mat2str(size(usdt));
      us_str = strrep(us_str(2:end-1),' ','x');
      if isa(usdt,'cell'),
         us_str = ['{',us_str,'} ',class(usdt)];
      else
         us_str = ['[',us_str,'] ',class(usdt)];
      end
      usdt=us_str;
    end
  elseif (length(usdt) > 25) | ~P_cLAss,
     usdt=sprintf('string of length %d',length(usdt));
  end;
  if P_cLAss,
     col2=char('CURRENT VALUE:   ',P.v,usdt);
     col3=char('  AVAILABLE VALUES:',[' ''s'',''p'',''z^-1'',''d'',''z'',''q'''],...
               '  arbitrary');
     disp([col1,col2,col3]);
  else,
     disp('The input is not a polynomial object.');
     disp(['It is a ',usdt,'.']);
  end;             
  disp(' ');
  return
end

% Now we have set(P,V1,V2, ...)
name = inputname(1);
if ~P_cLAss,
   error('1st argument is not a polynomial object.');
elseif isempty(name),
   error('1st argument must be a named variable.')
end;
for i=1:ni-1,
   Value = varargin{i}; [vs1,vs2,vd]=size(Value);
   if isa(Value,'pol')& all([vs1,vs2,vd]==1),
      Vcoef=Value.c;
      if ~any(Vcoef(:,:)-[0,1]), P.v=Value.v;
      else, P.u=Value;
      end;
   elseif isstr(Value),
      if strcmp(Value,'zi'), Value='z^-1'; end;
      I=strmatch(Value,{'s';'p';'z^-1';'d';'z';'q'},'exact');
      if ~isempty(I), P.v = Value;
      else, P.u = Value;
      end;
   else, P.u = Value;
   end;
end % for
if isempty(dP)|dP<=0, P.v=''; end;

% Assign P in caller's workspace
assignin('caller',name,P)

%end .. pprop

