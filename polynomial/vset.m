function V = vset(varargin)
%VSET  Value set of parametric polynomial for given parameters and frequencies.
%
% The commmand
%
%    V=VSET(Q1,...,QM,EXPRESSION-STRING,P0,P1,...PN [,OMEGA [,Q_TYPE]])
%
% computes values of the parametric polynomial family given by expression
%
% EXPRESSION-STRING = 'P0(s) + Expr1*P1(s) + Expr2*P2(s) + ... + ExprN*PN(s)'
%            with paramters Q1, Q2, ... QM  in generalized frequencies OMEGA
%            and where Expr1, Expr2, ...,ExprN are any expressions composed 
%            from the variable names used for Q1, ..., QM and those are 
%            consistent with Matlab syntax, such as
%            EXPRESSION-STRING='P0(s) + Q1*P1(s) + Q2*P2(s) + ...' or
%            EXPRESSION-STRING='P0(s) + (Q1*Q2-Q3)*P1(s)' etc. 
%
% Q1, Q2, ..., QM are vectors of (real or complex) parameters, existing in 
%                 the current workspace.
% P0, P1, ..., PN are polynomials existing in the current workspace.
% Please, note that the first M input arguments must be the parameter
% vector names rather than the vector values and also the next N+1
% input arguments (following EXPRESSION-STRING) must be the polynomial
% matrix names rather than the polynomial values. 
% 
% OMEGA is a vector of generalized frequencies.
%
% The last input argument Q_TYPE can be used to choose the type
% of construction the evaluate points from vectors Q1, ..., QM.
%   Q_TYPE='r' - DEFAULT, the rectangular M-dimensional grid 
%                defined by the vectors Q1,...,QM. It is  constructed
%                from all values combinations of vectors Q1,...,QM.
%                Number of points is NP=length(Q1)*length(Q2)*...*length(QM)
%   Q_TYPE='e' - exactly defined points of values; length of all
%                Q1,...,QM must be the same.
%                Number of points is NP = length(Q1)
%
% V is a matrix of complex numbers with values at particular frequecies
% organized columnwise. 
% The matrix has  NP*prod(size(P0)) rows and length(OMEGA) columns.  
% It can be ploted by macro VSETPLOT.

%       Author(s):  S. Pejchova 30-11-99
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 05-May-2000 14:42:34   $
%       $Revision: 3.0 $  $Date: 02-Nov-2000  S. Pejchova   $

n_I=nargin; eX_Str=''; n_Q=0; n_P=0; oMegA_val=0; P_iY=[];  q_tYPe='r';
N_iY_1=0;

global PGLOBAL;
eval('PGLOBAL.VERBOSE;', 'painit;');
if ~n_I, error('Not enough input arguments.'); end;

last_one=varargin{n_I};
if ischar(last_one),
   if (strcmp(last_one,'r')|strcmp(last_one,'e')),
      q_tYPe=last_one;  varargin(n_I)=[]; n_I=n_I-1;
   else,
      error('Inappropriate string argument for a type of grid construction.');
   end;
 end;
last_one=varargin{n_I};
if ~ischar(last_one)& ~isa(last_one,'pol'), 
   oMegA_val=last_one;  varargin(n_I)=[]; 
   n_I=n_I-1; end;
if isempty(oMegA_val)|(~isa(oMegA_val,'double')),
   error('The set of evaluation frequencies must be a nonempty vector.');
end;
oMegA_val=oMegA_val(:)';
C_iY_n={'.\';'.\(';'./';'./(';'.*';'.*('};

fl1=find(cellfun('isclass',varargin,'char'));
if length(fl1)~=1, error('Expression string is missing or too many strings.'); end;

eX_Str=varargin{fl1};
eX_Str=strrep(eX_Str,' ','');, eX_Str=strrep(eX_Str,'*','.*');
eX_Str=strrep(eX_Str,'/','./');, eX_Str=strrep(eX_Str,'\','.\');
eX_Str=strrep(eX_Str,'^','.^');, eX_Str=strrep(eX_Str,'..','.');
n_Q=fl1-1; n_P=n_I-fl1;
varargin(fl1)=[]; n_I=n_I-1;
if n_Q<1 |~n_P,
   error('Not enough input arguments.');
end;

if strcmp(q_tYPe,'r')& n_Q>1,
   [varargin{1:n_Q}]=ndgrid(varargin{1:n_Q});
end;

for i_I_i=1:n_Q,
   N_iY_n=inputname(i_I_i); N_iY_v=varargin{i_I_i}; N_iY_v=N_iY_v(:);
   if isempty(N_iY_n), error('Parameter vector must be a named variable.'); end;
   eval([N_iY_n,'=N_iY_v;'],'error(peel(lasterr));');
   if i_I_i==1,
      N_iY_1=length(N_iY_v);
   elseif length(N_iY_v)~=N_iY_1,
      error('The lengths of exactly defined evaluate vectors must agree.');
   end;
end;
for i_I_i=n_Q+1:n_I,
   N_iY_n=inputname(i_I_i+1);
   if isempty(N_iY_n), error('Polynomial input argument must be a named variable.'); end;
   eval(['P_iY=polyval(varargin{i_I_i},oMegA_val);',...
         N_iY_n,'=P_iY(:)'';'],'error(peel(lasterr));');
   N_iY_2=length(P_iY(:));
   for j_J_j=1:6,
      c_iX_n=C_iY_n{j_J_j};
      eX_Str=strrep(eX_Str,[c_iX_n,N_iY_n],[c_iX_n(2:end),N_iY_n]);
   end;
end;

eX_s1=cumsum((eX_Str=='(')-(eX_Str==')'));
%eX_s2 -> is 1 for + or - not inside parentheses.
eX_s2=((eX_Str=='-') | (eX_Str=='+')) & (eX_s1==0);
f_E_2=find(eX_s2(2:end));

f_E_2=[0,f_E_2];  V_auX=[];  V_auX1=[];

for kk=fliplr(f_E_2),
   str_XA=['V_auX1=',eX_Str(kk+1:end),';'];
   eval(str_XA,'error(peel(lasterr));');
   eX_Str(kk+1:end)=[];
   if size(V_auX1,1)~=N_iY_1, V_auX1=repmat(V_auX1,[N_iY_1,1]); end;
   if size(V_auX1,2)==1, V_auX1=repmat(V_auX1,[1,N_iY_2]); end;
   if ~isempty(V_auX), V_auX=V_auX+V_auX1; 
   else,  V_auX=V_auX1; 
   end;
end;
V=reshape(V_auX,[prod(size(V_auX))/length(oMegA_val),length(oMegA_val)]);

%end .. vset
