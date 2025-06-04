function S = cl2str(P)
%CL2STR  Returns the class and size of the object in a string.
%
% The commmand
%   S = CL2STR(P)
% returns the class and size of P as a short sentence in a string S.
% If  P is a short string (length<25) then S = P.
% If P is a vector of "double" and length(P)<25 then S = mat2str(P)

%       Author(s):  S. Pejchova  16-03-2000
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date:16-Mar-2000 13:50:50 $
%                         $Date:25-Jul-2002 J.Jezek arg check $

if nargin<1,
   error('Not enough input arguments.');
end;

S = P;
if ~ischar(S),
   if isa(S,'double') & ndims(S)==2 & (isempty(S) | (size(S,1)<=1 & size(S,2)<25)),
      if isempty(S) & ~isequal(size(S),[0 0]),
         S = sprintf('[%dx%d] double',size(S,1),size(S,2));
      else
         S = mat2str(S,3);
      end
   elseif isa(S,'cell') & isempty(S),
      if isequal(size(S),[0 0]),
         S = '{}';
      else
         S = sprintf('{%dx%d} cell',size(S,1),size(S,2));
      end
   else
      % Too big
      us_str = mat2str(size(S));
      us_str = strrep(us_str(2:end-1),' ','x');
      if isa(S,'cell'),
         us_str = ['{',us_str,'} ',class(S)];
      else
         us_str = ['[',us_str,'] ',class(S)];
      end
      S=us_str;
   end
elseif (length(S) > 25) ,
   S=sprintf('string of length %d',length(S));
end;

%end .. cl2str
