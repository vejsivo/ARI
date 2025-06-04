function s = nth(n)
%NTH  Converts the number into a string
%     like  1st, 2nd, 3rd, 4th ...

%    Author:  J.Jezek  1999
%    Copyright(c) 1999 by Polyx, Ltd.

if n==1, t = 'st';
elseif n==2, t = 'nd';
elseif n==3, t = 'rd';
else t = 'th';
end;
s = [int2str(n),t];

%end .. @tsp/private/nth
