function s2 = peelf(s1)
%PEELF    Peel off the first line in string
%            s2 = peel(s1)

%     Author:  J.Jezek  1999
%     Copyright(c) 1999 by Polyx, Ltd.

k = findstr(s1,char(10));
s2 = s1(k+1:end);

%end .. @tsp/private/peelf
