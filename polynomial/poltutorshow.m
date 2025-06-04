function slide=poltutorshow
% Getting started with Polynomial Toolbox 2.5
% ============================================
%
% This demo helps you develop some basic skills with Polynomial 
% Toolbox 2.5. The following topics are covered:
%
% 1.  Entering polynomials and polynomial matrices
% 2.  Concatenation and submatrices
% 3.  Extracting coefficients
% 4.  Algebraic operations: addition, subtraction, multiplication
% 5.  Determinants
% 6.  Unimodular polynomial matrices
% 7.  Adjoint and inverses
% 8.  Rank of a polynomial matrix.
% 9.  Null space and range of a polynomial matrix
% 10. Finite and infinite roots
% 11. Stability of polynomials and polynomial matrices
% 12. Division of polynomials and polynomial matrices
%
% See the Tutorial and Manual for the many operations not covered 
% in this demo.

% This is a slideshow file for use with playshow.m and makeshow.m
% To see it run, type 'playshow poltutorshow', 

% Author(s):    Z. Hurak, H. Kwakernaak   01-04-2001
% Copyright (c) 2001 by PolyX, Ltd.
% $Revision: 1.0.2 $    $Date: 01-10-2001 $
%-----------------------------------------------------------------------------------------

if nargout<1,
  playshow poltutorshow
else
  %========== Slide 1 ==========

  slide(1).code={
   'delete(findobj(gcf,''Position'',[ 0.04 0.55 0.65 0.4 ]));',
   'p0 = 0.5 + 2*s + 3*s^2 + 6*s^3 + 4*s^4 + 4*s^5 + s^6;',
   'spherplot(p0,0:.005:1.05,0.1);',
   'axis([-0.8 1.0 -0.5 1]);',
   'set(gca,''Visible'',''off'');',
   'hline = findobj(gca,''Type'',''line'',''Color'',[0 0 0],''LineStyle'','':'');',
   'delete(hline);',
   'htext = text(-0.95,0.9,''Getting started'');',
   'htext1 = text(-0.45,0.7,''with Polynomial Toolbox 2.5'');',
   'set(htext,''FontSize'',18);',
   'h_pa = patch([-0.64 -0.64 0.36 0.36],[-0.05 0.25 0.25 -0.05],[0.8 0.8 0.8]);',
   'set(h_pa,''EdgeColor'',''none'');',
   'htext2 = text(-0.61,0.19,''P = 1 + 2s^2   s^3'');',
   'htext3 = text(-0.41,0.0,''0.5 + s   3 + s'');',
   'set(htext2,''FontName'',''Courier New'',''FontSize'',12,''FontWeight'',''bold'',''Interpreter'',''none'',''Color'',[ 0 0.501960784313725 0 ]);',
   'set(htext3,''FontName'',''Courier New'',''FontSize'',12,''FontWeight'',''bold'',''Interpreter'',''none'',''Color'',[ 0 0.501960784313725 0 ]);',
   '',
   'P = [ 1+s 2*s^2; 2+s^3 4 ];',
   'h_CW = [];' };
  slide(1).text={
   'Polynomials and polynomial matrices are most easily entered using',
   'one of the indeterminate variables s, p, z, q, z^-1 or d that are',
   'recognized by the Polynomial Toolbox, combined with the usual MATLAB',
   'conventions for entering matrices.',
   '',
   'Thus you may easily define a simple polynomial matrix by typing',
   '',
   '>> P = [ 1+s 2*s^2; 2+s^3 4 ];',
   ''};

  %========== Slide 2 ==========

  slide(2).code={
   'delete(findobj(gcf,''Position'',[ 0.04 0.55 0.65 0.4 ]));' };
  slide(2).text={
   'In a system and control theoretic context, the indeterminates',
   's and p are usually associated with continuous time, and are',
   'closely related to the differential operator x(t) -> dx(t)/dt.',
   'The indeterminates z and q are associated with the discrete-time',
   'one-step shift operator x(t) ->  x(t+1), and the indeterminates',
   'd and z^-1 with the delay operator x(t) ->  x(t-1).',
   '',
   '>> Pa = 1 + d^2 - 3*d^3;',
   '>> Pb = 1 - q^4;'};

  %========== Slide 3 ==========

  slide(3).code={
   'P1 = 1 + z^-1 + z^-2;',
   'P2 = 1 + zi + zi^2;',
   'cla;',
   'set(gcf,''Units'',''normalized'');',
   'p1str = char(P1); p2str = char(P2);',
   'strDisp = sprintf(''%s\n%s\n\n%s\n%s'',''P1 ='',p1str,''P2 ='',p2str);',
   'h_CW = uicontrol(gcf,''Style'',''text'',''Units'',''normalized'',''Position'',[0.04 0.55 0.65 0.4],''String'',strDisp,''HorizontalAlignment'',''left'',''ForegroundColor'',[0 0 1],''Tag'',''CW'');',
   '',
   '' };
  slide(3).text={
   'The indeterminate z^-1 may be entered with arbitary nonnegative',
   'integral powers. For typing convenience z^-1 may be abbreviated to zi.',
   'Thus, the following  two commands return the same result:',
   '',
   '>> P1 = 1 + z^-1 + z^-2',
   '',
   '>> P2 = 1 + zi + zi^2',
   ''};

  %========== Slide 4 ==========

  slide(4).code={
   'P0 = [1 2; 3 4];',
   'P1 = [3 4; 5 1];',
   'P2 = [1 0; 0 1];',
   'P = pol([P0 P1 P2],2,''s'');',
   'Ps = char(P.'');',
   '[p1s,p2s,p3s,p4s] = deal(Ps{:});',
   'strDisp = sprintf(''%s\n%s\t%s\n%s\t\t%s'',''P ='',p1s,p2s,p3s,p4s);',
   'set(h_CW,''String'',strDisp);' };
  slide(4).text={
   'Polynomials and polynomial matrices can also be entered in terms',
   'of their coefficients or coefficient matrices. For this purpose the POL',
   'command is available.',
   '',
   '>> P0 = [1 2; 3 4];',
   '>> P1 = [3 4; 5 1];',
   '>> P2 = [1 0; 0 1];',
   '>> P = pol([P0 P1 P2],2,''s'')',
   '',
   ''};

  %========== Slide 5 ==========

  slide(5).code={
   'P = pol([P0 P1 P2],2);',
   'Ps = char(P.'');',
   '[p1s,p2s,p3s,p4s] = deal(Ps{:});',
   'strDisp = sprintf(''%s\n%s\t%s\n%s\t\t%s'',''P ='',p1s,p2s,p3s,p4s);',
   'set(h_CW,''String'',strDisp);',
   '' };
  slide(5).text={
   'After the Polynomial Toolbox has been started up, the default',
   'indeterminate variable is s. This implies, among other things, that',
   'the last command may also be entered as',
   '',
   '>> P = pol([P0 P1 P2],2)',
   '',
   ''};

  %========== Slide 6 ==========

  slide(6).code={
   'Q = [P; 1+s 3];',
   'p5s = char(1+s);',
   'p6s = char(pol(3));',
   'strDisp = sprintf(''%s\n%s\t%s\n%s\t\t%s\n%s\t\t%s'',''Q ='',p1s,p2s,p3s,p4s,p5s,p6s);',
   'set(h_CW,''String'',strDisp);',
   '' };
  slide(6).text={
   'Standard Matlab conventions may be used to concatenate polynomial',
   'and standard matrices.',
   '',
   '>> Q = [P; 1+s 3]',
   '',
   ''};

  %========== Slide 7 ==========

  slide(7).code={
   'Ps = char(Q(2:3,:).'');',
   '[p1s,p2s,p3s,p4s] = deal(Ps{:});',
   'strDisp = sprintf(''%s\n%s\t%s\n%s\t%s'',''ans ='',p1s,p2s,p3s,p4s);',
   'set(h_CW,''String'',strDisp);' };
  slide(7).text={
   'Submatrices may be selected in the usual way:',
   '',
   '>> Q(2:3,:)',
   ''};

  %========== Slide 8 ==========

  slide(8).code={
   'T = [ 1+s 2 s^2 s; 3 4 s^3 0 ];',
   '',
   'Ps = char(T.'');',
   '[p1s,p2s,p3s,p4s,p5s,p6s,p7s,p8s] = deal(Ps{:});',
   'T2s = T{2};',
   'T2s1 = num2str(T2s(1,:));',
   'T2s2 = num2str(T2s(2,:));',
   'strDisp = sprintf(''%s\n%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\n\n%s\n %s\n %s'',''T ='',p1s,p2s,p3s,p4s,p5s,p6s,p7s,p8s,''ans ='',T2s1,T2s2);',
   'set(h_CW,''String'',strDisp);' };
  slide(8).text={
   'It is easy to extract the coefficient matrices of a polynomial',
   'matrix. Consider the polynomial matrix',
   '',
   '>> T = [ 1+s 2 s^2 s; 3 4 s^3 0 ]',
   '',
   'The coefficient matrix of s^2 may be retrieved as',
   '',
   '>> T{2}'};

  %========== Slide 9 ==========

  slide(9).code={
   'T13s = num2str(T{:}(1,3));',
   'strDisp = sprintf(''%s\n %s'',''ans ='',T13s);',
   'set(h_CW,''String'',strDisp);',
   '' };
  slide(9).text={
   'The coefficients of the (1,3) entry of T follow as',
   '',
   '',
   '>> T{:}(1,3)',
   ''};

  %========== Slide 10 ==========

  slide(10).code={
   'P = [ 1+s 2; 3 4]; Q = [ s^2 s; s^3 0];',
   'S = P+Q;',
   'R = P*Q;',
   'Ps = char(S.'');',
   'Ks = char(R.'');',
   '[p1s,p2s,p3s,p4s] = deal(Ps{:});',
   '[p5s,p6s,p7s,p8s] = deal(Ks{:});',
   'strDisp = sprintf(''%s\n%s\t%s\n%s\t\t%s\n\n%s\n%s\t%s\n%s\t%s'',''S ='',p1s,p2s,p3s,p4s,''R ='',p5s,p6s,p7s,p8s);',
   'set(h_CW,''String'',strDisp);',
   '',
   '' };
  slide(10).text={
   'Algebraic operations with polynomial matrices also follow the Matlab syntax. Consider the polynomial matrices',
   '',
   '>> P = [ 1+s 2; 3 4]; Q = [ s^2 s; s^3 0];',
   '',
   'The sum and product of P and Q follow',
   '',
   '>> S = P + Q',
   '>> R = P*Q',
   ''};

  %========== Slide 11 ==========

  slide(11).code={
   'P = [1 s s^2; 1+s s 1-s; 0 -1 -s];',
   'detP = det(P);',
   'Ps = char(P.'');',
   'ds = char(detP);',
   '[p1s,p2s,p3s,p4s,p5s,p6s,p7s,p8s,p9s] = deal(Ps{:});',
   'strDisp = sprintf(''%s\n%s\t%s\t%s\n%s\t%s\t%s\n%s\t%s\t%s\n\n\n%s\n%s'',''P ='',p1s,p2s,p3s,p4s,p5s,p6s,p7s,p8s,p9s,''detP ='',ds);',
   'set(h_CW,''String'',strDisp);',
   '',
   '' };
  slide(11).text={
   'The determinant of a square polynomial matrix is defined exactly',
   'as its constant matrix counterpart. In fact, its computation is not',
   'much more difficult:',
   '',
   '>> P = [1 s s^2; 1+s s 1-s; 0 -1 -s]',
   '',
   '>> detP = det(P)'};

  %========== Slide 12 ==========

  slide(12).code={
   'U = [ 2-s-2*s^2 2-2*s^2 1+s; 1-s-s^2 1-s^2 s; -1-s -s 1];',
   'detU = det(U);',
   'Ps = char(U.'');',
   'ds = char(detU);',
   '[p1s,p2s,p3s,p4s,p5s,p6s,p7s,p8s,p9s] = deal(Ps{:});',
   'strDisp = sprintf(''%s\n%s\t%s\t%s\n%s\t%s\t\t%s\n%s\t\t%s\t\t%s\n\n\n%s\n%s'',''U ='',p1s,p2s,p3s,p4s,p5s,p6s,p7s,p8s,p9s,''detU ='',ds);',
   'set(h_CW,''String'',strDisp);',
   '' };
  slide(12).text={
   'If its determinant happens to be constant then the polynomial',
   'matrix is called unimodular:',
   '',
   '>> U = [ 2-s-2*s^2 2-2*s^2 1+s; 1-s-s^2 1-s^2 s; -1-s -s 1]',
   '>> detU = det(U)',
   ''};

  %========== Slide 13 ==========

  slide(13).code={
   'strDisp = sprintf(''%s\n %s'', ''ans ='',num2str(isunimod(U)));',
   'set(h_CW,''String'',strDisp);' };
  slide(13).text={
   'If a matrix is suspected of unimodularity then one may make',
   'sure by a special tester ISUNIMOD:',
   '',
   '>> isunimod(U)'};

  %========== Slide 14 ==========

  slide(14).code={
   'adjP = adj(P);',
   'Ps = char(adjP.'');',
   '[p1s,p2s,p3s,p4s,p5s,p6s,p7s,p8s,p9s] = deal(Ps{:});',
   'strDisp = sprintf(''%s\n%s\t%s\t%s\n%s\t\t%s\t%s\n%s\t\t%s\t%s'',''adjP ='',p1s,p2s,p3s,p4s,p5s,p6s,p7s,p8s,p9s);',
   'set(h_CW,''String'',strDisp);',
   '' };
  slide(14).text={
   'Also the adjoint matrix is defined as for constant matrices.',
   'The adjoint is a polynomial matrix and may be computed by typing',
   '',
   '>> adjP = adj(P)',
   ''};

  %========== Slide 15 ==========

  slide(15).code={
   '[num,den] = inv(P);',
   'I = P*num/den ;',
   'Is = char(I.'');',
   'Ps = char(num.'');',
   'ds = char(den);',
   '[p1s,p2s,p3s,p4s,p5s,p6s,p7s,p8s,p9s] = deal(Ps{:});',
   '[i1s,i2s,i3s,i4s,i5s,i6s,i7s,i8s,i9s] = deal(Is{:});',
   'strDisp = sprintf(''%s\n%s\t%s\t%s\n%s\t\t%s\t%s\n%s\t\t%s\t%s\n\n%s\n%s\n\n%s\n%s\t%s\t%s\n%s\t%s\t%s\n%s\t%s\t%s'',''num ='',p1s,p2s,p3s,p4s,p5s,p6s,p7s,p8s,p9s,''den ='',ds,''ans ='',i1s,i2s,i3s,i4s,i5s,i6s,i7s,i8s,i9s);',
   'set(h_CW,''String'',strDisp);',
   '',
   '' };
  slide(15).text={
   'The inverse of a polynomial matrix P is given by num/den,',
   'where num is a polynomial matrix and den is a scalar monic',
   'polynoimal. Hence, it follows by typing',
   '',
   '>> [num,den] = inv(P)',
   '',
   '>> P*num/den',
   '',
   'Note that this example involves polynomial division.',
   '',
   '',
   ''};

  %========== Slide 16 ==========

  slide(16).code={
   'Ps = char(P.'');',
   '[p1s,p2s,p3s,p4s,p5s,p6s,p7s,p8s,p9s] = deal(Ps{:});',
   'strDisp = sprintf(''%s\n%s\t%s\t%s\n%s\t%s\t%s\n%s\t%s\t%s\n\n%s\n %s'',''P ='',p1s,p2s,p3s,p4s,p5s,p6s,p7s,p8s,p9s,''ans='',num2str(isfullrank(P)));',
   'set(h_CW,''String'',strDisp);' };
  slide(16).text={
   'A polynomial matrix P(s) has full column rank (or full normal',
   'column rank) if it has full column rank everywhere in the complex',
   'plane except at a finite number of points. Similar definitions hold',
   'for full row rank and full rank. Recall that',
   '',
   '>> P',
   '',
   'The following rank test confirms that P has full rank:',
   '',
   '>> isfulrank(P)'};

  %========== Slide 17 ==========

  slide(17).code={
   'strDisp = sprintf(''%s\n %s'',''ans ='',num2str(rank(P)));',
   'set(h_CW,''String'',strDisp);' };
  slide(17).text={
   'The normal rank of a polynomial matrix P(s) equals max rank P(s)',
   'where s ranges over the complex numbers. The rank is calculated by',
   'typing',
   '',
   '>> rank(P)'};

  %========== Slide 18 ==========

  slide(18).code={
   'strDisp = sprintf(''%s\n %s'',''ans ='',num2str(issingular(P)));',
   'set(h_CW,''String'',strDisp);' };
  slide(18).text={
   'A polynomial matrix is nonsingular if it has full normal rank:',
   '',
   '>> issingular(P)'};

  %========== Slide 19 ==========

  slide(19).code={
   'A = P(1:2,:);',
   '',
   'N = null(A);',
   'Ps = char(A.'');',
   'Ns = char(N.'');',
   '[p1s,p2s,p3s,p4s,p5s,p6s] = deal(Ps{:});',
   '[p7s,p8s,p9s] = deal(Ns{:});',
   'strDisp = sprintf(''%s\n%s\t%s\t%s\n%s\t%s\t%s\n\n%s\n%s\n%s\n%s'',''A ='',p1s,p2s,p3s,p4s,p5s,p6s,''N ='',p7s,p8s,p9s);',
   'set(h_CW,''String'',strDisp);',
   '',
   '',
   '' };
  slide(19).text={
   'There are two important subspaces (more precisely, submodules)',
   'associated with a polynomial matrix A(s): its null space and its',
   'range (or span). The (right) null space is defined as the set of',
   'all polynomial vectors x(s) such that A(s)x(s) = 0.',
   '',
   '>> A = P(1:2,:)',
   '',
   '>> N = null(A)',
   ''};

  %========== Slide 20 ==========

  slide(20).code={
   'AN = A*N;',
   'ANs = char(AN.'');',
   '[a1s,a2s] = deal(ANs{:});',
   'strDisp = sprintf(''%s\n %s\n %s'',''ans ='',a1s,a2s);',
   'set(h_CW,''String'',strDisp);',
   '',
   '' };
  slide(20).text={
   'Hence the null space dimension is 1 and its basis has degree 3.',
   'Check:',
   '',
   '>> A*N'};

  %========== Slide 21 ==========

  slide(21).code={
   'MB = minbasis(A);',
   'Ps = char(MB.'');',
   '[p1s,p2s,p3s,p4s] = deal(Ps{:});',
   'strDisp = sprintf(''%s\n%s\t%s\n%s\t%s'',''P ='',p1s,p2s,p3s,p4s);',
   'set(h_CW,''String'',strDisp);' };
  slide(21).text={
   'The range of A(s) is the set of all polynomial vectors y(s)',
   'such that y(s) = A(s)x(s) for some polynomial vector x(s).',
   'In the Polynomial Toolbox, the minimal basis of the range is',
   'returned by the command',
   '',
   '>> minbasis(A)',
   '',
   'The columns of this matrix form a minimal basis.'};

  %========== Slide 22 ==========

  slide(22).code={
   'r = roots(P);',
   'rs  = num2str(r);',
   'strDisp = sprintf(''%s\n %s\n  %s'',''ans ='',rs(1,:),rs(2,:));',
   'set(h_CW,''String'',strDisp);' };
  slide(22).text={
   'The roots or zeros of a polynomial matrix P(s) are those points',
   's in the complex plane where P(s) loses rank:',
   '',
   '>> roots(P)'};

  %========== Slide 23 ==========

  slide(23).code={
   'r = roots(P,''all'');',
   'rs  = num2str(r);',
   'strDisp = sprintf(''%s\n %s\n %s\n %s'',''ans ='',rs(1,:),rs(2,:),rs(3,:));',
   'set(h_CW,''String'',strDisp);' };
  slide(23).text={
   'The roots can be both finite and infinite. The infinite',
   'roots are normally suppressed. To reveal them, type',
   '',
   '>> roots(P,''all'')'};

  %========== Slide 24 ==========

  slide(24).code={
   'r = roots(U,''all'');',
   'rs  = num2str(r);',
   'strDisp = sprintf(''%s\n %s\n %s\n %s'',''ans ='',rs(1,:),rs(2,:),rs(3,:));',
   'set(h_CW,''String'',strDisp);' };
  slide(24).text={
   'Unimodular matrices have no finite roots:',
   '',
   '>> roots(U,''all'')'};

  %========== Slide 25 ==========

  slide(25).code={
   'isS = isstable(s-0.5);',
   'isZ = isstable(z-0.5);',
   '',
   'isSs = num2str(isS);',
   'isSz = num2str(isZ);',
   '',
   'strDisp = sprintf(''%s\n %s\n\n\n%s\n %s'',''ans ='',isSs,''ans ='',isSz);',
   'set(h_CW,''String'',strDisp);',
   '',
   '' };
  slide(25).text={
   'The macro isstable checks stability according to the',
   'variable symbol. Thus,',
   '',
   '>> isstable(s-0.5)',
   '',
   '>> isstable(z-0.5)'};

  %========== Slide 26 ==========

  slide(26).code={
   'b = 1-s; c = 1+s; a = b*c;',
   'set(h_CW,''String'',[]);',
   '',
   '' };
  slide(26).text={
   'To understand when division of polynomials and polynomial',
   'matrices is possible, consider three polynomials a(s), b(s)',
   'and c(s) such that a(s) = b(s)c(s). We say that b(s) is a divisor',
   '(or factor) of a(s) or a(s) is a multiple of b(s), and write',
   'a(s)|b(s). This is sometimes also stated as b(s) divides a(s).',
   '',
   '>> b = 1-s; c = 1+s; a = b*c;',
   ''};

  %========== Slide 27 ==========

  slide(27).code={
   'ab = char(a/b);',
   'ac = char(a/c);',
   '',
   'strDisp = sprintf(''%s\n %s\n\n\n%s\n %s'',''ans ='',ab,''ans ='',ac);',
   'set(h_CW,''String'',strDisp);' };
  slide(27).text={
   'As b(s) and c(s) are both divisors of a(s), the following divisions',
   'both can be done:',
   '',
   '>> a/b',
   '',
   '>> a/c',
   ''};

  %========== Slide 28 ==========

  slide(28).code={
   'cb = char(c/b);',
   '',
   'strDisp = sprintf(''%s\n %s %s'',''ans ='',cb);',
   'set(h_CW,''String'',strDisp);' };
  slide(28).text={
   'Of course, the division by b(s) fails if b(s) is not a divisor:',
   '',
   '>> c/b'};

  %========== Slide 29 ==========

  slide(29).code={
   '[q,r] = rdiv(c,b);',
   '',
   'qs = char(q);',
   'rs = char(r);',
   '',
   'strDisp = sprintf(''%s\n %s\n\n\n%s\n %s'',''q ='',qs,''r ='',rs);',
   'set(h_CW,''String'',strDisp);',
   '' };
  slide(29).text={
   'The quotient and remainder of this division may be retrieved by typing',
   '',
   '>> [q,r] = rdiv(c,b)',
   '',
   'The Polynomial Toolbox performs these divisions not only',
   'for polynomials but also for polynomial matrices.',
   '',
   'See the Tutorial and Manual for the many operations not',
   'covered in this demo.',
   ''};
end

%end .. poltutorshow
