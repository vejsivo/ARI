function slide=robparshow
% PARAMETRIC UNCERTAINTIES - DEMO FOR POLYNOMIAL TOOLBOX 2.5
% ==========================================================
% This demo introduces some analytical tools for systems 
% with uncertain parameters. The following topics are covered:
%
% 1. Robust stability interval for one uncertain polynomial.
% 2. Continuous-time interval polynomials by Kharitonov theorem.
% 3. Continuous-time interval polynomials by graphical methods.
% 4. Discrete-time interval polynomials by graphical methods.
% 5. General polytope of polynomials.
% 6. Uncertain polynomials with multilinear uncertainty structure.
% 7. Polynomial uncertainty structure.
% 8. General uncertainty structure.
% 9. Spherical polynomial families.
%
% This is a slideshow file for use with playshow.m and makeshow.m
% To see it run, type 'playshow robparshow', 

% % Author(s):    Z. Hurak   01-10-2001
% Copyright (c) 2001 by PolyX, Ltd.
% $Revision: 1.0.0 $  $Date: 01-10-2001 $

if nargout<1,
  playshow robparshow
else
  %========== Slide 1 ==========

  slide(1).code={
   'set(findobj(''Tag'',''CW''),''Visible'',''off'');',
   'p0 = 0.5 + 2*s + 3*s^2 + 6*s^3 + 4*s^4 + 4*s^5 + s^6;',
   'spherplot(p0,0:.005:1.05,0.1);',
   'axis([-0.8 1.0 -0.5 1]);',
   'set(gca,''Visible'',''off'');',
   'hline = findobj(gca,''Type'',''line'',''Color'',[0 0 0],''LineStyle'','':'');',
   'delete(hline);',
   'htext = text(-0.95,0.9,''Parametric Uncertainties'');',
   'htext1 = text(-0.48,0.7,''demo for Polynomial Toolbox 2.5'');',
   'set(htext,''FontSize'',18);',
   'h_pa = patch([-0.64 -0.64 0.36 0.36],[-0.05 0.25 0.25 -0.05],[0.8 0.8 0.8]);',
   'set(h_pa,''EdgeColor'',''none'');',
   'htext2 = text(-0.61,0.19,''P = 1 + 2s^2   s^3'');',
   'htext3 = text(-0.41,0.0,''0.5 + s   3 + s'');',
   'set(htext2,''FontName'',''Courier New'',''FontSize'',12,''FontWeight'',''bold'',''Interpreter'',''none'',''Color'',[ 0 0.501960784313725 0 ]);',
   'set(htext3,''FontName'',''Courier New'',''FontSize'',12,''FontWeight'',''bold'',''Interpreter'',''none'',''Color'',[ 0 0.501960784313725 0 ]);' };
  slide(1).text={
   '1. ROBUST STABILITY INTERVAL FOR ONE UNCERTAIN PARAMETER',
   '============================================================',
   'Consider a characteristic polynomial whose coefficients depend affinely',
   'on one uncertain parameter. Such a polynomial can be written:',
   '',
   'p(s,q)=p0(s)+q*p1(s)',
   '',
   'where p0(s) is a stable polynomial,  q is the uncertain parameter.',
   'Now, find the maximum stability interval [qmin, qmax].',
   ''};

  %========== Slide 2 ==========

  slide(2).code={
   'cla;',  
   'p0 = 3 + 10*s + 12*s^2 + 6*s^3 + s^4;',
   'p1 = s + s^3;'
   'strDisp = sprintf(''%s\n      %s'',''isstable ='',num2str(isstable(p0)));',
   'set(gcf,''Units'',''normalized'');',
   'h_CW =uicontrol(gcf,''Style'',''Text'',''Units'',''normalized'',''Position'',[0.04 0.55 0.65 0.4],''String'',strDisp,''HorizontalAlignment'',''left'',''ForegroundColor'',[0 0 1],''Tag'',''CW'',''FontName'',''Courier'',''FontSize'',14);' };
  slide(2).text={
   'Enter the two polynomials p0(s) and p1(s)',
   '',
   '>> p0 = 3 + 10*s + 12*s^2 + 6*s^3 + s^4;',
   '>> p1 = s + s^3;',
   '',
   'First, stability of p0(s) is tested:',
   '',
   '>> isstable = isstable(p0)',
   ''};

  %========== Slide 3 ==========

  slide(3).code={
   'set(cla,''Visible'',''off'');',
   'cla;',
   '[qmin,qmax]=stabint(p0,p1);',
   'strDisp = sprintf(''%s\n   %s\n\n%s\n   %s'',''qmin ='',qmin,''qmax ='',qmax);',
   'set(findobj(''Tag'',''CW''),''String'',strDisp,''Visible'',''on'');',
   '',
   '' };
  slide(3).text={
   'Finding the stability interval is accomplished with the command STABINT:',
   '',
   '>> [qmin,qmax]=stabint(p0,p1)',
   ''};

  %========== Slide 4 ==========

  slide(4).code={
   'set(findobj(''Tag'',''CW''),''Visible'',''off'');',
   'rlocus(ss(p1,p0),qmin:.1:100);',
   '' };
  slide(4).text={
   'We can check this result by plotting the root locus of a (artificial) plant p1/p0',
   'and feedback gain q.',
   'This requires that the Control Toolbox 4.2 or higher be installed on your computer:',
   '',
   '>> rlocus(ss(p1,p0),qmin:.1:100)',
   '',
   ''};

  %========== Slide 5 ==========

  slide(5).code={
   'set(cla,''Visible'',''off'');',
   'cla;',
   '' };
  slide(5).text={
   '2. CONTINUOUS-TIME INTERVAL POLYNOMIALS BY KHARITONOV THEOREM',
   '============================================================',
   'Consider the continous-time characteristic polynomial, whose coefficients',
   'are only know to lie within some intervals:',
   '',
   ' p(s,q)=[0.45,0.55] + [1.95,2.05]*s + [2.95,3.05]*s^3 + [3.95,4.05]*s^4 + [3.95,4.05]*s^5 + s^6',
   '',
   'Check if it is robustly stable.',
   ''};

  %========== Slide 6 ==========

  slide(6).code={
   'pminus = 0.45+1.95*s+2.95*s^2+5.95*s^3+3.95*s^4+3.95*s^5+s^6;',
   'pplus = 0.55+2.05*s+3.05*s^2+6.05*s^3+4.05*s^4+4.05*s^5+s^6;',
   'strDisp = sprintf(''%s\n%s\n%s\n%s'',''pminus ='',char(pminus),''pplus ='',char(pplus));',
   'set(findobj(''Tag'',''CW''),''String'',strDisp,''FontSize'',12);' };
  slide(6).text={
   'Enter the interval polynomial via two ''lumped'' polynomials:',
   '',
   '>> pminus = 0.45+1.95*s+2.95*s^2+5.95*s^3+3.95*s^4+3.95*s^5+s^6;',
   '>> pplus = 0.55+2.05*s+3.05*s^2+6.05*s^3+4.05*s^4+4.05*s^5+s^6;',
   ''};

  %========== Slide 7 ==========

  slide(7).code={
   '[stability,K1,K2,K3,K4]=kharit(pminus,pplus);',
   'sta = num2str(stability);',
   'K1 = char(K1); K2 = char(K2); K3 = char(K3); K4 = char(K4);',
   'strDisp = sprintf(''%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s'',''stability ='',sta,''K1 ='',K1,''K2 ='',K2,''K3 ='',K3,''K4 ='',K4);',
   'set(findobj(''Tag'',''CW''),''String'',strDisp,''Visible'',''on'',''FontSize'',10);',
   '' };
  slide(7).text={
   'Theorem of Kharitonov says that it is sufficient to test for the stability',
   'of the four distinguished polynomials. To obtain these so-called Kharitonov',
   'polynomials we use the KHARIT function:',
   '',
   '>> [stability,K1,K2,K3,K4]=kharit(pminus,pplus)',
   ''};

  %========== Slide 8 ==========

  slide(8).code={
   'cla;',
   'set(gca,''Visible'',''off'');',
   'ispm = isstable(pminus);',
   'strDisp = sprintf(''%s\n   %s'',''ans ='',num2str(ispm));',
   'set(findobj(''Tag'',''CW''),''FontSize'',12,''String'',strDisp,''Visible'',''on'');' };
  slide(8).text={
   '3. CONTINUOUS-TIME INTERVAL POLYNOMIALS BY GRAPHICAL METHODS',
   '============================================================',
   'For the same data as in the previous example, use the graphical',
   'method of Zero Exclusion Principle to test the robust stability.',
   '',
   'First, check if there is at least one stable member in the polynomial family:',
   '',
   '>> isstable(pminus)',
   ''};

  %========== Slide 9 ==========

  slide(9).code={
   'set(findobj(''Tag'',''CW''),''Visible'',''off'');',
   'set(cla,''Visible'',''on'');',
   'khplot(pminus,pplus,0:.002:1);' };
  slide(9).text={
   'Plot the value sets (here Kharitonov rectangles) on imaginary axis.',
   'Let''s make the initial guess at he frequency range OMEGA between 0 and 1:',
   '',
   '>> khplot(pminus,pplus,0:.002:1)',
   ''};

  %========== Slide 10 ==========

  slide(10).code={
   ' p0 = 10 + 20*z + 128*z^2 + 260*z^3 + 168*z^4;',
   ' p1=pol(1); p2=z; p3=z^2; p4=z^3;',
   'cla;',
   'set(gca,''Visible'',''off'');',
   'strDisp = sprintf(''%s\n %s\n%s\n %s\n%s\n %s\n%s\n %s\n%s\n %s'',''p0 ='',char(p0),''p1 ='',char(p1),''p2 ='',char(p2),''p3 ='',char(p3),''p4 ='',char(p4));',
   'set(findobj(''Tag'',''CW''),''FontSize'',10,''Visible'',''on'',''String'',strDisp);' };
  slide(10).text={
   '4.  DISCRETE-TIME INTERVAL POLYNOMIALS BY GRAPHICAL METHODS',
   '===============================================================',
   'Consider a discrete-time interval polynomial given by:',
   '',
   'p(z,q) = [10,20] + [20,30]*z + [128,138]*z^2 + [260,270]*z^3 + [168,178]*z^4',
   '',
   'It can be described by the extreme polynomials:',
   '>> p0 = 10 + 20*z + 128*z^2 + 260*z^3 + 168*z^4;',
   '>> p1=1; p2=z; p3=z^2; p4=z^3;',
   ''};

  %========== Slide 11 ==========

  slide(11).code={
   'cla;',
   'set(gca,''Visible'',''off'');',
   'Qbounds=[0 10;0 10;0 10;0 10];',
   'ii = 1:4;',
   'Qstr(ii,:) = num2str(Qbounds(ii,:));',
   'strDisp = sprintf(''%s\n  %s\n  %s\n  %s\n  %s'',''Qbounds ='',Qstr(1,:),Qstr(2,:),Qstr(3,:),Qstr(4,:));',
   '',
   'set(findobj(''Tag'',''CW''),''FontSize'',12,''Visible'',''on'',''String'',strDisp);',
   '' };
  slide(11).text={
   'and the bounds on the uncertain parameters are given by',
   '',
   '>> Qbounds=[0 10;0 10;0 10;0 10]',
   '',
   ''};

  %========== Slide 12 ==========

  slide(12).code={
   'set(findobj(''Tag'',''CW''),''Visible'',''off'');',
   'set(gca,''Visible'',''on'');',
   'ptopplot(p0,p1,p2,p3,p4,Qbounds,exp(j*(0.3:0.002:0.7)*2*pi));' };
  slide(12).text={
   'Now, we can easily plot the value sets (octagons) for generalized frequencies',
   'on unit circle using the function PTOPPLOT:',
   '',
   '>> ptopplot(p0,p1,p2,p3,p4,Qbounds,exp(j*(0.3:0.002:0.7)*2*pi))',
   '',
   'we see that zero is excluded from the value set, so we conclude that the uncertain polynomial is robustly stable.',
   '',
   ''};

  %========== Slide 13 ==========

  slide(13).code={
   'cla;',
   'set(gca,''Visible'',''off'');' };
  slide(13).text={
   '5. GENERAL POLYTOPE OF POLYNOMIALS',
   '========================================',
   'We can define even more general polynomial families by entering the set of its extreme members, so-called generators of the polytope. This is especially useful, when the uncertain parameters enter the coefficients of the uncertain polynomial affinely.',
   ''};

  %========== Slide 14 ==========

  slide(14).code={
   'set(findobj(''Tag'',''CW''),''Visible'',''off'');' };
  slide(14).text={
   'Is the following uncertain polynomial robustly stable?',
   '',
   'p(s,q) = 3+q1-q2+2*q3 + (3+3*q1+q2+q3)*s + (3+3*q1-3*q2+q3)*s^2 + (1+2*q1-q2+2*q3)*s^3',
   '',
   'with',
   '           -0.245 <= q1 <= 0.245',
   '           -0.245 <= q2 <= 0.245',
   '           -0.245 <= q3 <= 0.245'};

  %========== Slide 15 ==========

  slide(15).code={
   'cla;',
   'set(gca,''Visible'',''off'');',
   'p0=pol([3 3 3 1],3);',
   'p1=pol([1 3 3 2],3);',
   'p2=pol([-1 1 -3 -1],3);',
   'p3=pol([2 1 1 2],3);',
   'Qbounds=[-0.245 0.245;-0.245 0.245;-0.245 0.245];',
   '',
   '' };
  slide(15).text={
   'First, enter the data:',
   '',
   '>> p0=pol([3 3 3 1],3);',
   '>> p1=pol([1 3 3 2],3);',
   '>> p2=pol([-1 1 -3 -1],3);',
   '>> p3=pol([2 1 1 2],3);',
   '>> Qbounds=[-0.245 0.245;-0.245 0.245;-0.245 0.245];',
   ''};

  %========== Slide 16 ==========

  slide(16).code={
   'set(findobj(''Tag'',''CW''),''Visible'',''off'');',
   'set(gca,''Visible'',''on'');',
   '',
   'ptopplot(p0,p1,p2,p3,Qbounds,j*(0:1.5/50:1.5));',
   '' };
  slide(16).text={
   'Plot the value sets (octagons) on imaginary axis  with the initial guess at the range of frequencies from 0 to 1.5',
   '',
   '>> ptopplot(p0,p1,p2,p3,Qbounds,j*(0:1.5/50:1.5))',
   '',
   'we see that zero is excluded from the value set (it can be verified that it holds',
   'even if we increase the frequency range), so we conclude that the uncertain',
   'polynomial is robustly stable.',
   ''};

  %========== Slide 17 ==========

  slide(17).code={
   'cla;',
   'set(gca,''Visible'',''off'');',
   '' };
  slide(17).text={
   '6. UNCERTAIN POLYNOMIALS WITH MULTILINEAR UNCERTAINTY STRUCTURE',
   '==============================================================',
   'If the uncertain parameters enter the coefficients in multilinear way,',
   'the nice theory from the previous examples doesn''t apply. There are essentially',
   'two approaches to robust stability analysis of multilinear systems.',
   '',
   '1) Overbounding by the polytope of polynomials and using the tools from the previous examples. This however, introduces some conservatism into analysis.',
   ''};

  %========== Slide 18 ==========

  slide(18).code={
   '' };
  slide(18).text={
   '2) Gridding the set of uncertain parameters and plotting the value set polygons for a range of frequencies. Visually checking if the zero is excluded from the value set, we can make nonconservative conclusion about robust stability. Since this is a ''brute force'' approach it will be very time consuming for high order systems.',
   ''};

  %========== Slide 19 ==========

  slide(19).code={
   'p0=1+s+s^2; p1=s; p2=1;',
   'q1=0:.01:1;',
   'q2=0:.01:1;',
   '' };
  slide(19).text={
   'Consider the uncertain polynomial with multilinear uncertainty structure described by p(s,q) = p0(s)+(q1+q2)*p2(s)+q1*q2*p1(s),  where',
   '',
   '>> p0=1+s+s^2; p1=s; p2=1;',
   '',
   'and the uncertain parameters range from 0 to 1:',
   '',
   '>> q1=0:.01:1;',
   '>> q2=0:.01:1;',
   ''};

  %========== Slide 20 ==========

  slide(20).code={
   'set(gca,''Visible'',''on'');',
   'V=vset(q1,q2,''p0+(q1+q2)*p2+q1*q2*p1'',p0,p1,p2,j*(0:.1:2));',
   'vsetplot(V);' };
  slide(20).text={
   'Computation of the value sets for a range of frequencies (chosen from 0 to 2)',
   'is accomplished by the VSET and VSETPLOT functions:',
   '',
   '>> V=vset(q1,q2,''p0+(q1+q2)*p2+q1*q2*p1'',p0,p1,p2,j*(0:.1:2));',
   '',
   '>> vsetplot(V);',
   ''};

  %========== Slide 21 ==========

  slide(21).code={
   'cla;',
   'set(gca,''Visible'',''off'');' };
  slide(21).text={
   '7. POLYNOMIAL UNCERTAINTY STRUCTURE',
   '========================================',
   'The polynomial family exhibits a polynomial uncertainty structure',
   'when the uncertain parameters enter in higher powers into the coefficients',
   'of the polynomial. We can use the same tool as for multilinear case - the two',
   'functions VSET and VSETPLOT accomplish the computation and plotting of the value sets by simple gridding of the set of uncertain parameters.',
   ''};

  %========== Slide 22 ==========

  slide(22).code={
   '' };
  slide(22).text={
   'Consider the uncertain polynomial',
   '',
   'p(s,q) = p0(s)+(q1^3+q1)*p1(s)+(q2^3-q2^2)*p2(s)+q1*q2*p3(s)',
   '',
   'with the two uncertain parameters q1 and q2 ranging from -1 to 1. Is the polynomial family robustly stable?',
   ''};

  %========== Slide 23 ==========

  slide(23).code={
   'cla;',
   'set(gca,''Visible'',''off'');',
   'uncertStruct=''p0+(q1^3+q1)*p1+(q2^3-q2^2)*p2+q1*q2*p3'';',
   'p0=1+s+s^2;',
   'p1=1+s;',
   'p2=1-s;',
   'p3=1-s+s^2;',
   'q1=-1:.02:1; q2=-1:.02:1;',
   '' };
  slide(23).text={
   'Enter the uncetainty structure as a string:',
   '',
   '>> uncertStruct=''p0+(q1^3+q1)*p1+(q2^3-q2^2)*p2+q1*q2*p3'';',
   '>> p0=1+s+s^2;',
   '>> p1=1+s;',
   '>> p2=1-s;',
   '>> p3=1-s+s^2;',
   '>> q1=-1:.02:1; q2=-1:.02:1;',
   ''};

  %========== Slide 24 ==========

  slide(24).code={
   'set(gca,''Visible'',''on'');',
   'V=vset(q1,q2,uncertStruct,p0,p1,p2,p3,j*(0:.4:4));',
   'vsetplot(V);',
   '',
   '' };
  slide(24).text={
   'Computing and plotting the value sets on imaginary axis for frequencies from 0 to 4 is easy with the VSET and VSETPLOT functions',
   '',
   '>> V=vset(q1,q2,uncertStruct,p0,p1,p2,p3,j*(0:.4:4));',
   '',
   '>> vsetplot(V);',
   ''};

  %========== Slide 25 ==========

  slide(25).code={
   'cla;',
   'set(gca,''Visible'',''off'');',
   '' };
  slide(25).text={
   '8. GENERAL UNCERTAINTY STRUCTURE',
   '=====================================',
   'The two functions for computing and plotting the values sets for a range of frequencies can be used for robust stability analysis of polynomial families with general uncertainty structures. Actually, the ''gridding'' approach of the two functions will be the only solution in such cases.',
   '',
   'Even though it might be rather time consuming, it is the only tool available for analysis of such uncertain systems.',
   ''};

  %========== Slide 26 ==========

  slide(26).code={
   'cla;',
   'set(gca,''Visible'',''off'');',
   'uncertStruct=''p0+sqrt(abs(q1))*p1+exp(q2)*p2+sin(q1*q2)*p3'';',
   'p0=1+s+s^2;',
   'p1=1+s;',
   'p2=1-s;',
   'p3=1-s+s^2;',
   'q1=-1:.02:1; q2=-1:.02:1;',
   '' };
  slide(26).text={
   'Consider the uncertain polynomial described:',
   '',
   '>> uncertStruct=''p0+sqrt(abs(q1))*p1+exp(q2)*p2+sin(q1*q2)*p3'';',
   '>> p0=1+s+s^2;',
   '>> p1=1+s;',
   '>> p2=1-s;',
   '>> p3=1-s+s^2;',
   '>> q1=-1:.02:1; q2=-1:.02:1;',
   ''};

  %========== Slide 27 ==========

  slide(27).code={
   'V=vset(q1,q2,uncertStruct,p0,p1,p2,p3,j*(0:.4:4));',
   'vsetplot(V);',
   '' };
  slide(27).text={
   'Computing and plotting the value sets on imaginary axis for frequencies from 0 to 4 is done by the VSET and VSETPLOT functions',
   '',
   '>> V=vset(q1,q2,uncertStruct,p0,p1,p2,p3,j*(0:.4:4));',
   '',
   '>> vsetplot(V);',
   ''};

  %========== Slide 28 ==========

  slide(28).code={
   'uncertStruct=''p0+sqrt(abs(q1))*p1+sin(q2)*p2+cos(q1*q2)*p3'';',
   'V=vset(q1,q2,uncertStruct,p0,p1,p2,p3,j*(0:.4:4));vsetplot(V)',
   '' };
  slide(28).text={
   'the uncertainty structure described in the string UNCERTSTRUCT is a very flexible way of entering whatever kind of dependence of coefficients on the uncertain parameters.',
   '',
   'Let''s illustrate this by presenting a few more examples of uncertainty structures,',
   'using the same data.',
   '',
   '>> uncertStruct=''p0+sqrt(abs(q1))*p1+sin(q2)*p2+cos(q1*q2)*p3'';',
   '>> V=vset(q1,q2,uncertStruct,p0,p1,p2,p3,j*(0:.4:4));vsetplot(V)',
   ''};

  %========== Slide 29 ==========

  slide(29).code={
   'uncertStruct=''p0+q2*p1+sinh(q2)*p2+cosh(q1*q2)*p3'';',
   'V=vset(q1,q2,uncertStruct,p0,p1,p2,p3,j*(0:.4:4));vsetplot(V)',
   '' };
  slide(29).text={
   'Another example:',
   '',
   '>> uncertStruct=''p0+q2*p1+sinh(q2)*p2+cosh(q1*q2)*p3'';',
   '>> V=vset(q1,q2,uncertStruct,p0,p1,p2,p3,j*(0:.4:4));',
   '>> vsetplot(V)',
   ''};

  %========== Slide 30 ==========

  slide(30).code={
   'uncertStruct=''p0+airy(q1)*p1+airy(q2)*p2+(q1*q2)*p3'';',
   'V=vset(q1,q2,uncertStruct,p0,p1,p2,p3,j*(0:.4:4));vsetplot(V)',
   '' };
  slide(30).text={
   'And the concluding example:',
   '',
   '>> uncertStruct=''p0+airy(q1)*p1+airy(q2)*p2+(q1*q2)*p3'';',
   '>> V=vset(q1,q2,uncertStruct,p0,p1,p2,p3,j*(0:.4:4));vsetplot(V)',
   ''};

  %========== Slide 31 ==========

  slide(31).code={
   'cla;',
   'set(gca,''Visible'',''off'');',
   '' };
  slide(31).text={
   '9. SPHERICAL POLYNOMIAL FAMILIES',
   '===================================',
   'If the admissible set of uncertain parameters is defined by bounds on L2 norm,',
   'such polynomials are called spherical uncertain polynomials. SPHERPLOT function is a tool for computing and plotting the value sets.'};

  %========== Slide 32 ==========

  slide(32).code={
   'cla;',
   'set(gca,''Visible'',''off'');' };
  slide(32).text={
   'Consider the polynomial family p(s,q) = p0(s) + q1*s + q2*s + q3*s,  where',
   '',
   'p0(s) = 0.5 + s + 2*s^2 + 4*s^3',
   '',
   'and the weighted L2 norm of the vector of uncertain parameters is less than 1',
   'The weighting matrix for theEuclidean  norm is given by   W = diag{2,5,3,1}.',
   '',
   'Is the polynomial family robustly stable?',
   ''};

  %========== Slide 33 ==========

  slide(33).code={
   'p0 = 0.5 + s + 2*s^2 + 4*s^3;',
   'weight = [2 5 3 1];',
   'bound = 1;',
   'omega = 0:.025:1.5;',
   'spherplot(p0,omega,bound,weight);',
   '',
   '',
   '' };
  slide(33).text={
   'First enter the nominal polynomial, vector of diagonal entries of the weighting matrix, the bound on the L2 norm and the frequency range:',
   '',
   '>> p0 = 0.5 + s + 2*s^2 + 4*s^3;  weight = [2 5 3 1];',
   '>> bound = 1;  omega = 0:.025:1.5;',
   '',
   'The value sets are obtained with the SPHERPLOT function:',
   '',
   '>> spherplot(p0,omega,bound,weight);',
   ''};
end

%end .. polrobustshow
