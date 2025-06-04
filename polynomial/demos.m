function tbxStruct=Demos
% DEMOS   Demo List information for Polynomial Toolbox

% Author(s):    Z. Hurak   12-7-2000
% Copyright (c) 2000 by PolyX, Ltd.
% $Revision: 1.0.3 $  $Date: 01-10-2000 $
%                     $Date: 04-07-2001  J.Jezek, cosmetics $ 

if nargout==0, demo toolbox; return; end

tbxStruct.Name='Polynomial';
tbxStruct.Type='toolbox';

tbxStruct.Help={
         ' The Polynomial Toolbox 3.0 is a package for systems,         '
         ' signals and control analysis and design based on advanced    '
         ' polynomial methods.                                          '  
         ' It features:                                                 '
         '                                                              '
         ' - object-oriented description of polynomials and polynomial  '
         '   matrices                                                   '
         ' - polynomial equation solvers                                '
         ' - new generation of numerical algorithms based on FFT: easy, '
         '   fast, reliable                                             '
         ' - system and signal models based on polynomial matrix        '
         '   fractions                                                  '
         ' - analysis and design tools for control  and  filters        '
         ' - classical, optimal and robust design:                      '
         '   deadbeat, pole-placement, LQG, H2, Hinfinity               '};

tbxStruct.DemoList={
   'TUTORIALS','',
   '   Getting started   [GUI]', 'playshow poltutorshow'
   'CONTROL DESIGN','',
   '   Deadbeat compensator   [command line]','poldemodebe',
   '   Parametric robust analysis   [GUI]','playshow polrobustshow',
   '   Mixed sensitivity Hinf design [command line]','poldemomixsens';
   'NUMERICAL VS. SYMBOLIC COMPUTATION','',
   '   Determinant of a polynomial matrix   [command line]','poldemodet'};

%end .. Demos



