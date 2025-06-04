% Polynomial Toolbox.
% Version 3.0.19 10-May-2003
%
% Global Structure
%
%   pinit           - Initialize the Polynomial Toolbox
%   pver            - Polynomial Toolbox version information
%   pversion        - Polynomial Toolbox version number
%   gprops          - Display/modify global properties
%   pformat         - Modify output format
%   tolerance       - Modify global relative tolerance
%   gensym          - Modify global variable symbol
%   verbose         - Modify global verbose level
%   checkpb         - Check conflicts of Polynomial Toolbox
%
% Polynomial(pol), Two-sided Polynomial(tsp) and Fraction Objects:
% Left-denominator Fractions(ldf), Right-denominator Fractions(rdf),
% Matrix-denominator Fractions(mdf) and Scalar-denominator
% Fractions(sdf)
%
%   props           - Display/modify properties of object
%   symbol          - Display/modify variable symbol of object
%   userdata        - Display/modify user data of object
%   pol             - Create polynomial
%   lop             - Create polynomial
%   tsp             - Create two-sided polynomial
%   rdf             - Create right-denominator fraction
%   ldf             - Create left-denominator fraction
%   mdf             - Create matrix-denominator fraction
%   sdf             - Create scalar-denominator fraction
%   size            - Size
%   length          - Length
%   deg             - Degree of pol or tsp
%   tdeg            - Trailing degree of pol or tsp
%   coef            - Coefficient of pol or tsp
%   lcoef           - Leading coefficient of pol or tsp
%   tcoef           - Trailing coefficient of pol or tsp
%
% Special matrices
%
%   s               - Basic variable (derivative operator)
%   p               - Basic variable (derivative operator)
%   z               - Basic variable (advance shift operator)
%   q               - Basic variable (advance shift operator)
%   zi              - Basic variable (delay shift operator)
%   d               - Basic variable (delay shift operator)
%   v               - Basic variable
%   mono            - Monomial
%
% Convertors
%
%   double          - Convert to double (standard Matlab matrix)
%   pol             - Convert to polynomial
%   tsp             - Convert to two-sided polynomial
%   rdf             - Convert to right-den fraction
%   ldf             - Convert to left-den fraction
%   mdf             - Convert to matrix-den fraction
%   sdf             - Convert to scalar-den fraction
%   coldeg2pol      - Convert column degree expanded form into the POL object.
%   pol2coldeg      - Decomposition acording to column degrees
%   pol2rowdeg      - Decomposition acording to row degrees
%   rowdeg2pol      - Convert row degree expanded form into the POL object.
%   declass         - Declass (convert to possible lower class)
%   defract         - Remove denominator if possible
%   reverse         - Reverse variable of fraction
%   abcd            - Convert fraction to state space
%   abcde           - Convert fraction to descriptor state space
%   ss              - Convert fraction to Control Toolbox state space
%   dss             - Convert fraction to Control Toolbox descriptor state space
%   tf              - Convert fraction to Control Toolbox transfer function
%   zpk             - Convert fraction to Control Toolbox zero-pole-gain
%   sym             - Convert to Symbolic Toolbox
%   
% Overloaded operations
%
%   uplus           - (+)   Unary plus
%   uminus          - (-)   Unary minus
%   plus            - (+)   Add
%   minus           - (-)   Subtract
%   times           - (.*)  Element-wise multiply
%   mtimes          - (*)   Matrix multiply
%   power           - (.^)  Element-wise power
%   mpower          - (^)   Matrix power
%   rdivide         - (./)  Element-wise right divide
%   ldivide         - (.\)  Element-wise left divide
%   mrdivide        - (/)   Matrix right divide
%   mldivide        - (\)   Matrix left divide
%   kron            -       Kronecker tensor product
%   transpose       - (.')  Transpose
%   ctranspose      - (')   Conjugated transpose
%   eq              - (==)  Test if equal
%   ne              - (~=)  Test if unequal
%   display         -       Display
%   char            -       Convert to character string
%   horzcat         - ([,]) Horizontal concatenate
%   vertcat         - ([;]) Vertical concatenate
%   cat             -       Concatenate
%   subsref         -       Subscripted reference
%   subsasgn        -       Subscripted assignment
%
% Overloaded functions
%
%   sum             - Element-wise sum
%   prod            - Element-wise product
%   inv             - Inverse
%   pinv            - Pseudoinverse
%   real            - Real part
%   imag            - Imaginary part
%   conj            - Complex conjugate
%   diag            - Diagonal matrix or diagonal
%   trace           - Sum of diagonal elements
%   triu            - Upper triangular part
%   tril            - Lower triangular part
%   det             - Determinant
%   adj             - Adjoint of pol or tsp
%   isempty         - Test if empty
%   isreal          - Test if real
%   isfinite        - Test if pol or tsp is finite
%   isinf           - Test if pol or tsp is infinite
%   isnan           - Test if pol or tsp is Not-a-Number
%   isprime         - Test if pol matrix is left or right prime
%   rank            - Rank
%   null            - Nullspace of pol, tsp or fraction
%   norm            - Norm
%   repmat          - Replicate and tile
%   reshape         - Change size
%   flipud          - Flip rows up-down
%   fliplr          - Flip columns left-right
%   rot90           - Rotate 90 degrees
%   roots           - Roots
%   compan          - Block companion matrix
%   lu              - LU factorization of pol matrix
%   sylv            - Sylvester matrix
%
% Basic functions (other than overloaded)
%
%   pzer            - Perform zeroing on a polynomial matrix
%   value           - Value
%   mvalue          - Matrix value
%   rdiv            - Right divide with remainder
%   ldiv            - Left divide with remainder
%   pgcd            - Greatest common divisor of pols
%   gld             - Greatest left divisor of pol matrices
%   grd             - Greatest right divisor of pol matrices
%   plcm            - Least common multiple of pols
%   llm             - Least left multiple of pol matrices
%   lrm             - Least right multiple of pol matrices
%   fact            - Factor exctraction
%   deriv           - Derivative
%   integral        - Integral of pol or tsp
%   iscolred        - Test for column reducedness of a polynomial matrix
%   isfullrank      - Test if full rank
%   isrowred        - Test for row reducedness of a polynomial matrix
%   issingular      - Test if singula
%   isproper        - Test if fraction is proper
%   evenpart        - Even part of pol
%   oddpart         - Odd part of pol
%   mirror          - Mirror image of pol or tsp
%   pos             - Positive part of tsp
%   neg             - Negative part of tsp
%   npos            - Nonpositive part of tsp
%   nneg            - Nonnegative part of tsp
%   shift           - Shift pol or tsp
%   scale           - Scale pol
%   linvt           - Linear transform of variable
%   coprime         - Make fraction coprime
%   longldiv        - Long left divide
%   longrdiv        - Long right divide
%   laurent         - Laurent series of fraction
%   polfit          - Fit pol matrix element-by-element to data
%   polpart         - Polynomial matrix part extraction
%   polyvalm        - Evaluation of a polynomial in a matrix
%
% Random functions
%
%   prand           - Random pol
%   trand           - Random tsp
%   rrand           - Random rdf
%   lrand           - Random ldf
%   mrand           - Random mdf
%   srand           - Random sdf
%
% Advanced functions
%
%   charact         - Characteristic vectors of pol matrix
%   complete        - Complete to unimodular
%   ellsta          - Ellipsoidal approximation of the stability domain
%   h2norm          - H2 norm of fraction
%   hinfnorm        - H-infinity norm of fraction
%   hurwitz         - Create Hurwitz matrix of polynomial
%   inertia         - Inertia of polynomial matrix
%   isminph         - Test if minimum phase
%   isstable        - Test if stable
%   isunimod        - Test if unimodular
%   laplace         - Inverse Laplace transform of fraction
%   minbasis        - Minimal polynomial basis
%   moments         - Moment of roots of a polynomial
%   jury            - Jury matrix 
%   kharit          - Kharitonov polynomial
%   gram            - Gramian of polynomial matrix fraction
%
% Canonical and reduced forms
%
%   rowred          - Row reduced form
%   colred          - Column reduced form
%   diagred         - Diagonal reduced form
%   pdg             - Diagonal form
%   tri             - Triangular or staircase form
%   hermite         - Hermite form
%   smith           - Smith form
%   echelon         - Echelon form
%   bhf             - Block Hessenberg form
%   bezout          - Bezoutian matrix
%   hermfuji        - Hermite-Fujiwara matrix
%   schurcohn       - Schur-Cohn matrix
%   reduce          - Make fraction reduced
%   
% Sampling period functions
%   
%   samp            - Sample continuous-time object
%   samph           - Sample and zero order hold
%   samph1          - Sample and first order hold
%   samph2          - Sample and second order hold
%   resamp          - Resample discrete-time object
%   resamph         - Resample and zero order hold
%   resamph1        - Resample and first order hole
%   resamph2        - Resample and second order hold
%   unsamp          - Unsample discrete-time object
%   unsamph         - Un- sample and zero order hold
%   unsamph1        - Un- sample and first order hold
%   unsamph2        - Un- sample and second order hold
%   chsamp          - Change sampling of discrete-time object
%   dilate          - Dilate discrete-time object
%
% Equation solvers
%
%   axb             - Solution of  AX = B
%   xab             - Solution of  XA = B
%   axbc            - Solution of  AXB = C
%   axby0           - Solution of  AX + BY = 0
%   xayb0           - Solution of  XA + YB = 0
%   axbyc           - Solution of  AX + BY = C
%   xaybc           - Solution of  XA + YB = C
%   axybc           - Solution of  AX + YB = C
%   axxab           - Solution of  A'X + X'A = B
%   xaaxb           - Solution of  XA' + AX' = B
%   axyab           - Solution of  AX' + Y'A = B
%   hqr             - Hyperbolic QR factorization Q*A = R of a square matrix A 
%   spf             - Spectral factorization     X'JX = B
%   spcof           - Spectral co-factorization  XJX' = B
%   gare            - Generalized algebraic Riccati equation
%
% Linear-matrix-inequality functions
%
%   lmianalysis     - Analysis of stability of pol matrix
%   lmihinfnorm     - H-inf norm of fraction
%   lmipolytope     - Robust stability analysis of polytope
%   lmirank         - Solve LMI for a matrix of specified rank
%   lmisimstab      - Simultaneous stabilization of family of SISO systems
%   lmispf          - Spectral factorization
%   lmispcof        - Spectral co-factorization
%
% Matrix pencil functions
%
%   clements        - Conversion to Clements standard form
%   clements1       - Conversion to Clements standard form
%   pencan          - Conversion to real Kronecker canonical form
%   plyap           - Solution of the pencil equation  A*X + Y*B = C
%   
% Numerical functions
%
%   cgivens1        - Calculates Givens rotation
%   qzord           - Ordered QZ transformation
%   schurst         - Ordered complex Schur decomposition of a matrix

% Control functions
%
%   debe            - Deadbeat controllers of discrete-time linear systems
%   dssh2           - Descriptor solution of the H2 problem
%   dsshinf         - H-inf suboptimal compensator for descriptor systems
%   dssmin          - Minimize dimension of pseudo state descriptor system
%   dssrch          - Search Optimal Solution descriptor H-inf problem
%   dssreg          - "Regularizes" a standard descriptor plant
%   h2              - H2-optimization
%   mixeds          - Solution SISO mixed sensitivity problem
%   plqg            - Polynomial solution of a MIMO LQG problem
%   pplace          - Polynomial pole placement
%   psseig          - Polynomial approach to eigenstructure assignment 
%                     for state-space system
%   psslqr          - Polynomial approach to linear-quadratic regulator 
%                     design for state-space system
%   splqg           - Polynomial solution of a SISO LQG problem
%   stab            - Stabilizing controllers of linear systems
%
% Robustness functions
%
%   edgetest        - Test for robust stability of a polytope of polynomials
%   isdiscstable    - Test for robust stability of a disc polynomial
%   stabint         - Stability interval of a polynomial matrix
%   sarea           - Stability test for a family of polynomials
%                     with parametric uncertainties
%   tsyp            - Robustness margin for a continuous interval polynomial
%   vset            - Value set of parametric polynomial
%   ptopex          - Extreme polynomials for a polytope of polynomials
%
% 2-d functions
%
%   det2d           - Determinant of 2-d polynomial matrix
%
% Simulink
%
%   polblock        - Simulink mdl-file
%   polymask        - Initialize mask of Polynomial Library
%   update          - Update polynomial library block
%
% Visualization
%
%   pdisp           - Display of a polynomial matrix without printing the name
%   picture         - Graphic picture of polynomial or fraction
%   khplot          - Plot of Kharitonov rectangles for interval polynomials
%   pplot           - 2-D plot of polynomial matrix
%   pplot3          - 3-D plot of polynomial matrix
%   ptopplot        - Plot polygonal values set for polytope of polynomials
%   sareaplot       - Plot stability area 
%   spherplot       - Plot the value set ellipses for a spherical polynomial family
%   sphspectrum     - Plot the spectral set for a spherical polynomial family
%   vsetplot        - Plots value set for a parametric polynomial
%   zpplot          - Plot of zero-pole map
%
% Graphic User Interface
%
%   pme             - Polynomial Matrix Editor
%
% Demonstrations and helps
%
%   covf            - Covariance function of an ARMA process
%   demoB           - Script file for the demo "Control of a batch process"
%   demoM           - Script file for the demo "Polynomial solution of the SISO
%		                mixed sensitivity H-infinity problem"
%   demos           - Demo List information for Polynomial Toolbox
%   minsens         - Minimum peak sensitivity
%   poldemo         - Run Polynomial Toolbox demonstrations
%   poldemodebe     - Run demonstrations of deadbeat compensator design
%   poldemodet      - Comparision of numerical and symbolic computation 
%                     of polynomial matrix determinant
%   poldemomixsens  - Mixed sensitivity
%   poldemorobpar   - Stability analysis of uncertain systems
%   poldesk         - Comprehensive hypertext documentation
%   polrobustshow   - Demo introduces some analytical tools for systems 
%                     with uncertain parameters
%   poltutorshow    - Getting started with Polynomial Toolbox 2.5

%
% Copyright(c) 2003 by Polyx, Ltd.

