% POLDEMO Script to demonstrate some of the polynomial toolbox oparations.
%
% Type
%    poldemo 
% to run the script.

% Huibert Kwakernaak, January 5, 2000
% Copyright PolyX Ltd, 2000
% Modified by:
%       M. Hromcik  09-Aug-05:  call inv for V3

echo off
pinit
disp('This demo reviews several of the functions and operations defined')
disp('in the Polynomial Toolbox for polynomials and polynomial matrices.') 
disp(' ')
disp('Detailed information is available in the Manual.')
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')

disp('The following topics are covered:')
disp(' ')
disp('   Entering polynomials and polynomial matrices')
disp('   Concatenation and submatrices')
disp('   Extracting coefficients')
disp('   Algebraic operations: addition, subtraction, multiplication')
disp('   Determinants')
disp('   Unimodular polynomial matrices')
disp('   Adjoint and inverses')
disp('   Rank')
disp('   Null space and range of a polynomial matrix')
disp('   Finite and infinite roots')
disp('   Stability of polynomials and polynomial matrices')
disp('   Division of polynomials and polynomial matrices')
disp('  ')
disp('See the Tutorial and Manual for the many operations not')
disp('included in this demo.')
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')

disp('Polynomials and polynomial matrices are most easily entered using')
disp('one of the indeterminate variables s, p, z, q, z^-1 or d that are')
disp('recognized by the Polynomial Toolbox, combined with the usual MATLAB')
disp('conventions for entering matrices.')
disp(' ')
disp('Thus you may easily define a simple polynomial matrix by typing')
disp(' ')
echo on
P = [ 1+s 2*s^2; 2+s^3 4 ]
echo off 
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')

disp('In a system and control theoretic context, the indeterminates')
disp('s and p are usually associated with continuous time, and are ')
disp('closely related to the differential operator x(t) -> dx(t)/dt.')
disp('The indeterminates z and q are associated with the discrete-time ')
disp('one-step shift operator x(t) ->  x(t+1), and the indeterminates ')
disp('d and z^-1 with the delay operator x(t) ->  x(t-1).')
disp(' ')
disp('The indeterminate z^-1 may be entered with arbitary nonnegative ')
disp('integral powers. For typing convenience z^-1 may be abbreviated to zi.')
disp('Thus, the following  two commands return the same result:')
disp(' ')
echo on
P = 1 + z^-1 + z^-2

P = 1 + zi + zi^2
echo off 
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')

disp('Polynomials and polynomial matrices may also be entered in terms ')
disp('of their coefficients or coefficient matrices. For this purpose the ')
disp('pol command is available.')
disp(' ')
disp('With this command we may for instance define the following polynomial ')
disp('matrix')
disp(' ')
echo on
P0 = [1 2; 3 4];
P1 = [3 4; 5 1];
P2 = [1 0; 0 1];
P = pol([P0 P1 P2],2,'s')
echo off 
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')
disp('After the Polynomial Toolbox has been started up the default ')
disp('indeterminate variable is s. This implies, among other things, that ')
disp('the last command may also be entered as')
disp(' ')
echo on
P = pol([P0 P1 P2],2)
echo off 
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')

disp('Standard MATLAB conventions may be used to concatenate polynomial ')
disp('and standard matrices. Here is an example:')
disp(' ')
echo on
Q = [P; 1+s 3]
echo off 
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')
disp('Submatrices may be selected in the usual way:')
disp(' ')
echo on
Q(2:3,:)
echo off 
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')

disp('It is easy to extract the coefficient matrices of a polynomial ')
disp('matrix. Consider the polynomial matrix')
disp(' ')
echo on
T = [ 1+s 2 s^2 s; 3 4 s^3 0 ]
echo off 
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')
disp('The coefficient matrix of s^2 may be retrieved as')
disp(' ')
echo on
T{2}
echo off 
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')
disp('The coefficients of the (1,3) entry of T follow as')
disp(' ')
echo on
T{:}(1,3)
echo off 
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')

disp('Define the polynomial matrices')
disp(' ')
echo on
P = [ 1+s 2; 3 4], Q = [ s^2 s; s^3 0]
echo off 
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')
disp('The sum and product of P and Q follow easily:')
disp(' ')
echo on
S = P+Q

R = P*Q
echo off 
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')

disp('The determinant of a square polynomial matrix is defined exactly ')
disp('as its constant matrix counterpart. In fact, its computation is not ')
disp('much more difficult:')
disp(' ')
echo on
P = [1 s s^2; 1+s s 1-s; 0 -1 -s]

detP = det(P)
echo off 
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')

disp('If its determinant happens to be constant then the polynomial')
disp('matrix is called unimodular:')
disp(' ')
echo on
U = [ 2-s-2*s^2 2-2*s^2 1+s; 1-s-s^2 1-s^2 s; -1-s -s 1]

det(U)

echo off 
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')

disp('If a matrix is suspected of unimodularity then one may make')
disp('sure by a special tester isunimod:')
disp(' ')
echo on
isunimod(U)
echo off 
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')

disp('Also the adjoint matrix is defined as for constant matrices.')
disp('The adjoint is a polynomial matrix and may be computed by typing')
disp(' ')
echo on
adj(P)
echo off 
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')
disp('The inverse of a polynomial matrix P is given by adj(P)/det(P)')
disp('and, hence, follows by typing')
disp(' ')
echo on
%[adjP,detP] = inv(P)
Pinv = inv(P)
echo off 
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')
disp('As expected')
disp(' ')
echo on
P*adj(P)/det(P)
echo off 
disp('Note that this example involves polynomial division.')
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')

disp('A polynomial matrix P(s) has full column rank (or full normal')
disp('column rank) if it has full column rank everywhere in the complex')
disp('plane except at a finite number of points. Similar definitions hold')
disp('for full row rank and full rank. Recall that')
disp(' ')
echo on
P
echo off 
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')
disp('The following rank test confirms that P has full rank:')
disp(' ')
echo on
isfullrank(P)
echo off 
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')
disp('The normal rank of a polynomial matrix P(s) equals max rank P(s)')
disp('where s ranges over the complex numbers. The rank is calculated by')
disp('typing')
disp(' ')
echo on
rank(P)
echo off 
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')
disp('A polynomial matrix is nonsingular if it has full normal rank:')
disp(' ')
echo on
issingular(P)
echo off
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')

disp('There are two important subspaces (more precisely, submodules)')
disp('associated with a polynomial matrix A(s): its null space and its')
disp('range (or span). The (right) null space is defined as the set of ')
disp('all polynomial vectors x(s) such that A(s)x(s) = 0. It is computed')
disp('by typing')
disp(' ')
echo on
A = P(1:2,:)

N = null(A)
echo off
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')
disp('Hence the null space dimension is 1 and its basis has degree 3.')
disp('Check:')
disp(' ')
echo on
A*N
echo off
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')
disp('The range of A(s) is the set of all polynomial vectors y(s)')
disp('such that y(s) = A(s)x(s) for some polynomial vector x(s).')
disp('In the Polynomial Toolbox, the minimal basis of the range is')
disp('returned by the command')
disp(' ')
echo on
minbasis(A)
echo off
disp('The columns of this matrix form a minimal basis.')
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')

disp('The roots or zeros of a polynomial matrix P(s) are those points')
disp('s in the complex plane where P(s) loses rank:')
disp(' ')
echo on
roots(P)
echo off
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')
disp('The roots can be both finite and infinite. The infinite')
disp('roots are normally suppressed. To reveal them, type')
disp(' ')
echo on
roots(P,'all')
echo off
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')
disp('Unimodular matrices have no finite roots:')
disp(' ')
echo on
roots(U,'all')
echo off
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')

disp('The macro isstable checks stability according to the')
disp('variable symbol. Thus,')
echo off
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')
disp(' ')
echo on
isstable(s-0.5)

isstable(z-0.5)
echo off
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')

disp('To understand when division of polynomials and polynomial ')
disp('matrices is possible, consider three polynomials a(s), b(s)')
disp('and c(s) such that a(s) = b(s)c(s). We say that b(s) is a divisor')
disp('(or factor) of a(s) or a(s) is a multiple of b(s), and write')
disp('a(s)|b(s). This is sometimes also stated as b(s) divides a(s).')
disp('For instance, take')
disp(' ')
echo on
b = 1-s; c = 1+s; a = b*c
echo off
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')

disp('As b(s) and c(s) are both divisors of a(s), the following divisions')
disp('both can be done:')
disp(' ')
echo on
a/b

a/c
echo off
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')
disp('Of course, the division by b(s) fails if b(s) is not a divisor:')
disp(' ')
echo on
c/b
echo off
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')
disp('The quotient and remainder of this division may be retrieved by typing')
disp(' ')
echo on
[q,r] = rdiv(c,b)
echo off
disp(' '); disp('>> Hit any key to continue'); pause; disp(' ')
disp('The Polynomial Toolbox performs these divisions not only')
disp('for polynomials but also for polynomial matrices.')
disp(' ')
disp('See the Tutorial and Manual for the many operations not')
disp('covered in this demo.')

%end .. poldemo

