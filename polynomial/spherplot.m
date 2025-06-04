function spherplot(p0,omega,r,weight)
%SPHERPLOT   Plot the value set ellipses for a spherical polynomial family. 
%
%  SPHERPLOT(P0,OMEGA,R,WEIGHT) Plot the value set ellipses for a spherical  
%  polynomial family given by 
%
%  p(s,q)=p0(s)+ q0 + q1*s + q2*s^2 + ... + qk*s^k 
%
%  where the vector of uncertain parameters q0, q1, ..., qk is subject 
%  to norm(q,2,WEIGHT) <= R, which reads that the weighted Euclidian norm 
%  of the vector of uncertain parameters is bounded by the robustness margin R.
%  WEIGHT is a vector of diagonal entries of the weighting matrix, OMEGA is 
%  a vector of frequencies for which the value set ellipses are computed and plotted.
%  
%  See also: KHPLOT, PTOPPLOT, VSETPLOT

%  See B.R.Barmish: New Tools for Robustness of Linear Systems, Macmillan,1994, pp.270.

%  Author: Zdenek Hurak 04-26-2000
%  Copyright (c) 2000 by Polyx, Ltd.
%  $ Revision: 2.5.1 $     $ Date: 12-13-2000 $
%  $ Revision: 3.0.0 $     $ Date: 08-19-2002 $

global PGLOBAL;
try,
    PGLOBAL.FORMAT;
catch,
    painit
end

%-------------------------------------------------------------------------------------
% Input tests:
%-------------------------------------------------------------------------------------

if nargin < 2
   error('Not enough input arguments; at least two required.');
end

try,
    p0 = pol(p0);
catch,
    error(peel(lasterr));
end

if ~isnumeric(omega) | ndims(omega)>2 | ~any(size(omega)==1),
   error('Invalid 2nd argument; must be a nonempty numerical vector.');
end

if nargin<3 | isempty(r),
   r = 1;
else
   if ~isnumeric(r) | length(r)~=1 | ~isreal(r) | r<0,
      error('Invalid 3rd argument; must be a nonnegative number.');
   end;
end;

if nargin<4 | isempty(weight),
   weight = ones(1,p0.deg+1);
else
   if ~isnumeric(weight) | ~any(size(weight)==1),
      error('Invalid 4th argument; must be a numerical vector.')
   end
   [rows,columns] = size(weight);
   if rows*columns ~= (p0.deg+1)
      error('Invalid size of the vector of weighting factors.');
   end
end

%-------------------------------------------------------------------------------------
% Computation of weighted sums of (doubled) powers of frequencies (even and odd):
%-------------------------------------------------------------------------------------
Mon = sum(mono(0:p0.deg));      %creates the polynomial 1+x+x^2+x^3+x^4+...+x^n
reMon = evenpart(Mon);			%separates the polynomial into even
imMon = oddpart(Mon);           %and odd parts
reK = 2*find(reMon{:}>0)-2;     %doubles the powers of the even part (real part)
imK = 2*find(imMon{:}>0)-2;     %and the odd part (imaginary part)
                                %we thus have two arrays of polynomials: [1 x^4 x^8 ...]
                                %and [x^2 x^6 x^10 ...]

wEven = weight(1:2:length(weight));       %separates the vector of weighting coefficients into
wOdd = weight(2:2:length(weight));        %even and odd parts 

reMon2w = sum(wEven.*mono(reK));%performs scaling with the scaling factors WEIGHT and puts 
imMon2w = sum(wOdd.*mono(imK)); %together to get the two polynomials: weight0+weight2*x^4+
                                %+weight4*x^8+... 
                                %and weight1*x^2+weight3*x^6+...
                                %Note that the elements of the vector WEIGHT are already
                                %squared, i.e. they are actually the entries of the 
                                %weighting matrix. So it is not necessary to square above.

S_even = squeeze(polyval(reMon2w,omega)); %evaluates them at the frequencies OMEGA
S_odd = squeeze(polyval(imMon2w,omega));

%-------------------------------------------------------------------------------------
% Computation of the major axes of the value set ellipses (in Re and Im directions):
%-------------------------------------------------------------------------------------
R0 = r*sqrt(S_even); %length of half the major axis of the ellipse in the Re direction
I0 = r*sqrt(S_odd);  %length of half the major axis of the ellipse in the Im direction


jj = 0:.01:1;          
circ = exp(j*2*pi*jj); %computes the ring (a two-dimensional sphere in a complex plane)
reCirc = real(circ);
imCirc = imag(circ);
reCirc = repmat(reCirc,length(omega),1);
imCirc = repmat(imCirc,length(omega),1);

R0 = repmat(R0,1,length(jj)); 
I0 = repmat(I0,1,length(jj)); 

reElip = reCirc .* R0; %scales the ring with R0 at all frequencies in Re dir.
imElip = imCirc .* I0; %scales the ring with I0 at all frequencies in Im dir.


elipse = reElip + j*imElip;
shiftElipse = elipse + repmat(squeeze(polyval(p0,j*omega)),1,length(jj));
%positiones the ellipse at the place given by the nominal p0

Rmax = max(max(abs(real(shiftElipse))));
Imax = max(max(abs(imag(shiftElipse))));

plot(shiftElipse.');
hold on;
plot([-Rmax,Rmax,NaN,0,0],[0,0,NaN,-Imax,Imax],':k');
hold off;
axis([-Rmax Rmax -Imax Imax]);
title('Value set of a spherical polynomial family');
xlabel('Real');
ylabel('Imaginary');

% end .. spherplot







