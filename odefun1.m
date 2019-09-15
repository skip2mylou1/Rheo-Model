function dAdt = odefun1(t, A)
% Loading the configuration tensor A(6x1), then use the scalar form of the
% tensor to calculate the dA/dt. Next, this relationship is used to solve
%  the ODE and obtain A as a function of t

dAdt = zeros(6,1);
global Rate
% (Shear) Rate is a global variable. It is a parameter of the model that
% needs to be studied.

% Below are different parameters of the polymer based on the "tube" model of
% polymer melts.
tD = 3.206; tR = 0.0697; Bccr = 1.25; delta = -0.5;


% Using the scalar form of the matrix A to conduct calculation on each
% independent variables of the matrix.

% dAdt(1)=a11; dAdt(2)=a22; dAdt(3)=a33; 
% dAdt(4)=a12; dAdt(5)=a23; dAdt(6)=a13;
dAdt(1) = 2 * Rate * A(4) - 1 / tD * (A(1) - 1) - 2 / tR * (1 - sqrt(3 / (A(1) + A(2) + A(3)))) * (A(1) + Bccr * ((A(1) + A(2) + A(3)) / 3)^delta * (A(1) - 1));
dAdt(2) = -1 / tD * (A(2) - 1) - 2 / tR * (1 - sqrt(3 / (A(1) + A(2) + A(3)))) * (A(2) + Bccr * ((A(1) + A(2) + A(3)) / 3)^delta * (A(2) - 1));
dAdt(3) = -1 / tD * (A(3) - 1) - 2 / tR * (1 - sqrt(3 / (A(1) + A(2) + A(3)))) * (A(3) + Bccr * ((A(1) + A(2) + A(3)) / 3)^delta * (A(3) - 1));
dAdt(4) = Rate * A(2) - 1 / tD * A(4) - 2 / tR * (1 - sqrt(3 / (A(1) + A(2) + A(3)))) * (A(4) + Bccr * ((A(1) + A(2) + A(3)) / 3)^delta * A(4));
dAdt(5) = -1 / tD * A(5) - 2 / tR * (1 - sqrt(3 / (A(1) + A(2) + A(3)))) * (A(5) + Bccr * ((A(1) + A(2) + A(3)) / 3)^delta * A(5));
dAdt(6) = Rate * A(5) - 1 / tD * A(6) - 2 / tR * (1 - sqrt(3 / (A(1) + A(2) + A(3)))) * (A(6) + Bccr * ((A(1) + A(2) + A(3)) / 3)^delta * A(6));   
end
