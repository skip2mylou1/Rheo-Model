% In Startup Shear, the constitutive equation is used to calculate the
% stress as a function of time under different shear rate
clear all;
% G is the modulus of the polymer, the unit is Pa. This is an adjustable parameter.
G = 320.33;

% (Shear) Rate is a global variable. It is a parameter of the model that
% needs to be studied.
% global Rate
% global Rate_matrix
% create a series of shear rate ranging from 10^-3 ~ 10^3 s-1
i_number = 31;
sigma_max = zeros(1,i_number);
t_max1 = zeros(1,i_number);   
strain_max1 = zeros(1,i_number);
sigma_max_steady = zeros(1,i_number);
N1_max = zeros(1,i_number);
t_max2 = zeros(1,i_number); 
strain_max2 = zeros(1,i_number);
N1_max_steady = zeros(1,i_number);
A_ss = zeros(i_number,9);
rate = zeros(1,i_number);
SR = zeros(1,i_number); 
tic;
for i = 1 : i_number
    Rate = 10^(-3 + 0.2*(i - 1));
    rate(i) = Rate;     %store the shear rate
%     matrix of shear rate, for simple steady shear, the result is as below
%     Rate_matrix = [[0,0,0];[Rate,0,0];[0,0,0]]
    % created the time range used in ODE
    tspan = [0 10^(3 - 0.1*(i - 1))];
    % Use ODE solver to calculate A as a function of t
    A0 = [1;0;0;0;1;0;0;0;1];
    [t,A] = ode15s(@(t, A) odefunmatrix(t, A, Rate),tspan,A0);
    sigma = G * A(:,4); % sigma = stress
    N1 = G * (A(:,1) - A(:,5)); % N1 is the first normal stress difference
    [t_steady,steady] = max(t);

    %ploting the shear stress and normal stress
    figure(i);
    subplot(1,2,1);
    plot(t,sigma);
    title(['Shear Stress ~ t (' num2str(i) ')'] )
    
    subplot(1,2,2);
    plot(t,N1);
    title(['Normal Stress Difference ~ t ' num2str(i) ')'])
    %the following results are documented to make a plot in order to find
    %the scaling of properties
    
    %shear stress
    [sigma_max(i),POS1] = max(sigma);
    t_max1(i) = t(POS1);    % sigma_max ~ t_max
    strain_max1(i) = t_max1(i) * rate(i);
    sigma_max_steady(i) = sigma_max(i) / sigma(steady);
    
    %normal stress difference
    [N1_max(i),POS2] = max(N1);
    t_max2(i) = t(POS2);
    strain_max2(i) = t_max2(i) * rate(i);
    N1_max_steady(i) = N1_max(i) / N1(steady);
    
%%  stiffness of the ode and matrix
    A_ss(i,:) = A(end,:);
    A_s = reshape(A(end,:),3,3);
    e = eig(A_s);
    SR(i) = max(abs(e))/min(abs(e));

end;
toc;
% save all the files above in start_up_shear_result.csv
columns = {'t_max1', 'sigma_max', 't_max2', 'N1_max', 'rate', 'strain_max1', 'strain_max2', 'sigma_max_steady', 'N1_max_steady'};
data = table(t_max1', sigma_max', t_max2', N1_max', rate', strain_max1', strain_max2', sigma_max_steady', N1_max_steady','VariableNames', columns);
writetable(data, 'start_up_shear_result.csv')



% ploting the results
figure(32);loglog(t_max1,sigma_max,'o');
xlabel('t max1') 
ylabel('sigma max') 
title('sigma max ~ t max1')
figure(33);loglog(t_max2,N1_max,'o');   %result not good
xlabel('t max2') 
ylabel('N1 max') 
title('N1 max ~ t max2')
figure(34);loglog(rate,sigma_max,'o');
xlabel('rate') 
ylabel('sigma max') 
title('sigma max ~ rate')
figure(35);loglog(rate,N1_max,'o');
xlabel('rate') 
ylabel('N1 max') 
title('N1 max ~ rate')
% figure(36);loglog(strain_max1,sigma_max,'o');   %result not good
% xlabel('strain max1') 
% ylabel('sigma max') 
% title('sigma max ~ strain max1')
% figure(37);loglog(strain_max2,N1_max,'o');  %result not good
% xlabel('strain max2') 
% ylabel('N1 max') 
% title('N1 max ~ strain max2')
% figure(38);loglog(rate,strain_max1,'o');    %result not good
% xlabel('rate') 
% xlabel('strain max1')
% title('strain max1 ~ rate')
% figure(39);loglog(rate,strain_max2,'o');    %result not good
% xlabel('rate') 
% xlabel('strain max2')
% title('strain max2 ~ rate')
% figure(40);loglog(rate,sigma_max_steady,'o');
% xlabel('rate') 
% xlabel('sigma max steady')
% title('sigma max steady ~ rate')
% figure(41);loglog(rate,N1_max_steady,'o');
% xlabel('rate') 
% xlabel('N1 max_steady')
% title('N1 max steady ~ rate')







