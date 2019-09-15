% In Startup Shear, the constitutive equation is used to calculate the
% stress as a function of time under different shear rate
tic;
clear all;
% G is the modulus of the polymer, the unit is Pa. This is an adjustable parameter.
G = 320.33;

% (Shear) Rate is a global variable. It is a parameter of the model that
% needs to be studied.
global Rate

% create a series of shear rate ranging from 10^-3 ~ 10^3 s-1
for i = 1 : 31
    Rate = 10^(-3 + 0.2*(i - 1));
    rate(i) = Rate;     %store the shear rate
    A0 = [1;1;1;0;0;0];     %initialize the configuration to isotropic state
    % created the time range used in ODE
    tspan = [0 10^(3 - 0.1*(i - 1))];
    % Use ODE solver to calculate A as a function of t
    [t,A] = ode45('odefun1',tspan,A0);
    sigma = G * A(:,4); % sigma = stress
    N1 = G * (A(:,1) - A(:,2)); % N1 is the first normal stress difference
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
    
    %for shear stress
    [sigma_max(i),POS1] = max(sigma);
    t_max1(i) = t(POS1);    % sigma_max ~ t_max
    strain_max1(i) = t_max1(i) * rate(i);
    sigma_max_steady(i) = sigma_max(i) / sigma(steady);
    %for normal stress difference
    [N1_max(i),POS2] = max(N1);
    t_max2(i) = t(POS2);
    strain_max2(i) = t_max2(i) * rate(i);
    N1_max_steady(i) = N1_max(i) / N1(steady);
end;
% save all the files above in start_up_shear_result.csv
columns = {'t_max1', 'sigma_max', 't_max2', 'N1_max', 'rate', 'strain_max1', 'strain_max2', 'sigma_max_steady', 'N1_max_steady'};
data = table(t_max1', sigma_max', t_max2', N1_max', rate', strain_max1', strain_max2', sigma_max_steady', N1_max_steady','VariableNames', columns);
writetable(data, 'start_up_shear_result.csv')

toc;

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
figure(36);loglog(strain_max1,sigma_max,'o');   %result not good
xlabel('strain max1') 
ylabel('sigma max') 
title('sigma max ~ strain max1')
figure(37);loglog(strain_max2,N1_max,'o');  %result not good
xlabel('strain max2') 
ylabel('N1 max') 
title('N1 max ~ strain max2')
figure(38);loglog(rate,strain_max1,'o');    %result not good
xlabel('rate') 
xlabel('strain max1')
title('strain max1 ~ rate')
figure(39);loglog(rate,strain_max2,'o');    %result not good
xlabel('rate') 
xlabel('strain max2')
title('strain max2 ~ rate')
figure(40);loglog(rate,sigma_max_steady,'o');
xlabel('rate') 
xlabel('sigma max steady')
title('sigma max steady ~ rate')
figure(41);loglog(rate,N1_max_steady,'o');
xlabel('rate') 
xlabel('N1 max_steady')
title('N1 max steady ~ rate')










