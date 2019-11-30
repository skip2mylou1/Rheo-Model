clear all;
% create a series of shear rate ranging from 10^-3 ~ 10^3 s-1
% for low frequency, use ode15 to calculate (make it smaller to speed up)
loop_i = 20;  %i is for different shear rate
loop_j = 31;  %j is for different frequency
run_time = zeros(loop_i,1);  %keep track of all the runtime of different rate

G = 320.33;
freq = zeros(1,loop_j);
f = zeros(100,loop_j);
strain = zeros(loop_j,loop_i);
stressAmplitude = zeros(loop_j,loop_i);
I1 = zeros(loop_j,loop_i);
I3 = zeros(loop_j,loop_i);
I31 = zeros(loop_j,loop_i);
rate = zeros(1,loop_i);
nfft = 201;
%% recording the FFT frequency spectrum
for k = 1 : loop_j
        Freq = 10^(-3 + 0.2*(k - 1)); % frequency from 10^-3 ~ 10^-0.4
        freq(k) = Freq;
        Fs = (10*Freq)/pi;
        f(:,k) = (0 : nfft/2 - 1) * Fs / nfft;
end

%% parallel loop of LAOS calculation
parfor i = 1 : loop_i   %the original is 36
    Rate = 10^(-3 + 0.2*(i - 1));
    rate(i) = Rate;
    tic;
    st = zeros(1,loop_j);
    sta = zeros(1,loop_j);
    i31 = zeros(1,loop_j);

    for j = 1 : 14
        Freq = 10^(-3 + 0.2*(j - 1)); % frequency from 10^-3 ~ 10^-0.4
        %freq(j) = Freq;
        st(j) = Rate / Freq;
        %strain(j,i) = Rate / Freq;
        A0 = [1;1;1;0;0;0]; %initialize the configuration to isotropic state
        tspan = [0 : pi/(10*Freq) : 80 * pi / Freq]; %time range for ODE solver
        [t,A] = ode15s(@(t, A) Lodefun_par(t, A, Rate, Freq),tspan,A0);
        sigma = G * A(:,4);
        % parameters for FFT
        Fs = (10*Freq)/pi;   %sampling frequency
        %When shear rate is low and time is small(tspan(1:600)), the
        %system is in SAOS region, the constitutive equation does not work.
        %The result is not accurate.
        ts = tspan(601 : 801);  %time vector use the time range when it is stable
        sigmas = sigma(601 : 801);  %select FFT signal
        sigma_max = max(sigmas);
        sigma_min = min(sigmas);
        sta(j) = (sigma_max - sigma_min) / 2;
        %stressAmplitude(j,i) = (sigma_max - sigma_min) / 2; % amplitude
        %subplot(3,9,2*j-1);
        %plot(ts,sigmas);
        %title('stress signal');xlabel('t');ylabel('stress(Pa)'); 
%% FFT
        nfft = length(sigmas);
        Fsigma = fft(sigmas,nfft);
        Fsigma = Fsigma(1:nfft/2);
        Mag = abs(Fsigma);
        %f(:,j) = (0 : nfft/2 - 1) * Fs / nfft; 
        %subplot(3,9,2*j);
        %plot(f(:,j),Mag);
        %title('power spectrum');xlabel('f(Hz)');ylabel('I');
        pos1 = (nfft - 1)/20 + 1;
        %I1(j,i) = Mag(pos1);
%         basefreq(j,i) = f(pos1);
        pos3 = 3*(nfft - 1)/20 + 1;
        i31(j) = Mag(pos3) / Mag(pos1);
%         thirdfreq(j,i) = f(pos3);
        %I31(j,i) = I3(j,i)/I1(j,i);
    end;
    
    
%% For high frequency, use ode45 to calculate
    for j = 15 : loop_j
        Freq = 10^(-3 + 0.2*(j - 1));
       % freq(j) = Freq;
        st(j) = Rate / Freq;
        %strain(j,i) = Rate / Freq;
        A0 = [1;1;1;0;0;0];
        tspan = [0 : pi / (10*Freq) : 80 * j * pi / Freq];%time range for ODE solver
        [t,A] = ode45(@(t, A) Lodefun_par(t, A, Rate, Freq),tspan,A0);
        sigma = G * A(:,4);
        Fs = (10*Freq)/pi;   %sampling frequency
        ts = tspan(800*j-199 : 800*j+1);    %time vector
        sigmas = sigma(800*j-199 : 800*j+1);    %select FFT signal
        sigma_max = max(sigmas);
        sigma_min = min(sigmas);
        sta(j) = (sigma_max - sigma_min) / 2;
        %stressAmplitude(j,i) = (sigma_max - sigma_min) / 2;
        %subplot(3,9,2*j-1);
        %title('stress signal');xlabel('t');ylabel('stress(Pa)');
        %plot(ts,sigmas);
%% FFT
        nfft = length(sigmas);
        Fsigma = fft(sigmas,nfft);
        Fsigma = Fsigma(1:nfft/2);
        Mag = abs(Fsigma);
       % f(:,j) = (0 : nfft/2 - 1) * Fs / nfft; 
        %subplot(3,9,2*j);
        %plot(f(:,j),Mag);
        %title('power spectrum');xlabel('f(Hz)');ylabel('I');
        pos1 = (nfft - 1)/20 + 1;
        %I1(j,i) = Mag(pos1);
%         basefreq(j,i) = f(pos1);
        pos3 = 3*(nfft - 1)/20 + 1;
        %I3(j,i) = Mag(pos3);
%         thirdfreq(j,i) = f(pos3);
        %I31(j,i) = I3(j,i)/I1(j,i);
        i31(j) = Mag(pos3) / Mag(pos1);
    end;
    
    strain(:,i)= st;
    stressAmplitude(:,i)= sta;
    I31(:,i) = i31;
    t_run = toc;
    run_time(i) = t_run;
end;
%% store the data calculated 
csvwrite('LAOS_strain_2core.csv', strain);
csvwrite('LAOS_stressAmplitude_2core.csv', stressAmplitude);
csvwrite('LAOS_I31_2core.csv', I31);
csvwrite('LAOS_runtime_2core.csv', run_time);