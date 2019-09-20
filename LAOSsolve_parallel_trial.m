%
clear all;
% G is the modulus of the polymer, the unit is Pa. This is an adjustable parameter.
G = 320.33;
global Rate
global Freq
% create a series of shear rate ranging from 10^-3 ~ 10^3 s-1

% for low frequency, use ode15 to calculate
for i = 1 : 36
    Rate = 10^(-3 + 0.2*(i - 1));
    rate(i) = Rate;
    %figure(i);
    for j = 1 : 14
        Freq = 10^(-3 + 0.2*(j - 1)); % frequency from 10^-3 ~ 10^-0.4
        freq(j) = Freq;
        strain(j,i) = rate(i) / freq(j);
        A0 = [1;1;1;0;0;0]; %initialize the configuration to isotropic state
        tspan = [0 : pi/(10*freq(j)) : 80 * pi / freq(j)]; %time range for ODE solver
        [t,A] = ode15s('Lodefun',tspan,A0);
        sigma = G * A(:,4);
        % parameters for FFT
        Fs = (10*freq(j))/pi;   %sampling frequency
        %When shear rate is low and time is small(tspan(1:600)), the
        %system is in SAOS region, the constitutive equation does not work.
        %The result is not accurate.
        ts = tspan(601 : 801);  %time vector
        sigmas = sigma(601 : 801);  %select FFT signal
        sigma_max = max(sigmas);
        sigma_min = min(sigmas);
        stressAmplitude(j,i) = (sigma_max - sigma_min) / 2; % amplitude
        %subplot(3,9,2*j-1);
        %plot(ts,sigmas);
        %title('stress signal');xlabel('t');ylabel('stress(Pa)');
        nfft = length(sigmas);
        Fsigma = fft(sigmas,nfft);
        Fsigma = Fsigma(1:nfft/2);
        Mag = abs(Fsigma);
        f(:,j) = (0 : nfft/2 - 1) * Fs / nfft; 
        %subplot(3,9,2*j);
        %plot(f(:,j),Mag);
        %title('power spectrum');xlabel('f(Hz)');ylabel('I');
        pos1 = (nfft - 1)/20 + 1;
        I1(j,i) = Mag(pos1);
        basefreq(j,i) = f(pos1);
        pos3 = 3*(nfft - 1)/20 + 1;
        I3(j,i) = Mag(pos3);
        thirdfreq(j,i) = f(pos3);
        I31(j,i) = I3(j,i)/I1(j,i);
    end;
    % for high frequency, use ode15 to calculate
        for j = 15 : 31
        Freq = 10^(-3 + 0.2*(j - 1));
        freq(j) = Freq;
        strain(j,i) = rate(i) / freq(j);
        A0 = [1;1;1;0;0;0];
        tspan = [0 : pi / (10*freq(j)) : 80 * j * pi / freq(j)];%time range for ODE solver
        [t,A] = ode45('Lodefun',tspan,A0);
        sigma = G * A(:,4);
        Fs = (10*freq(j))/pi;   %sampling frequency
        ts = tspan(800*j-199 : 800*j+1);    %time vector
        sigmas = sigma(800*j-199 : 800*j+1);    %select FFT signal
        sigma_max = max(sigmas);
        sigma_min = min(sigmas);
        stressAmplitude(j,i) = (sigma_max - sigma_min) / 2;
        %subplot(3,9,2*j-1);
        %title('stress signal');xlabel('t');ylabel('stress(Pa)');
        %plot(ts,sigmas);
        nfft = length(sigmas);
        Fsigma = fft(sigmas,nfft);
        Fsigma = Fsigma(1:nfft/2);
        Mag = abs(Fsigma);
        f(:,j) = (0 : nfft/2 - 1) * Fs / nfft; 
        %subplot(3,9,2*j);
        %plot(f(:,j),Mag);
        %title('power spectrum');xlabel('f(Hz)');ylabel('I');
        pos1 = (nfft - 1)/20 + 1;
        I1(j,i) = Mag(pos1);
        basefreq(j,i) = f(pos1);
        pos3 = 3*(nfft - 1)/20 + 1;
        I3(j,i) = Mag(pos3);
        thirdfreq(j,i) = f(pos3);
        I31(j,i) = I3(j,i)/I1(j,i);
    end;
%     filename = 'LAOS_result.xlsx';
%     xlswrite(filename,strain,1);
%     xlswrite(filename,stressAmplitude,2);
%     xlswrite(filename,I31,3);
end;

%store the data calculated 
csvwrite('LAOS_strain.csv', strain);
csvwrite('LAOS_stressAmplitude.csv', stressAmplitude);
csvwrite('LAOS_I31.csv', I31);
