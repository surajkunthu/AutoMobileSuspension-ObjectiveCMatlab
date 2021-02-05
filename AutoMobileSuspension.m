% Suraj Kunthu           %%%
% Auto-mobile Suspension %%%
% Mass Spring Damper     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To test a part of the project uncheck comments '%' for specified sections

% Clear all ---------------------------------------------------------------
clear all;
clc;
close all;

% Given Data --------------------------------------------------------------
mb = 650;                             % [lbf]
ms = 60;                              % [lbf]
kt = 1500;                            % [lbf/in]
ct = 4;                               % [lbf-s/in]

% Unit Conversion ---------------------------------------------------------
mb = mb/(32.2*12);                    % [lbf] -> [blobs]
ms = ms/(32.2*12);                    % [lbf] -> [blobs]

% Unknown Criteria --------------------------------------------------------
% ks                                  % [lbf/in]
% cs                                  % [lbf-s/in]

%% Part 2 -----------------------------------------------------------------
% Specify the suspension spring constant, such that the undamped natural 
% frequency of the automobile body mass is approximately 1 Hz.
% -------------------------------------------------------------------------

% Set cs = 0 so that there is no damping acting on the mass ms
% When ks = 69.5, undamped nat. freq. wn1 = 1 [Hz] OR 2*pi (6.28318)[rad/s]
%cs = 0;                              % [lbf-s/in]
%ks = 69.5;                           % [lbf/in]

% Transfer Function -------------------------------------------------------
%s = tf('s');
%A = [(mb*s^2)+(cs*s)+ks, -(cs*s+ks); -(cs*s+ks), (ms*s^2)+(cs*s)+(ks+kt)];
%B = [0; kt];
%X_F = A\B;

% Function that returns natural frequencies, damping ratios, and poles ----
%[wn,zeta,P] = damp(X_F(1))

% Displacement due to Step Response ---------------------------------------
%step_size = 1;                        % [in]
%figure
%step(step_size * X_F(1))
%hold on
%step(step_size * X_F(2))
%title('Displacement Due to Road Step')
%ylabel('Displacement (in)')
%legend('Body Mass','Suspension Mass')

%% Part 3 -----------------------------------------------------------------
% Investigate the effects of the suspension damping constant, cs, on both 
% the maximum force transmitted to the automobile body mass and the 
% settling time of the automobile body mass in response to hitting a 1? pot 
% hole modeled as a the step function z(t) = u_s(t).
% -------------------------------------------------------------------------
% Keeping ks constant at 69.5 [rad/s]

%cs = 10000;                           % [lbf-s/in]
% cs values tested: 0, 1, 5, 10, 50, 100, 500, 1000, 5000, 10000
%ks = 69.5;                            % [lbf/in]

% State-Space Model -------------------------------------------------------
%A = [0, 0, 1, 0; 0, 0, 0, 1; -ks/mb, ks/mb, -cs/mb, cs/mb;...
%    ks/ms, -(ks+kt)/ms, cs/ms, -cs/ms];
%B = [0; 0; 0; kt/ms];
%C = [1, 0, 0, 0; 0, 1, 0, 0; -ks, ks, -cs, cs; ks, -(ks+kt), cs, -cs];
%D = [0; 0; 0; kt];
%sys = ss(A, B, C, D);
%Roots = eig(A)

% Step Function -----------------------------------------------------------
%stepsize = 1;                         % [in]
%[y,t] = step(stepsize*sys);
%figure
%plot(t,y(:,3:4))
%legend('f_b', 'f_s')
%xlabel('Time(s)')
%ylabel('Force (lbf)')
%title('Force Transmission in Quarter Car Suspension')

%% Part 4 -----------------------------------------------------------------
% Specify the suspension damping constant, cs, such that the maximum force
% and settling time in response to the 1? pot hole are approximately 525 
% lbf and 2 s, respectively.
% -------------------------------------------------------------------------
% Keeping ks constant at 69.5 [rad/s]

%cs = 7.17726;                         % [lbf-s/in]
%ks = 69.5;                            % [lbf/in]

% State-Space Model -------------------------------------------------------
%A = [0, 0, 1, 0; 0, 0, 0, 1; -ks/mb, ks/mb, -cs/mb, cs/mb;...
%    ks/ms, -(ks+kt)/ms, cs/ms, -cs/ms];
%B = [0; 0; 0; kt/ms];
%C = [1, 0, 0, 0; 0, 1, 0, 0; -ks, ks, -cs, cs; ks, -(ks+kt), cs, -cs];
%D = [0; 0; 0; kt];
%sys = ss(A, B, C, D);
%Roots = eig(A)

% Step Function -----------------------------------------------------------
%stepsize = 1;                         % [in]
%[y,t] = step(stepsize*sys);
%figure
%plot(t,y(:,3))
%legend('f_b')
%xlabel('Time(s)')
%ylabel('Force (lbf)')
%title('Force Transmission in Quarter Car Suspension')
%MaxForce = max(y);
%disp('Maximum Force');
%disp(MaxForce(3))

%% Part 5 -----------------------------------------------------------------
% Plot the frequency response of the force transmissibility transfer 
% function. Note the resonant frequency and amplification factor 
% at resonance.
% -------------------------------------------------------------------------
%cs = 7.17726;                         % [lbf-s/in]
%ks = 69.5;                            % [lbf/in]

% State-Space Model -------------------------------------------------------
%A = [0, 0, 1, 0; 0, 0, 0, 1; -ks/mb, ks/mb, -cs/mb, cs/mb;...
%    ks/ms, -(ks+kt)/ms, cs/ms, -cs/ms];
%B = [0; 0; 0; kt/ms];
%C = [1, 0, 0, 0; 0, 1, 0, 0; -ks, ks, -cs, cs; ks, -(ks+kt), cs, -cs];
%D = [0; 0; 0; kt];
%sys = ss(A, B, C, D);
%Roots = eig(A);

% Bode Plot Frequency Response of the Force Transmissibility TF -----------
%[mag, phase, wout] = bode(sys);
%f = wout/2/pi;                         % [Hz]
%M = squeeze(mag(1,1,:));
%figure
%plot(f, M)
%xlim([0,50]);
%title('Frequency Response of the Force Transmissibility Transfer Function')
%xlabel('Frequency [Hz]')
%ylabel('M  - Amplitude Ratio')
%Mmax = max(M);
%freq = f;
%disp('Max. Amplitude Ratio')
%disp(Mmax)
%disp('Frequency at Max. Amplitude Ratio [Hz]')
%disp(freq(24))

% Resonant Frequency Calculation ------------------------------------------
%w = 2*pi*freq(24);                      % [rad/s]
%wn = 2*pi;                              % [rad/s]
%r = w/wn;                               
%zeta = sqrt((-((r)^2)+1)/(2));
%wr = wn*sqrt(1-2*(zeta^2));             % [rad/s]
%disp('       w        wn         r       zeta       wr')
%disp([w, wn, r, zeta, wr])
%disp('Amplification Factor')
%disp([0, r*wn])

