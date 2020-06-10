%% Term 3

C_a2_3 = 200/37; % Constant 3
K_a2_3 = Wn^2/C_a2_3; % Controller Gain 3

% Plant TF, 'Gp3'
Gp3 = zpk(minreal(C_a2_3/s^2))

% Chosen Youla Parameter, 'Y3' -> Y3(0) = 0
Y3 = zpk(minreal(((K_a2_3*s^2)*(tz*s + 1)/((s^2 + 2*Z*Wn*s + Wn^2)*(tp*s + 1))),1e-05))

% Complementary Sensitivity TF, 'T3' -> T3(0) = 1 (1st interpolation
% condition)
T3 = zpk(minreal((Y3*Gp3),1e-05))

% Sensitivity TF, 'S3'
S3 = zpk(minreal((1-T3),1e-05))

% Controller TF, 'Gc3'
Gc3 = zpk(minreal((Y3/S3),1e-05))

% Return Ratio, 'L3'
L3 = zpk(minreal((Gc3*Gp3),1e-05))

GpS3 = zpk(minreal((Gp3*S3),1e-05))

% Internal stability check
Y3_stability = isstable(Y3)
T3_stability = isstable(T3)
S3_stability = isstable(S3)
GpS3_stability = isstable(GpS3)

M2_3 = 1/getPeakGain(S3) % M2-margin
BW_3 = bandwidth(T3) % Bandwidth of the closed-loop
AE_3 = getPeakGain(Y3) % Maximum actuator effort

figure(1)
bodemag(Y3, S3, T3);
legend('Y3','S3','T3');