% Astrobee Model

%% System

% A Matrix

load Matrices/A_matrix.mat
A = A_matrix

% B Matrix: Stowed

load Matrices/B_stowed.mat
B = B_stowed

% Full-State Feedback

Cf = eye(12);
Df = [zeros(12, 6)];

sys_full = ss(A, B, Cf, Df);

tf_full = tf(sys_full);

syms s

tf_full_sym = simplify(Cf * inv(s * eye(12) - A) * B + Df);
pretty(tf_full_sym)

%% Section 1: Attitude Controller Design -> Marginally Stable Pole at the Origin

% Top-half matrix (a1): 

% Satisfaction of the first interpolation condition

attitude_full = [tf_full_sym(4:6, 4:6); tf_full_sym(10:12, 4:6)];
pretty(attitude_full)

Gp_a1 = attitude_full(1:3, 1:3);
P_a1 = Gp_a1 * s;

[UL_a1, UR_a1, S_a1] = smithForm(P_a1, s)

Mp_a1 = S_a1/s

K_a1 = 1;
tp_a1 = 100;

Y1 = (K_a1 * s)/(tp_a1 * s + 1);
Y2 = Y1;
Y3 = Y1;

My_a1 = diag([Y1 Y2 Y3])

Mt = Mp_a1 * My_a1

Gc_a1_sym = simplify((UR_a1 * inv(eye(size(My_a1 * Mp_a1)) - My_a1 * Mp_a1) * My_a1 * UL_a1))

% % Convert to a string
% Gc_a1_arr = [];
% Gc_a1_str = [];
% for i = 1:size(Gc_a1_sym, 1)
%         Gc_a1_str = char(Gc_a1_sym(i, i));
%         % Define ?s? as transfer function variable 
%         s = tf('s');
%         % Evaluate the expression:
%         eval(['Gc_a1_arr(i) = ', Gc_a1_str])
% end
% 
% Gc_a1 = diag(Gc_a1_arr)

%% Section 2: Translation Controller Design -> Unstable Double-Pole at the Origin

% Bottom-half matrix (a2):

% Run this section first to calculate 'tz' to ensure that the second interpolation condition is satisfied

% d^k(T)/ds^k|(s=0) = 0, where k = 1 (since there is a double unstable pole
% (multiplicity ap = 2) in the plant at s = 0; k = ap - 1) -> 2nd
% interpolation condition

Wn = 0.01; % Natural Frequency of the Control System
Z = 2^-0.5; % Damping Ratio
tp = 1/(10*Wn); % Time Constant (of the included pole) 

C_a2_1 = 200/37; % Constant 1
K_a2_1 = Wn^2/C_a2_1; % Controller Gain 1

C_a2_2 = 500/101; % Constant 2
K_a2_2 = Wn^2/C_a2_2; % Controller Gain 2

C_a2_3 = 250/47; % Constant 3
K_a2_3 = Wn^2/C_a2_3; % Controller Gain 3

syms s tz1 tz2 tz3

% tz1

TF1 = ((K_a2_1*C_a2_1)*(tz1*s + 1))/((s^2 + 2*Z*Wn*s + Wn^2)*(tp*s + 1))
dTF1 = diff(TF1,s)
eqn1 = subs(dTF1,s,0) == 0;
tz1 = solve(eqn1,tz1)

% tz2

TF2 = ((K_a2_2*C_a2_2)*(tz2*s + 1))/((s^2 + 2*Z*Wn*s + Wn^2)*(tp*s + 1))
dTF2 = diff(TF2,s)
eqn2 = subs(dTF2,s,0) == 0;
tz2 = solve(eqn2,tz2)

% tz3

TF3 = ((K_a2_3*C_a2_3)*(tz3*s + 1))/((s^2 + 2*Z*Wn*s + Wn^2)*(tp*s + 1))
dTF3 = diff(TF3,s)
eqn3 = subs(dTF3,s,0) == 0;
tz3 = solve(eqn3,tz3)

%% Section 3: Translation Controller Design -> Unstable Double-Pole at the Origin

% Youla Control Design

s = tf('s');

% Constants & Design Parameters
Wn = 0.01; % Natural Frequency of the Control System
Z = 2^-0.5; % Damping Ratio
tp = 1/(10*Wn); % Time Constant of the added pole 
tz = 100*2^(1/2) + 10;

%% Term 1

C_a2_1 = 200/37; % Constant 1
K_a2_1 = Wn^2/C_a2_1; % Controller Gain 1

% Plant TF, 'Gp1'
Gp1 = zpk(minreal(C_a2_1/s^2))

% Chosen Youla Parameter, 'Y1' -> Y1(0) = 0
Y1 = zpk(minreal(((K_a2_1*s^2)*(tz*s + 1)/((s^2 + 2*Z*Wn*s + Wn^2)*(tp*s + 1))),1e-05))

% Complementary Sensitivity TF, 'T1' -> T1(0) = 1 (1st interpolation
% condition)
T1 = zpk(minreal((Y1*Gp1),1e-05))

% Sensitivity TF, 'S1'
S1 = zpk(minreal((1-T1),1e-05))

% Controller TF, 'Gc1'
Gc1 = zpk(minreal((Y1/S1),1e-05))

% Return Ratio, 'L1'
L1 = zpk(minreal((Gc1*Gp1),1e-05))

GpS1 = zpk(minreal((Gp1*S1),1e-05))

% Internal stability check
Y1_stability = isstable(Y1)
T1_stability = isstable(T1)
S1_stability = isstable(S1)
GpS1_stability = isstable(GpS1)

M2_1 = 1/getPeakGain(S1) % M2-margin
BW_1 = bandwidth(T1) % Bandwidth of the closed-loop
AE_1 = getPeakGain(Y1) % Maximum actuator effort

figure(1)
bodemag(Y1, S1, T1);
legend('Y1','S1','T1');


%% Term 2

C_a2_2 = 500/101; % Constant 2
K_a2_2 = Wn^2/C_a2_2; % Controller Gain 2

% Plant TF, 'Gp2'
Gp2 = zpk(minreal(C_a2_2/s^2))

% Chosen Youla Parameter, 'Y2' -> Y2(0) = 0
Y2 = zpk(minreal(((K_a2_2*s^2)*(tz*s + 1)/((s^2 + 2*Z*Wn*s + Wn^2)*(tp*s + 1))),1e-05))

% Complementary Sensitivity TF, 'T2' -> T2(0) = 1 (1st interpolation
% condition)
T2 = zpk(minreal((Y2*Gp2),1e-05))

% Sensitivity TF, 'S2'
S2 = zpk(minreal((1-T2),1e-05))

% Controller TF, 'Gc2'
Gc2 = zpk(minreal((Y2/S2),1e-05))

% Return Ratio, 'L2'
L2 = zpk(minreal((Gc2*Gp2),1e-05))

GpS2 = zpk(minreal((Gp2*S2),1e-05))

% Internal stability check
Y2_stability = isstable(Y2)
T2_stability = isstable(T2)
S2_stability = isstable(S2)
GpS2_stability = isstable(GpS2)

M2_2 = 1/getPeakGain(S2) % M2-margin
BW_2 = bandwidth(T2) % Bandwidth of the closed-loop
AE_2 = getPeakGain(Y2) % Maximum actuator effort

figure(1)
bodemag(Y2, S2, T2);
legend('Y2','S2','T2');

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

%%

Gc_a2 = [tf(Gc1) 0 0; 0 tf(Gc2) 0; 0 0 tf(Gc3)]

% Convert to symbolic matrix:
syms s

% Term 1
[Num1,Den1] = tfdata(tf(Gc1), 'v')
Gc_a2_sym_term_1 = poly2sym(Num1, s)/poly2sym(Den1, s)

% Term 2
[Num2,Den2] = tfdata(tf(Gc2), 'v')
Gc_a2_sym_term_2 = poly2sym(Num2, s)/poly2sym(Den2, s)

% Term 3
[Num3,Den3] = tfdata(tf(Gc3), 'v')
Gc_a2_sym_term_3 = poly2sym(Num3, s)/poly2sym(Den3, s)

Gc_a2_sym = diag([Gc_a2_sym_term_1 Gc_a2_sym_term_2 Gc_a2_sym_term_3])

Gc_a = [Gc_a1_sym Gc_a2_sym]

%% Coordinate Feedback

% Cc = [zeros(6, 12)];
% Cc(1:6, 1:6) = eye(6);
% 
% Dc = [zeros(6, 6)];
% 
% sys_coord = ss(A, B, Cc, Dc);
% 
% tf_coord = tf(sys_coord);
% 
% syms s
% 
% tf_coord_sym = simplify(Cc * inv(s * eye(12) - A) * B + Dc);
% pretty(tf_coord_sym)
% 
% translation_coord = [tf_coord_sym(1:3, 1:3); tf_coord_sym(7:9, 1:3)];
% pretty(translation_coord)
% 
% attitude_coord = [tf_coord_sym(4:6, 4:6); tf_coord_sym(10:12, 4:6)];
% pretty(attitude_coord)



