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

%% Translation Controller Design

translation_full = [tf_full_sym(1:3, 1:3); tf_full_sym(7:9, 1:3)];
pretty(translation_full);

% translation 1 (top half matrix)

Gp_t1 = translation_full(1:3, 1:3);
P_t1 = Gp_t1 * s;

[UL_t1, UR_t1, S_t1] = smithForm(P_t1, s)

Mp_t1 = S_t1/s

K_t1 = 1;
tp_t1 = 100;

Y1 = (K_t1 * s)/(tp_t1 * s + 1);
Y2 = Y1;
Y3 = Y1;

My_t1 = diag([Y1 Y2 Y3])

Mt = Mp_t1 * My_t1

Gc_t1_sym = simplify((UR_t1 * inv(eye(size(My_t1 * Mp_t1)) - My_t1 * Mp_t1) * My_t1 * UL_t1))

% % Convert to a string
% Gc_t1_arr = [];
% Gc_t1_str = [];
% for i = 1:size(Gc_t1_sym, 1)
%         Gc_t1_str = char(Gc_t1_sym(i, i));
%         % Define ?s? as transfer function variable 
%         s = tf('s');
%         % Evaluate the expression:
%         eval(['Gc_t1_arr(i) = ', Gc_t1_str])
% end
% 
% Gc_t1 = diag(Gc_t1_arr)

%% Run this section first to calculate 'tz' to ensure that the second interpolation condition is satisfied

% d^k(T)/ds^k|(s=0) = 0, where k = 1 (since there is a double unstable pole
% (multiplicity ap = 2) in the plant at s = 0; k = ap - 1) -> 2nd
% interpolation condition

C_t2 = 7939/500; % Constant
Wn = 0.01; % Natural Frequency of the Control System
K = Wn^2/C_t2; % Controller Gain
Z = 2^-0.5; % Damping Ratio
tp = 1/(10*Wn); % Time Constant of the added pole 

syms s tz

TF = ((K*C_t2)*(tz*s + 1))/((s^2 + 2*Z*Wn*s + Wn^2)*(tp*s + 1))
dTF = diff(TF,s)
eqn = subs(dTF,s,0) == 0;
tz = solve(eqn,tz)

%% Youla Control Design (Translation)

s = tf('s');

% Constants & Design Parameters
C_t2 = 500/7939; % Constant
Wn = 0.01; % Natural Frequency of the Control System
K = Wn^2/C_t2; % Controller Gain
Z = 2^-0.5; % Damping Ratio
tp = 1/(10*Wn); % Time Constant of the added pole 
tz = 100*2^(1/2) + 10;

% Plant TF, 'Gp'
Gp = zpk(minreal(C_t2/s^2))

% Chosen Youla Parameter, 'Y' -> Y(0) = 0
Y = zpk(minreal(((K*s^2)*(tz*s + 1)/((s^2 + 2*Z*Wn*s + Wn^2)*(tp*s + 1))),1e-05))

% Complementary Sensitivity TF, 'T' -> T(0) = 1 (1st interpolation
% condition)
T = zpk(minreal((Y*Gp),1e-05))

% Sensitivity TF, 'S'
S = zpk(minreal((1-T),1e-05))

% Controller TF, 'Gc'
Gc = zpk(minreal((Y/S),1e-05))

% Return Ratio, 'L'
L = zpk(minreal((Gc*Gp),1e-05))

GpS = zpk(minreal((Gp*S),1e-05))

% Internal stability check
Y_stability = isstable(Y)
T_stability = isstable(T)
S_stability = isstable(S)
GpS_stability = isstable(GpS)

M2 = 1/getPeakGain(S) % M2-margin
BW = bandwidth(T) % Bandwidth of the closed-loop
AE = getPeakGain(Y) % Maximum actuator effort

figure(1)
bodemag(Y, S, T);
legend('Y','S','T');

Gc_t2 = [tf(Gc) 0 0; 0 tf(Gc) 0; 0 0 tf(Gc)]

% Convert to symbolic matrix
[Num,Den] = tfdata(tf(Gc), 'v')
syms s
Gc_t2_sym_term = poly2sym(Num, s)/poly2sym(Den, s)
Gc_t2_sym = diag([Gc_t2_sym_term Gc_t2_sym_term Gc_t2_sym_term])

Gc_t = [Gc_t1_sym Gc_t2_sym]

%% Attitude Controller Design

attitude_full = [tf_full_sym(4:6, 4:6); tf_full_sym(10:12, 4:6)];
pretty(attitude_full)

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
% attitude_full = [tf_coord_sym(4:6, 4:6); tf_coord_sym(10:12, 4:6)];
% pretty(attitude_coord)



