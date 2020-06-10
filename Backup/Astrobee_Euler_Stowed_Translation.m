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

%% Section 1: Translation Controller Design -> Marginally Stable Pole at the Origin

% Top-half matrix (t1): 

% Satisfaction of the first interpolation condition

translation_full = [tf_full(1:3, 1:3); tf_full(7:9, 1:3)];
% pretty(translation_full);

G_t1 = translation_full(1:3, 1:3);

Gp_t1 = translation_full(1, 1);

K_t1 = 500/7939;

% P_t1 = Gp_t1 * s;
% 
% [UL_t1, UR_t1, S_t1] = smithForm(P_t1, s)
% 
% Mp_t1 = S_t1/s
% 
% K_t1 = 1;
% tp_t1 = 100;
% 
% Y1 = (K_t1 * s)/(tp_t1 * s + 1);
% Y2 = Y1;
% Y3 = Y1;
% 
% My_t1 = diag([Y1 Y2 Y3])
% 
% Mt = Mp_t1 * My_t1
% 
% Gc_t1_sym = simplify((UR_t1 * inv(eye(size(My_t1 * Mp_t1)) - My_t1 * Mp_t1) * My_t1 * UL_t1))
% 
% Gc_t1 = tf(double(Gc_t1_sym));

s = tf('s');

tp_t1 = 100;

% Chosen Youla Parameter, 'Y' -> Y(0) = 0
Y_t1 = s/(K_t1*(tp_t1 * s + 1));

% Complementary Sensitivity TF, 'T' -> T(0) = 1 (1st interpolation
% condition)
T_t1 = minreal((Y_t1*Gp_t1),1e-04)

% Sensitivity TF, 'S'
S_t1 = minreal((1-T_t1),1e-04)

% Controller TF, 'Gc'
Gc_t1 = minreal((Y_t1/S_t1),1e-04)

% Return Ratio, 'L'
L_t1 = minreal((Gc_t1*Gp_t1),1e-04)

GpS_t1 = minreal((Gp_t1*S_t1),1e-04)

% Internal stability check
Y_t1_stability = isstable(Y_t1)
T_t1_stability = isstable(T_t1)
S_t1_stability = isstable(S_t1)
GpS_t1_stability = isstable(GpS_t1)

M2_t2 = 1/getPeakGain(S_t1) % M2-margin
BW_t2 = bandwidth(T_t1) % Bandwidth of the closed-loop
AE_t2 = getPeakGain(Y_t1) % Maximum actuator effort

figure(1)
bodemag(Y_t1, S_t1, T_t1);
legend('Y_t1','S_t1','T_t1');

Gc_1 = Gc_t1 * eye(3);


%% Section 2: Translation Controller Design -> Unstable Double-Pole at the Origin

% Bottom-half matrix (t2):

% Run this section first to calculate 'tz' to ensure that the second interpolation condition is satisfied

% d^k(T)/ds^k|(s=0) = 0, where k = 1 (since there is a double unstable pole
% (multiplicity ap = 2) in the plant at s = 0; k = ap - 1) -> 2nd
% interpolation condition

C_t2 = 500/7939; % Constant
Wn = 3.25; % Natural Frequency of the Control System
K_t2 = Wn^2/C_t2; % Controller Gain
Z = 2^-0.5; % Damping Ratio
tp_t2 = 1/(10*Wn); % Time constant (of the included pole)

syms s tz

TF = ((K_t2*C_t2)*(tz*s + 1))/((s^2 + 2*Z*Wn*s + Wn^2)*(tp_t2*s + 1))
dTF = diff(TF,s)
eqn = subs(dTF,s,0) == 0;
tz = double(solve(eqn,tz))

%% Section 3: Translation Controller Design -> Unstable Double-Pole at the Origin

% Youla Control Design

s = tf('s');

% % Constants & Design Parameters
% C_t2 = 500/7939; % Constant
% Wn = 3.25; % Natural Frequency of the Control System
% K = Wn^2/C_t2; % Controller Gain
% Z = 2^-0.5; % Damping Ratio
% tp = 1/(10*Wn); % Time Constant of the added pole 
% tz = (4*2^(1/2))/13 + 2/65; % 100*2^(1/2) + 10;

% Plant TF, 'Gp'
Gp_t2 = zpk(minreal(C_t2/s^2))

% Chosen Youla Parameter, 'Y' -> Y(0) = 0
Y_t2 = zpk(minreal(((K_t2*s^2)*(tz*s + 1)/((s^2 + 2*Z*Wn*s + Wn^2)*(tp_t2*s + 1))),1e-04))

% Complementary Sensitivity TF, 'T' -> T(0) = 1 (1st interpolation
% condition)
T_t2 = zpk(minreal((Y_t2*Gp_t2),1e-04))

% Sensitivity TF, 'S'
S_t2 = zpk(minreal((1-T_t2),1e-04))

% Controller TF, 'Gc'
Gc_t2 = zpk(minreal((Y_t2/S_t2),1e-04))

% Return Ratio, 'L'
L_t2 = zpk(minreal((Gc_t2*Gp_t2),1e-04))

GpS_t2 = zpk(minreal((Gp_t2*S_t2),1e-04))

% Internal stability check
Y_t2_stability = isstable(Y_t2)
T_t2_stability = isstable(T_t2)
S_t2_stability = isstable(S_t2)
GpS_t2_stability = isstable(GpS_t2)

M2_t2 = 1/getPeakGain(S_t2) % M2-margin
BW_t2 = bandwidth(T_t2) % Bandwidth of the closed-loop
AE_t2 = getPeakGain(Y_t2) % Maximum actuator effort

figure(1)
bodemag(Y_t2, S_t2, T_t2);
legend('Y_t2','S_t2','T_t2');

Gc_2 = Gc_t2 * eye(3);

%% Simulation

Gp = minreal([tf_full(1:3, 1:3); tf_full(7:9, 1:3)], 1e-04);
Gc = [Gc_1 Gc_2]
Lu = minreal(Gc * Gp, 1e-04);
Ly = minreal(Gp * Gc, 1e-04);
Y = minreal(inv(eye(3) + Lu) * Gc); 
Ty = minreal(inv(eye(6) + Ly) * Ly); 
Sy = minreal(inv(eye(6) + Ly), 1e-04);
Su = minreal(inv(eye(3) + Lu), 1e-04);

figure
step(Ty);

figure
step(Y);

figure
sigma(Y, Ty, Sy, Su)
[l, hObj] = legend('$Y$', '$T_{y}$', '$S_{y}$', '$S_{u}$','Interpreter','latex','FontSize', 12);
set(l,'string',{'$Y$', '$T_{y}$', '$S_{y}$', '$S_{u}$'});
hL = findobj(hObj,'type','line');
set(hL,'linewidth', 2); 

figure
sigma(Gc, Gp, Ly, Y)
[l, hObj] = legend('$G_{c}$', '$G_{p}$', '$L_{y}$', '$Y$','Interpreter','latex','FontSize', 12);
set(l,'string',{'$G_{c}$', '$G_{p}$', '$L_{y}$', '$Y$'});
hL = findobj(hObj,'type','line');
set(hL,'linewidth', 2);

figure
sigma(Gc, Gp, Y)
[l, hObj] = legend('$G_{c}$', '$G_{p}$', '$Y$','Interpreter','latex','FontSize', 12);
set(l,'string',{'$G_{c}$', '$G_{p}$', '$Y$'});
hL = findobj(hObj,'type','line');
set(hL,'linewidth', 2);

figure
sigma(Ly, Sy, Ty)
[l, hObj] = legend('$L_{y}$', '$S_{y}$', '$T_{y}$','Interpreter','latex','FontSize', 12);
set(l,'string',{'$L_{y}$', '$S_{y}$', '$T_{y}$'});
hL = findobj(hObj,'type','line');
set(hL,'linewidth', 2);

figure
sigma(Sy, Su)
[l, hObj] = legend('$S_{y}$', '$S_{u}$','Interpreter','latex','FontSize', 12);
set(l,'string',{'$S_{y}$', '$S_{u}$'});
hL = findobj(hObj,'type','line');
set(hL,'linewidth', 2); 

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



