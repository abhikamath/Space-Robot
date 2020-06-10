% Astrobee Model

%% System

% A Matrix

load Matrices/A_matrix.mat
A = A_matrix

% B Matrix: Stowed

load Matrices/B_stowed.mat
B = B_stowed

% % Full-State Feedback
% 
% Cf = eye(12);
% 
% Df = [zeros(12, 6)];
% 
% sys_full = ss(A, B, Cf, Df);
% 
% tf_full = tf(sys_full);
% 
% syms s
% 
% tf_full_sym = simplify(Cf * inv(s * eye(12) - A) * B + Df);
% pretty(tf_full_sym)

% Coordinate Feedback

Cc = [zeros(6, 12)];
Cc(1:6, 7:12) = eye(6);

Dc = [zeros(6, 6)];

sys_coord = ss(A, B, Cc, Dc);

tf_coord = tf(sys_coord);

syms s

tf_coord_sym = simplify(Cc * inv(s * eye(12) - A) * B + Dc);
pretty(tf_coord_sym)

% translation_coord_sym = [tf_coord_sym(4:6, 7:9); tf_coord_sym(7:9, 7:9)];
% pretty(translation_coord_sym)

% attitude_coord_sym = [tf_coord_sym(4:6, 4:6); tf_coord_sym(10:12, 4:6)];
% pretty(attitude_coord_sym)

%% Section 1: Translation Controller Design -> Unstable Double-Pole at the Origin

% Run this section first to calculate 'tz' to ensure that the second interpolation condition is satisfied

% d^k(T)/ds^k|(s=0) = 0, where k = 1 (since there is a double unstable pole
% (multiplicity ap = 2) in the plant at s = 0; k = ap - 1) -> 2nd
% interpolation condition

C = 250/47; % Constant
Wn = 1; % Natural Frequency of the Control System
K = Wn^2/C; % Controller Gain
Z = 2^-0.5; % Damping Ratio
tp = 1/(10*Wn); % Time constant (of the included pole)
tpx = 0.5; % Time constant (of the pole included to drop Youla at high frequencies)

syms s tz

T_eqn = ((K*C)*(tz*s + 1)/((s^2 + 2*Z*Wn*s + Wn^2)*(tp*s + 1)*(tpx*s + 1)^2));
dT_eqn = diff(T_eqn,s);
eqn = subs(dT_eqn,s,0) == 0;
tz = double(solve(eqn,tz))

%% Section 3: Translation Controller Design -> Unstable Double-Pole at the Origin

% Youla Control Design

s = tf('s');

% % Constants & Design Parameters
% C = 500/7939; % Constant
% Wn = 3.25; % Natural Frequency of the Control System
% K = Wn^2/C; % Controller Gain
% Z = 2^-0.5; % Damping Ratio
% tp = 1/(10*Wn); % Time Constant of the added pole 
% tz = (4*2^(1/2))/13 + 2/65; % 100*2^(1/2) + 10;

% Plant TF, 'Gp'
Gp = minreal(C/s^2)

% Chosen Youla Parameter, 'Y' -> Y(0) = 0
Ys = minreal(((K*s^2)*(tz*s + 1)/((s^2 + 2*Z*Wn*s + Wn^2)*(tp*s + 1)*(tpx*s + 1)^2)),1e-04)

% Complementary Sensitivity TF, 'T' -> T(0) = 1 (1st interpolation
% condition)
T = minreal((Ys*Gp),1e-04)

% Sensitivity TF, 'S'
S = minreal((1-T),1e-04)

% Controller TF, 'Gc'
Gc = minreal((Ys/S),1e-04)

% Return Ratio, 'L'
L = minreal((Gc*Gp),1e-04)

GpS = minreal((Gp*S),1e-04)

% Internal stability check
Ys_stability = isstable(Ys)
T_stability = isstable(T)
S_stability = isstable(S)
GpS_stability = isstable(GpS)

M2 = 1/getPeakGain(S) % M2-margin
BW = bandwidth(T) % Bandwidth of the closed-loop
AE = getPeakGain(Ys) % Maximum actuator effort

figure(1)
bodemag(Ys, S, T);
legend('Ys','S','T');

%% Simulation

load Gcs.mat

Gp = minreal(tf_coord, 1e-04);
Gc = minreal([Gc1 0 0 0 0 0; 0 Gc2 0 0 0 0; 0 0 Gc3 0 0 0; 0 0 0 Gc4 0 0; 0 0 0 0 Gc5 0; 0 0 0 0 0 Gc6], 1e-04);
Lu = minreal(Gc * Gp, 1e-04);
Ly = minreal(Gp * Gc, 1e-04);
Y = minreal(inv(eye(6) + Lu) * Gc); 
Ty = minreal(inv(eye(6) + Ly) * Ly); 
Sy = minreal(inv(eye(6) + Ly), 1e-04);
Su = minreal(inv(eye(6) + Lu), 1e-04);

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



