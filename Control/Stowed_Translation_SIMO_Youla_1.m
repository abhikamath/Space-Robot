%% Astrobee Stowed Translation SIMO Youla Controller (I/O Channel 1)

% Run sections sequentially

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

tf_full = minreal(tf(sys_full));

tf_translation = minreal([tf_full(1:3, 1:3); tf_full(7:9, 1:3)]);

tf_full_sym = simplify(tf2sym(tf_full));
disp('tf_full_sym = ');
pretty(tf_full_sym);

tf_translation_sym = simplify(tf2sym(tf_translation));
disp('tf_translation_sym = ');
pretty(tf_translation_sym);

%% Smith-McMillan Form

% help smform

s=tf('s');

Gp = 500/7939*[1/s; 0; 0; 1/s^2; 0; 0];

Mp = minreal(smform(Gp));

syms s

Gp_sym = tf2sym(Gp);

Mp_sym = tf2sym(Mp);

[UL_sym, H] = hermiteForm(Gp_sym * s^2)

UR_sym = (UL_sym * Gp_sym)\Mp_sym

% UL_sym = [7939/(500*s), 0, 0,  0,  0,  0;
%           0, 7939/(500*s), 0,  0,  0,  0;
%           0, 0, 7939/(500*s),  0,  0,  0;
%           1/s, 0, 0, -1,  0,  0;
%           0, 1/s, 0,  0, -1,  0;
%           0, 0, 1/s,  0,  0, -1];
% 
% UR_sym = eye(3);
      
% UL_sym = [0, 0, 0,  1,  0,  0;
%           0, 0, 0,  0,  1,  0;
%           0, 0, 0,  0,  0,  1;
%           -1, 0, 0,  s,  0,  0;
%           0, -1, 0,  0,  s,  0;
%           0, 0, -1,  0,  0,  s]
%       
% UR_sym = [7939/500, 0, 0;
%           0, 7939/500, 0;
%           0, 0, 7939/500]

% UL_sym = [0 1; 1 -s];
% UR_sym = 7939/500;

disp('UL_sym = ');
pretty(UL_sym);

% Mp_sym = tf2sym(Mp);
disp('Mp_sym = ');
pretty(Mp_sym);

disp('UR_sym = ');
UR_sym

UL = sym2tf(UL_sym);
UR = double(UR_sym);

%% Interpolation Conditions

% Run this section first to calculate 'tz' to ensure that the second interpolation condition is satisfied

% d^k(T)/ds^k|(s=0) = 0, where k = 1 (since there is a double unstable pole
% (multiplicity ap = 2) in the plant at s = 0; k = ap - 1) -> 2nd
% interpolation condition

% Constants & Design Parameters

C = 1; % Constant
Wn = 1; % Natural Frequency of the Control System
K = Wn^2/C; % Controller Gain
Z = 2^-0.5; % Damping Ratio
tp = 1/(10*Wn); % Time constant (of the first included pole)
tpx = 0.5; % Time constant (of the pole included to drop Youla at high frequencies)
% tzx = 999;

syms s tz

T_eqn = ((K*C)*(tz*s + 1)/((s^2 + 2*Z*Wn*s + Wn^2)*(tp*s + 1)*(tpx*s + 1)^2));
dT_eqn = diff(T_eqn,s);
eqn = subs(dT_eqn,s,0) == 0;
tz = double(solve(eqn,tz))

%% Control Design

s = tf('s');

% Plant TF, 'Gp'
G = minreal(C/s^2, 1e-04) % Nonzero terms of Mp

% Chosen Youla Parameter, 'Y' -> Y(0) = 0
Ys = minreal(((K*s^2)*(tz*s + 1)/((s^2 + 2*Z*Wn*s + Wn^2)*(tp*s + 1)*(tpx*s + 1)^2)),1e-04)

% Complementary Sensitivity TF, 'T' -> T(0) = 1 (1st interpolation
% condition)
T = minreal((Ys*G),1e-04)

% Sensitivity TF, 'S'
S = minreal((1-T),1e-04)

% Controller TF, 'Gc'
Gc_term = minreal((Ys/S),1e-04)

% Return Ratio, 'L'
L = minreal((Gc_term*G),1e-04)

GS = minreal((G*S),1e-04)

% Internal stability check
Y_stability = isstable(Ys)
T_stability = isstable(T)
S_stability = isstable(S)
GS_stability = isstable(GS)

M2 = 1/getPeakGain(S) % M2-margin
BW = bandwidth(T) % Bandwidth of the closed-loop
AE = getPeakGain(Ys) % Maximum actuator effort

figure(1)
bodemag(Ys, S, T);
legend('Ys','S','T');

% Gc = minreal([tf(Gc_term) 0 0 0 0 0; 0 tf(Gc_term) 0 0 0 0; 0 0 tf(Gc_term) 0 0 0]);
% 
% Gc_sym = expand(tf2sym(Gc));
% disp('Gc_sym = ');
% pretty(Gc_sym);

My = minreal([Ys 0 0 0 0 0], 1e-04);
Mt = minreal(Mp * My, 1e-04);
% Mt = minreal(T * eye(6), 1e-04);
Y = minreal(UR * My * UL, 1e-04);

%% Simulation

% Y = minreal(inv(eye(3) + Lu) * Gc); 
% Ty = minreal(inv(eye(6) + Ly) * Ly);
% Sy = minreal(inv(eye(6) + Ly), 1e-04);

Ty = minreal(inv(UL) * Mp * My * UL, 1e-04);

Sy = minreal(eye(size(Ty)) - Ty, 1e-04);

Gc = minreal(UR * inv(eye(size(My * Mp)) - (My * Mp)) * My * UL, 1e-04);

SyGp = minreal(inv(UL) * (eye(size(Mp * My)) - (Mp * My)) * Mp * inv(UR), 1e-04);

% MIMO Internal Stability Check:
Ty_stability = isstable(Ty)
Sy_stability = isstable(Sy)
Gc_stability = isstable(Gc)
SyGp_stability = isstable(SyGp)

SV_Gp = sigma(Gp); 
k_Gp = max(max(SV_Gp))/min(min(SV_Gp)) % condition-number check for Gp

SV_Gc = sigma(Gc); 
k_Gc = max(max(SV_Gc))/min(min(SV_Gc)) % condition-number check for Gc

Lu = minreal(Gc * Gp, 1e-04);
Ly = minreal(Gp * Gc, 1e-04);
Su = minreal(inv(eye(size(Lu)) + Lu), 1e-04);

figure
step(Ty);

figure
step(Y);

figure
sigma(Y, Ty, Sy, Su)
[l, hObj] = legend('$Y$', '$T_{y}$', '$S_{y}$', '$S_{u}$','Interpreter','latex','FontSize', 30);
set(l,'string',{'$Y$', '$T_{y}$', '$S_{y}$', '$S_{u}$'});
hL = findobj(hObj,'type','line');
set(hL,'linewidth', 5); 

% figure
% sigma(Gc, Gp, Ly, Y)
% [l, hObj] = legend('$G_{c}$', '$G_{p}$', '$L_{y}$', '$Y$','Interpreter','latex','FontSize', 12);
% set(l,'string',{'$G_{c}$', '$G_{p}$', '$L_{y}$', '$Y$'});
% hL = findobj(hObj,'type','line');
% set(hL,'linewidth', 2);
% 
% figure
% sigma(Gc, Gp, Y)
% [l, hObj] = legend('$G_{c}$', '$G_{p}$', '$Y$','Interpreter','latex','FontSize', 12);
% set(l,'string',{'$G_{c}$', '$G_{p}$', '$Y$'});
% hL = findobj(hObj,'type','line');
% set(hL,'linewidth', 2);
% 
% figure
% sigma(Ly, Sy, Ty)
% [l, hObj] = legend('$L_{y}$', '$S_{y}$', '$T_{y}$','Interpreter','latex','FontSize', 12);
% set(l,'string',{'$L_{y}$', '$S_{y}$', '$T_{y}$'});
% hL = findobj(hObj,'type','line');
% set(hL,'linewidth', 2);
% 
% figure
% sigma(Sy, Su)
% [l, hObj] = legend('$S_{y}$', '$S_{u}$','Interpreter','latex','FontSize', 12);
% set(l,'string',{'$S_{y}$', '$S_{u}$'});
% hL = findobj(hObj,'type','line');
% set(hL,'linewidth', 2); 

%% Save Controller

Gc_Stowed_Translation_SIMO_Youla_1 = Gc;
save('Matrices/Gc_Stowed_Translation_SIMO_Youla_1.mat', 'Gc_Stowed_Translation_SIMO_Youla_1');





