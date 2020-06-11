%% Astrobee Stowed Translation SIMO Youla Control Script

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

Gp = tf_translation;

load Matrices/Gc_Stowed_Translation_SIMO_Youla_1.mat
load Matrices/Gc_Stowed_Translation_SIMO_Youla_2.mat
load Matrices/Gc_Stowed_Translation_SIMO_Youla_3.mat

Gc = [Gc_Stowed_Translation_SIMO_Youla_1; Gc_Stowed_Translation_SIMO_Youla_2; Gc_Stowed_Translation_SIMO_Youla_3];

Lu = minreal(Gc * Gp, 1e-04);
Ly = minreal(Gp * Gc, 1e-04);
Y = minreal(inv(eye(size(Lu)) + Lu) * Gc, 1e-04); 
Ty = minreal(inv(eye(size(Ly)) + Ly) * Ly, 1e-04); 
Sy = minreal(inv(eye(size(Ly)) + Ly), 1e-04);
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