% Astrobee Model

%% A Matrix

load Matrices/A_matrix.mat
A = A_matrix

%% B Matrix: Original

load Matrices/B_original.mat
B = B_original

%% B Matrix: Stowed

load Matrices/B_stowed.mat
B = B_stowed

%% B Matrix: Deployed

load Matrices/B_deployed.mat
B = B_deployed

%% Full-State Feedback

Cf = eye(12);

Df = [zeros(12, 6)];

sys_full = ss(A, B, Cf, Df);

tf_full = tf(sys_full)

%% Coordinate Feedback

Cc = [zeros(6, 12)];
Cc(1:6, 1:6) = eye(6);

Dc = [zeros(6, 6)];

sys_coord = ss(A, B, Cc, Dc);

tf_coord = tf(sys_coord)