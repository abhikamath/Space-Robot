%% Astrobee Model

A = [zeros(12, 12)];
A(7:12, 1:6) = eye(6);

B = [zeros(12, 6)];
B(1:6, :) = eye(6);

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