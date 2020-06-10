% Astrobee Model

%% A Matrix

load Matrices/A_matrix.mat
A = A_matrix

%% B Matrix: Original

load Matrices/B_original.mat
B = B_original

%% Full-State Feedback

Cf = eye(12);

Df = [zeros(12, 6)];

sys_full = ss(A, B, Cf, Df);

tf_full = tf(sys_full);

syms s

tf_full_sym = simplify(Cf * inv(s * eye(12) - A) * B + Df);
pretty(tf_full_sym)

translation_full = [tf_full_sym(1:3, 1:3); tf_full_sym(7:9, 1:3)];
pretty(translation_full)

attitude_full = [tf_full_sym(4:6, 4:6); tf_full_sym(10:12, 4:6)];
pretty(attitude_full)

%% Coordinate Feedback

Cc = [zeros(6, 12)];
Cc(1:6, 1:6) = eye(6);

Dc = [zeros(6, 6)];

sys_coord = ss(A, B, Cc, Dc);

tf_coord = tf(sys_coord);

syms s

tf_coord_sym = simplify(Cc * inv(s * eye(12) - A) * B + Dc);
pretty(tf_coord_sym)

translation_coord = [tf_coord_sym(1:3, 1:3); tf_coord_sym(7:9, 1:3)];
pretty(translation_coord)

attitude_full = [tf_coord_sym(4:6, 4:6); tf_coord_sym(10:12, 4:6)];
pretty(attitude_coord)

