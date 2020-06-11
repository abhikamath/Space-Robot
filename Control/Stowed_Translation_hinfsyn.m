%% Astrobee Stowed H-infinity Control Design for Translation (hinfsyn)

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

%% Hinf

w=logspace(-3,3);

s=tf('s');

% 
Gtf=Gp;

Gtf=minreal(Gtf);
Gtfss=ss(Gtf);
[Ag,Bg,Cg,Dg]=ssdata(Gtfss);
[Ag,Bg,Cg,Dg]=minreal(Ag,Bg,Cg,Dg);

% Hinf does not work well, as Hinf norm is not defined on the jw-axis, with poles and zeros of Gp on the jw-axis,
% So we use a simple shifting trick to make the Ag BIBO stable.

% Ag1=Ag-0.05*eye(6);

%define Wd (Robust uncertainty filter)

Wdtf=[makeweight(0.01,5,1000) 0 0 0 0 0; 0 makeweight(0.01,5,1000) 0 0 0 0; 0 0 makeweight(0.01,5,1000) 0 0 0; 0 0 0 makeweight(0.01,5,1000) 0 0; 0 0 0 0 makeweight(0.01,5,1000) 0; 0 0 0 0 0 makeweight(0.01,5,1000)];

%% define a new robust uncertanity filter

Wdss=ss(Wdtf);
[Ad,Bd,Cd,Dd]=ssdata(Wdss);
[Ad,Bd,Cd,Dd]=minreal(Ad,Bd,Cd,Dd);

%define Wu (Actuator filter)

Wutf=[0.01 0 0; 0 0.01 0; 0 0 0.01];
Wuss=ss(Wutf);
[Au,Bu,Cu,Du]=ssdata(Wuss);
[Au,Bu,Cu,Du]=minreal(Au,Bu,Cu,Du);

%define Wp (Performance filter)

Wptf=[makeweight(1000,0.5,0.1) 0 0 0 0 0; 0 makeweight(1000,0.5,0.1) 0 0 0 0; 0 0 makeweight(1000,0.5,0.1) 0 0 0; 0 0 0 makeweight(1000,0.5,0.1) 0 0; 0 0 0 0 makeweight(1000,0.5,0.1) 0; 0 0 0 0 0 makeweight(1000,0.5,0.1)];

%% define a new performance filter

Wpss=ss(Wptf); [Ap,Bp,Cp,Dp]=ssdata(Wpss);
[Ap,Bp,Cp,Dp]=minreal(Ap,Bp,Cp,Dp);

%compute augmented plant
[A,B1,B2,C1,C2,D11,D12,D21,D22]=...
    augss(Ag,Bg,Cg,Dg,...
          Ap,Bp,Cp,Dp,...
          Au,Bu,Cu,Du,...
          Ad,Bd,Cd,Dd);

% P_1 = [A,B1,B2,C1,C2,D11,D12,D21,D22];
B = [B1 B2];
C = [C1; C2];
D = [D11 D12; D21 D22];

P = ss(A, B, C, D);

%Compute Hinf controller
[K, CL, Gamma]=...
    hinfsyn(P, 6, 3);

% We have to shift back the A matrix of the controller to make this
% approach work properly.
% acp1=acp+0.05*eye(size(acp));

%% compute Tzw using lft command

acp = K.A; bcp = K.B; ccp = K.C; dcp = K.D;

% open loop
[aly,bly,cly,dly]=series(acp,bcp,ccp,dcp,Ag,Bg,Cg,Dg);
[aly,bly,cly,dly]=minreal(aly,bly,cly,dly);

%closed loop
[at,bt,ct,dt]=cloop(aly,bly,cly,dly,-1);
[at,bt,ct,dt]=minreal(at,bt,ct,dt);
as=at;
bs=bt;
cs=-ct;
ds=eye(6)-dt;

[ay,by,cy,dy]=series(as,bs,cs,ds,acp,bcp,ccp,dcp);
[ay,by,cy,dy]=minreal(ay,by,cy,dy);

figure(1);
sigma(aly,bly,cly,dly)
hold on
sigma(inv(Wdtf(2,2)),w)
hold on
sigma(Wptf(1,1),w)
legend('Ly singular values','1/Wd','Wp')


figure(2)
sigma(at,bt,ct,dt,w)
hold on
sigma(as,bs,cs,ds,w)
hold on
sigma(inv(Wdtf(2,2)),w)
hold on
sigma(inv(Wptf(1,1)),w)
legend('T','S','1/Wd','1/Wp')
hold off


figure(3)
sigma(ay,by,cy,dy,w)
legend('Youla singular values')

%closed loop time response
t=0:0.01:10;
figure(4)
step(at,bt,ct,dt,1,t)
legend('step responses to the first input')
figure(5)
step(at,bt,ct,dt,2,t)
legend('step responses to the second input')


figure(6)
sigma(Gtf)
legend('Plant singular values')


%% By students, compute the maximum sigular value of Tzw from lft computation

figure (7)
%sigma(Tzw)
hold on
sigma(CL/Gamma)
legend('Inf norm of Tzw')

Gc_ss=ss(acp,bcp,ccp,dcp);
Gc_tf=tf(Gc_ss);
Gc=minreal(Gc_tf);


% save('Gc_Hinf_ss.mat','Gc_ss')  % function form
% 
% save('Gc_Hinf_tf.mat','Gc_tf')  % function form

% Y = minreal(inv(eye(3) + Lu) * Gc); 
% Ty = minreal(inv(eye(6) + Ly) * Ly);
% Sy = minreal(inv(eye(6) + Ly), 1e-02);

%% Save Controller

Gc_Stowed_Translation_hinfsyn = Gc;
save('Matrices/Gc_Stowed_Translation_hinfsyn.mat', 'Gc_Stowed_Translation_hinfsyn');




