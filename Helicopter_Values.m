    % clear all

%% Table 1: Parameters of the system: Part 1
lcm = 0.015;   % Distance between the pivot point and the center of mas of helicopter [m]
mheli = 0.479; % Total moving mass of the helicopter [kg]
jeqp = 0.0172; % Moment of inertia about the pitch axis [kg-m^2]
jeqy = 0.021;  % Moment of inertia about the yaw axis 
g = 9.81;      % Acceleration due to gravity on planet earth [m-s^-2]

%% Table 2: Parameters of the system: Part 2
kpp = 0.0556;%0.063805;%0.0496; %0.050093; %0.0556;  % Torque constant on pitch axis from pitch motor/propeller [Nm/V]

kyy =0.21084;%0.21151 ;% 0.21084; % Torque constant on yaw axis from yaw motor/propeller [Nm/V]

kpy = 0.005;   % Torque constant on pitch axis from yaw motor/propeller [Nm/V]
kyp = 0.15;    % Toruqe constant on yaw axis from´pitch motor/propeller [Nm/V]
Bp = 0.01;%0.028807;%0.292865;%0.2199; %0.028807;%

%0.01;     % Damping friction factor about pitch axis [N/V]
%By = 0.21625;     % Damping friction factor about yaw axis [N/V]
By = 0.08;
%% Try to tune Bp, By, kpp, kyy 

%% Liniarization parameters
p_r_op = -10*pi/180;
y_r_op = pi/2;
wp_op = 0;
wy_op = 0;


Vmp_op = (kyy*mheli*g*cos(p_r_op)*lcm)/(kyy*kpp -kyp*kpy);
Vmy_op = kyp*Vmp_op/kyy;
%% A and B
%A
A1_3 = 1;

A2_4 = 1;

A3_1 = mheli*g*lcm*sin(p_r_op)/(jeqp +mheli*lcm*lcm);
A3_3 = -Bp/(jeqp +mheli*lcm*lcm);

A4_1_n = (kyp*Vmp_op -kyy*Vmy_op)*2*mheli*lcm*lcm*sin(p_r_op)*cos(p_r_op);
A4_1_d = (jeqy +mheli*cos(p_r_op)*cos(p_r_op)*lcm*lcm)^2;
A4_1 = A4_1_n/A4_1_d;

A4_4_d = (jeqy +mheli*cos(p_r_op)*cos(p_r_op)*lcm*lcm);
A4_4 = -By/A4_4_d;

A1 = [0, 0, A1_3, 0];
A2 = [0, 0, 0, A2_4];
A3 = [A3_1, 0, A3_3, 0];
A4 = [A4_1, 0, 0, A4_4];


%B
B3_d = (jeqp +mheli*lcm*lcm);
B3_1 = kpp/(B3_d);
B3_2 = -kpy/(B3_d);
B4_d = jeqy +mheli*cos(p_r_op)*cos(p_r_op)*lcm*lcm;
B4_1 = kyp/B4_d;
B4_2 = -kyy/B4_d;

B1 = [0,0];
B2 = [0,0];
B3 = [B3_1,B3_2];
B4 = [B4_1,B4_2];



A =[A1;A2;A3;A4];
B = [B1;B2;B3;B4];
C = [1,0,0,0;0,1,0,0];
D = [0,0;0,0];


%% Creating te discrete matrixes

Ts = 0.1;
sys_t = ss(A,B,C,D);
sys_d = c2d(sys_t,Ts);
Ad = sys_d.a;
Bd = sys_d.b;
Cd = sys_d.c;
Dd = sys_d.d;


G = eye(4);% Identity matrix of Number of states
H = zeros(2,4); %Nmber of measurements and states

sys = ss(Ad,[Bd G],Cd,[Dd H],Ts);

%% Kalman filter (K) calculation

%Q = diag([10,20,40,30]);  % Covariance matrixes Q for process noise( states) 
R = diag([25,30]);  % R for measuremets, 2 measurements, smaller beacuse we trust our sensors
Q = diag([20,35,25,15]);

%Because of matlab Kalman function we need to define Gw and Hw and Vk
[lqrgain,K,P] = kalman(sys,Q,R);
K

[k_lqr,s_lqr,e_lqr] = lqr(A,B,Q,R,0);
k_lqr

%% Kalman Filter with error integration
Ac_e = [A,zeros(4,2);C,zeros(2,2)];
Bc_e =[B;zeros(2,2)];
Cc_e = eye(size(Ac_e));
% Dc_e = D;
% Ts = 0.1;
% sys_t_e = ss(Ac_e,Bc_e,Cc_e,Dc_e);
% sys_d_e = c2d(sys_t,Ts);
% Ad_e = sys_d.a;
% Bd_e = sys_d.b;
% Cd_e = sys_d.c;
% Dd_e = sys_d.d;

Q_e = diag([1.0,1.4,0.8,1.5,0.4,0.5]);  % Covariance matrixes Q for process noise( states) 
R_e = diag([50,30]);  % R for measuremets, 2 measurements, smaller beacuse we trust our sensors

%Q_e=diag([20,35,25,15,5,10]);
%R_e=diag([50,30]);


[k_lqr_e,s_lqr_e,e_lqr_e] = lqr(Ac_e,Bc_e,Q_e,R_e,0);
k_lqr_e


