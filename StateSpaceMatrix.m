%定义参量
mb = 0.128;  % mass of the float-based body (kg)
mw = 0.15;  % mass of the wheel (kg)
rw = 0.035;  % radius of the wheel (m)
d = 0.15;    % 左右两轮的距离(m)
h = 0.15;   % height of the body (m)
w = 0.1;   % width of the body (m)

l = 0.07;   % length to the 质心 (m)
g = 9.81;  % gravitational acceleration (m/s²)

jp = mb*l^2/3;
jd = mb*d^2/12;
I = mw*rw^2/2;

Qeq = jp*mb + (jp+mb*l^2)*(2*mw+2*l/rw^2);

a23 = -1*(mb^2*l^2*g)/Qeq;
a43 = (mb*l*g*(mb+2*mw+2*l/rw^2))/Qeq;
b21 = (jp+mb*l^2+mb*l*rw)/(Qeq*rw);
b22 = b21;
b41 = -1*(mb*l/rw+mb+2*mw+2*I/rw^2)/Qeq;
b42 = b41;
b61 = 1/(rw*(mw*d+I*d/rw^2+2*jd/d));
b62 = -b61;

% Matrix A
A = zeros(6, 6);
A(1,2)=1;
A(2,3)=a23;
A(3,4)=1;
A(4,3)=a43;
A(5,6)=1;

% Matrix B
B = zeros(6, 2);
B(2,1)=b21;
B(2,2)=b22;
B(4,1)=b41;
B(4,2)=b42;
B(6,1)=b61;
B(6,2)=b62;


Q = diag([10,1,1,1,1,1]);
R = diag([1,2]);


% Calculate LQR gain matrix K
K_LQR = lqr(A,B,Q,R);
tK_LQR = 0.1 * [0:29]';
Simulation_K_LQR = repmat(K_LQR,[1 1 length(tK_LQR)]);
dataK_LQR.time=tK_LQR;
dataK_LQR.signals.values = Simulation_K_LQR;
dataK_LQR.signals.dimensions=[size(K_LQR, 1) size(K_LQR, 2)];








































