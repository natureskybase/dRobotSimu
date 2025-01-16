% 此模型根据哈尔滨工业大学单腿模型进行建模
% 资料链接为(https://zhuanlan.zhihu.com/p/563048952)
function K = get_LQR_K(length)
%变量定义
% length =0.065;

%车体常量定义
L = 0.15;       %车体长
W = 0.1;        %车体宽
H = 0.03;        %车体高

g = 9.81 ;

RW = 0.035;     %驱动轮半径
LC = 0.01;      %机体质心到腿部关节的距离
MW = 0.24;      %驱动轮质量
MB = 0.128;     %车体质量
ML = 0.02;      %腿部质量

LW = length/2;  %虚拟摆杆重心到驱动轮转轴的距离
LB = length/2;  %虚拟摆杆重心到转轴的距离

IW = MW * RW^2;                 %驱动轮转动惯量
IL = ML*((LW+LB)^2 + H^2)/12.0; %摆杆绕质心转动惯量
IB = MB*(L^2+W^2)/12.0;         %机体绕质心的转动惯量

% todo 文案书写
% x:驱动轮位移
syms x(t)  theta(t)  phi(t)    % x矩阵变量
syms dd_theta  dd_x  dd_phi    % 二阶导
syms d_theta  d_x  d_phi       % 一阶导
syms T  Tp                     % u矩阵变量
syms theta0  x0  phi0          % 各变量初值

% 已经根据哈工大公式解得变量
NM = MB*diff(x + (LB + LW )*sin(theta)-LC*sin(phi),t,2);
N = NM + ML*diff(x + LB*sin(theta),t,2);
PM = MB*g + MB*diff((LB+LW)*cos(theta)+LC*cos(phi),t,2);
P = PM +ML*g+ML*diff(LB*cos(theta),t,2);

% 需求解的非线性方程组
eqn1 = diff(x,t,2) == (T -N*RW)/(IW/RW + MW*RW);
eqn2 = IL*diff(theta,t,2) == (P*LW + PM*LW)*sin(theta)-(N*LW+NM*LW)*cos(theta)-T+Tp;
eqn3 = IB*diff(phi,t,2) == Tp +NM*LC*cos(phi)+PM*LC*sin(phi);

eqn10 = subs(subs(subs(subs(subs(subs(subs(subs(subs(eqn1,diff(theta,t,2),dd_theta),diff(x,t,2),dd_x),diff(phi,t,2),dd_phi),diff(theta,t),d_theta),diff(x,t),d_x),diff(phi,t),d_phi),theta,theta0),x,x0),phi,phi0);
eqn20 = subs(subs(subs(subs(subs(subs(subs(subs(subs(eqn2,diff(theta,t,2),dd_theta),diff(x,t,2),dd_x),diff(phi,t,2),dd_phi),diff(theta,t),d_theta),diff(x,t),d_x),diff(phi,t),d_phi),theta,theta0),x,x0),phi,phi0);
eqn30 = subs(subs(subs(subs(subs(subs(subs(subs(subs(eqn3,diff(theta,t,2),dd_theta),diff(x,t,2),dd_x),diff(phi,t,2),dd_phi),diff(theta,t),d_theta),diff(x,t),d_x),diff(phi,t),d_phi),theta,theta0),x,x0),phi,phi0);

[dd_theta,dd_x,dd_phi] = solve(eqn10,eqn20,eqn30,dd_theta,dd_x,dd_phi);

A=subs(jacobian([d_theta,dd_theta,d_x,dd_x,d_phi,dd_phi],[theta0,d_theta,x0,d_x,phi0,d_phi]),[theta0,d_theta,d_x,phi0,d_phi,T,Tp],[0,0,0,0,0,0,0]);
A=double(A);
B=subs(jacobian([d_theta,dd_theta,d_x,dd_x,d_phi,dd_phi],[T,Tp]),[theta0,d_theta,d_x,phi0,d_phi,T,Tp],[0,0,0,0,0,0,0]);
B=double(B);

%theta   d_theta   x   d_x   phi   d_phi
Q=diag([100 1 500 10 5000 1]);
R=[240 0;0 25];

K=lqr(A,B,Q,R);

end
