%计算不同腿长下适合的K矩阵，再进行多项式拟合，得到2*6矩阵每个参数对应的多项式参数
tic
j=1;
leg=0.04:0.001:0.13;
for i=leg
    k=get_LQR_K_var(i);
    k11(j) = k(1,1);
    k12(j) = k(1,2);
    k13(j) = k(1,3);
    k14(j) = k(1,4);
    k15(j) = k(1,5);
    k16(j) = k(1,6);

    k21(j) = k(2,1);
    k22(j) = k(2,2);
    k23(j) = k(2,3);
    k24(j) = k(2,4);
    k25(j) = k(2,5);
    k26(j) = k(2,6);
    j=j+1;
end
a11=polyfit(leg,k11,3);
a12=polyfit(leg,k12,3);
a13=polyfit(leg,k13,3);
a14=polyfit(leg,k14,3);
a15=polyfit(leg,k15,3);
a16=polyfit(leg,k16,3);

a21=polyfit(leg,k21,3);
a22=polyfit(leg,k22,3);
a23=polyfit(leg,k23,3);
a24=polyfit(leg,k24,3);
a25=polyfit(leg,k25,3);
a26=polyfit(leg,k26,3);

a = [a11; a12; a13; a14; a15; a16; a21; a22; a23; a24; a25; a26]; % 2×6 矩阵
for i = 1:size(a,1)
    fprintf('float a%d%d[4] = {%.5f, %.5f, %.5f, %.5f};\n', ceil(i/6), mod(i-1,6)+1, a(i,:));
end


