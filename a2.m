%   循环搜索特征向量正交的时候
idx = 1;
[SP,SP_inv] = deal(zeros(181,1));
scale = -90:90;
for angle_degree = scale
    a =exp(1j.*2*pi*d*sind(angle_degree)/lambda).';     %构建信号导向矢量，用共轭转至全部加负号
    En = EV(:,1:end-M);                                 %用前面的几个小特征值的特征向量
    SP(idx) = (a'*En)*(En'*a);                          %利用前面讲的正交来判断结果
    SP_inv(idx) = 1/abs((a'*En)*(En'*a));               %用倒数翻转一下 变成峰值

    idx = idx + 1;
end
SP_db = db(SP);                                         %转换为dB
SP_inv_db = db(abs(SP_inv));                            %转换为dB

figure(8)
subplot(211);plot(scale,SP_db);
xlabel('入射角/(degree)');ylabel('空间谱/(dB)');
grid on;title('正交表示目标信号来向')
subplot(212);plot(scale,SP_inv_db);
xlabel('入射角/(degree)');ylabel('空间谱/(dB)');
grid on;title('用伪谱表示目标信号来向')
