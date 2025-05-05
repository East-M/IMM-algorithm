clear all;
clc;
clf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �����������״��һ��ֱ�߷��еķɻ���Ŀ��ĸ����˲�
% �״������������Բ��ɨ���״�,�������״�Ӧ�����,�����һ���򻯰�����,��Ϊx,y����Ĳ����Ƕ�����
% Ŀ����������ٻ��ȼ��ٷ��У���ʼ����\��ʼ�ٶ�\��ʼ���پ�����
% �˲�ģ�ͣ�Kalman�˲�������ģ�ͣ���CVģ�ͣ�
% ע�����ע��������������ã�����̫С����Ԥ��Ϊ�����Ƚ�ƽ������Ҳ�п��ܷ�ɢ������̫����������Ϊ����ƽ��Ч������
% ��ҵҪ��:
% (1) ������,���й۲�Ч��
% (2) �ı��������sig_w��ȡֵ,������С,�۲��˲���ƽ���̶�,������ԭ��
% (3) �ı��ʱ���ٶ�ax��ay��ȡֵ,����ax=0.1,�۲��˲����Ƿ�ɢ,������ԭ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ŀ����״����
T = 1;     % ɨ������s
k = 1000;  % �۲�����
x0 = 0;    % Ŀ���ʼλ��(0km,500km)
y0 = 500e3;
vx = 300;  % Ŀ���ʱ�ٶ�(m\s)
vy = 1;
ax = 0;    % Ŀ���ʱ���ٶ�(m\s^2)
ay = 0;
sig_x = 100;            % ��ྫ��m
sig_y = 100;            % ��Ǿ���rad\
t = (0:k-1) * T;        % ʱ�伯��

% ������ʵ�켣
x=  x0 + vx*t + 0.5*ax*t.^2;
y = y0 + vy*t + 0.5*ay*t.^2;

% �����״�ԭʼ��������(measure)
xm = x  + normrnd(0,sig_x,1,k);
ym = y  + normrnd(0,sig_y,1,k);

% Kalman�˲�����
n = 4;              % ״̬����ά��(x xv y yv)
m = 2;              % �������ά��(xm,ym)
Phi = [ 1 T 0 0;    % ת�ƾ���4x4,Ҳ����Ԥ�����
        0 1 0 0;
        0 0 1 T; 
        0 0 0 1 ];    
G = [ T/2 0;        % ��������
      1   0;
      0 T/2; 
      0    1]; 
H=[ 1 0 0 0;        % �������,[mxn]ά
    0 0 1 0];  
sig_w = [1e-2];     % �����������Ʋ���,�ò����ǳ���Ҫ,�������۲�Ч��.
Q = [sig_w^2,    0; % ��������Э����
    0,   sig_w^2];
R = [sig_x^2        0; 
     0         sig_y^2]; %����Э����
X0 = [x(2)  (x(2)-x(1))/T   y(2) (y(2)-y(1))/T]'; %�˲���ʼֵ,�����ֵõ��ٶ�
P0 = [ sig_x^2        sig_x^2/T    0    0;        % ��ʼ�˲�Э����,���ݲ�ֹ�ʽ�õ�,��Ϊxy��������ֲ�
      sig_x^2/T     2*sig_x^2/T^2  0    0; 
      0    0    sig_y^2      sig_y^2/T; 
      0    0    sig_y^2/T    2*sig_y^2/T^2];
Z = [xm; ym]; %����ʸ��

%% Kalman�˲����ļ���ģ��
for i=1:k
   [X1,P1] = Fun_KF_Predict(X0,P0,Phi,Q,G);     % Ԥ��
   [XX,PP] = Fun_KF_Update(X1,P1,Z(:,i),H,R);   % ����(��ʵ���״����ݴ�����Ҫ�������������м��жϵ㼣z�Ƿ�����Ԥ�Ⲩ��,���������ݹ���,����Ϊ�����,�������Ӳ����龯)
    
    X0=XX; %Ϊ�˱�֤ѭ��,��Ҫ���¸�ֵ
    P0=PP;
    
    %������,����������˲�ֵ
    xf(i)  = XX(1,1); 
    xvf(i) = XX(2,1);
    yf(i)  = XX(3,1);
    yvf(i) = XX(4,1);
    
    px(i) = PP(1,1);
    py(i) = PP(3,3);
end

%% ���¿�ʼ������
figure(1)
% ������ͼ
plot(x/1e3,y/1e3,':','LineWidth',1);   % ��ʵ
hold on;
plot(xm/1e3,ym/1e3,'g.','LineWidth',2);% ����
plot(xf/1e3,yf/1e3,'k','LineWidth',2); % �˲�
hold off;
xlabel('X(km)');
ylabel('Y(km)');
title('�˲�����ͼ');
legend('��ʵ','�۲�','�˲�');
axis tight

figure(2)
% �״�P����ʾ
rou =sqrt( x.^2+y.^2);
theta = atan2(y,x);
rou_m = sqrt(xm.^2+ym.^2);
theta_m = atan2(ym,xm);
rou_f = sqrt(xf.^2+yf.^2);
theta_f = atan2(yf,xf);
polar(theta,rou,':');%��ʵ
hold on
polar(theta_m,rou_m,'g.');%��ʵ
polar(theta_f,rou_f,'k');%��ʵ
hold off
title('�״�PPI��ʾ');

figure(3)
% ��x�����˲����
plot(t,xm-x,'g:','LineWidth',2);
hold on
plot(t,xf-x,'k','LineWidth',2);
plot(t,3*sqrt(px),'r','LineWidth',2);
plot(t,-3*sqrt(px),'r','LineWidth',2);
hold off
legend('�۲����','�˲����','3\sigma���������');
xlabel('ʱ��(s)');
ylabel('Xλ�����(m)');
title('X�����˲����');

figure(4)
% ��x������Ƶ��ٶ�
plot(t,xvf,'k-','linewidth',2); %����ֵ
hold on
plot(t,vx,'b'); %��ʵֵ
hold off
xlabel('ʱ��(s)');
ylabel('�ٶ�(m/s)');
title('X��������ٶ�(m/s)');

figure(5)
% ��y������Ƶ��ٶ�
plot(t,yvf,'k-','linewidth',2); %����ֵ
hold on
plot(t,vy,'b'); %��ʵֵ
hold off
xlabel('ʱ��(s)');
ylabel('�ٶ�(m/s)');
title('Y��������ٶ�(m/s)');

