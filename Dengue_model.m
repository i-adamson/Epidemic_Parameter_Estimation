%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%parameters%%

b1 = 104;     %infection rate of strain 1
b2 = 110;   %infection rate of strain 2
g =  52;   %recovery rate
m = 1/65;    %demography rate
a = 2;    %temporary immunity rate
p = 0.1;    %ratio of contribution to force of initial infection (both strains)
p1 = 0.01;   %ratio of contribution to force of secondary infection of 1
p2 = 0.01;   %ratio of contribution to force of secondary infection of 2


%initial conditions
S0 = 10000;
I10 =20;
I20 = 20;
R10 = 0;
R20 = 0;
S10 = 0;
S20 = 0;
I120 = 10;
I210 = 10;
R0 = 0;
N = S0+I10+I20+R10+R20+S10+S20+I120+I210+R0;

t = 0:0.1:500;

%%model%%
%y = [S I1 I2 R1 R2 S1 S2 I12 I21 R];

[S,I1,I2,R1,R2,S1,S2,I12,I21,R] = two_strain(b1,b2,g,m,a,p,p1,p2,S0,I10,I20,R10,R20,S10,S20,I120,I210,R0,N,t);

% %plotting all the variables separately
figure
plot(S, 'LineWidth',2)
hold on, grid on
plot(I1, 'LineWidth',2)
plot(I2, 'LineWidth',2)
plot(R1, 'LineWidth',2)
plot(R2, 'LineWidth',2)
plot(S1, 'LineWidth',2)
plot(S2, 'LineWidth',2)
plot(I12, 'LineWidth',2)
plot(I21, 'LineWidth',2)
plot(R, 'LineWidth',2)
xlabel('Time (Years)','fontsize',14,'Interpreter','latex')
ylabel('Population (people)','fontsize',14,'Interpreter','latex')
legend('S','I1','I2','R1','R2','S1','S2','I12','I21','R', 'fontsize',14,'Interpreter','latex') 

%%plotting compiled S, I, R
figure
plot(S+S1+S2, 'LineWidth',2)
hold on, grid on
xlabel('Time (years)','fontsize',14,'Interpreter','latex')
ylabel('$S = S+S_1+S_2$','fontsize',14,'Interpreter','latex')
legend('S', 'fontsize',14,'Interpreter','latex')
% 
% 
figure
plot(I1+I2+I12+I21, 'LineWidth',2, 'Color', [0.88 0.45 0.08])
hold on, grid on
xlabel('Time (years)','fontsize',14,'Interpreter','latex')
ylabel('$I = I_1+I_2+I_{21}+I_{12}$','fontsize',14,'Interpreter','latex')
legend('I', 'fontsize',14,'Interpreter','latex')
% 
figure
plot(R+R1+R2, 'LineWidth',2, 'Color', [0.16 0.68 0.27])
hold on, grid on
% % title('Dengue Multi-Strain Model','fontsize',14,'Interpreter','latex')
xlabel('Time (years)','fontsize',14,'Interpreter','latex')
ylabel('$R = R+R_1+R_2'$,'fontsize',14,'Interpreter','latex')
legend('R', 'fontsize',14,'Interpreter','latex')

times = t(:);
counts = (I1+I2+I12+I21);
Idata = table(times, counts, 'VariableNames',{'Time', 'I'});

save('I.mat', 'Idata');

function[S, I1, I2, R1, R2, S1, S2, I12, I21, R] = two_strain(b1,b2,g,m,a,p,p1,p2,S0,I10,I20,R10,R20,S10,S20,I120,I210,R0,N,t)
dydt = @(t,y)[-(b1/N)*y(1)*(y(2)+p*y(9))-(b2/N)*y(1)*(y(3)+p*y(8))+m*(N-y(1));
               (b1/N)*y(1)*(y(2)+p*y(9))-(g+m)*y(2);
               (b2/N)*y(1)*(y(3)+p*y(8))-(g+m)*y(3);
               g*y(2)-(a+m)*y(4);
               g*y(3)-(a+m)*y(5);
               -(b2/N)*y(6)*(y(3)+p2*y(8))+a*y(4)-m*y(6);
               -(b2/N)*y(7)*(y(2)+p1*y(9))+a*y(5)-m*y(7);
               (b2/N)*y(6)*(y(3)+p2*y(8))-(g+m)*y(8);
               (b2/N)*y(7)*(y(2)+p1*y(9))-(g+m)*y(9);
               g*(y(8)+y(9))-m*y(10)];
  
opts = odeset('abstol',1e-25, 'reltol', 1e-25, 'maxstep', 0.01);

[~,y] = ode15s(dydt, t, [S0 I10 I20 R10 R20 S10 S20 I120 I210 R0], opts);

%parse outputs
S = y(:,1);
I1 = y(:,2);
I2 = y(:,3);
R1 = y(:,4);
R2 = y(:,5);
S1 = y(:,6);
S2 = y(:,7);
I12 = y(:,8);
I21 = y(:,9);
R = y(:,10);

end
