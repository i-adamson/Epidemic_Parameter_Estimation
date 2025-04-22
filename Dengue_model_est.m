%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading the data to use from the PMCMC samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dat = readtable('C:/Users/nwwt55/Documents/post_dat.csv');
aux_params = readtable('C:/Users/nwwt55/Documents/aux_params.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N= 10000;

sample1_data = dat(dat.np==0,:);
sample1_aux = aux_params(aux_params.np==0,:);

g = table2array(sample1_aux(:,2));   %recovery rate
a = table2array(sample1_aux(:,3));    %temporary immunity rate

b1 = @(t) interp1(sample1_data.time, sample1_data.beta1, t, 'linear', 'extrap');
b2 = @(t) interp1(sample1_data.time, sample1_data.beta2, t, 'linear', 'extrap');
s = @(t) interp1(sample1_data.time, sample1_data.S, t, 'linear', 'extrap');
i1 = @(t) interp1(sample1_data.time, sample1_data.I1, t, 'linear', 'extrap');
i2 = @(t) interp1(sample1_data.time, sample1_data.I2, t, 'linear', 'extrap');



t = table2array(sample1_data(:,2));

%initial conditions
S0 = 10;
I10 = 10;
I20 = 10;
R10 = 0;
R20 = 0;
S10 = 10;
S20 = 10;
I120 = 0;
I210 = 0;
R0 = 0;

[S,I1,I2,R1,R2,S1,S2,I12,I21,R] = two_strain2(b1,b2,g,a,s,i1,i2,S0,I10,I20,R10,R20,S10,S20,I120,I210,R0,N,t);


figure
plot(S, 'LineWidth',2)
% hold on, grid on
% plot(I1, 'LineWidth',2)
% plot(I2, 'LineWidth',2)
% plot(R1, 'LineWidth',2)
% plot(R2, 'LineWidth',2)
% plot(S1, 'LineWidth',2)
% plot(S2, 'LineWidth',2)
% plot(I12, 'LineWidth',2)
% plot(I21, 'LineWidth',2)
% plot(R, 'LineWidth',2)
% xlabel('Time (Years)','fontsize',14,'Interpreter','latex')
% ylabel('Population (people)','fontsize',14,'Interpreter','latex')
% legend('S','I1','I2','R1','R2','S1','S2','I12','I21','R', 'fontsize',14,'Interpreter','latex')

%%then the goal would be to run the chaos test on the number of infected
%%individuals for each sample.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The adapted Multi-strain model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S, I1, I2, R1, R2, S1, S2, I12, I21, R] = two_strain2(b1, b2, g, a, s, i1, i2, S0, I10, I20, R10, R20, S10, S20, I120, I210, R0, N, t)
    % Define the system of ODEs
    dydt = @(t, y) [-(b1(t)*s(t)*i1(t)*y(1))/((y(1)+y(7))*N)-(b2(t)*s(t)*i2(t)*y(1))./((y(1)+y(6))*N);
                     (b1(t)*s(t)*i1(t)*y(1))/((y(1)+y(7))*N)-g*y(2);
                     (b2(t)*s(t)*i2(t)*y(1))/((y(1)+y(6))*N)-g*y(3);
                     g*y(2)-a*y(4);
                     g*y(3)-a*y(5);
                     a*y(4)-(b2(t)*s(t)*i2(t)*y(6))/((y(1)+y(6))*N);
                     a*y(5)-(b1(t)*s(t)*i1(t)*y(7))/((y(1)+y(7))*N);
                     (b2(t)*s(t)*i2(t)*y(6))/((y(1)+y(6))*N)-g*y(8);
                     (b1(t)*s(t)*i1(t)*y(7))/((y(1)+y(7))*N)-g*y(9);
                     g*(y(8)+y(9))];
    
    opts = odeset('abstol',1e-16, 'reltol', 1e-16, 'maxstep', 0.01);

    % Solve the ODE system
    [~, y] = ode23tb(dydt, t, [S0 I10 I20 R10 R20 S10 S20 I120 I210 R0], opts);

    % Parse the outputs (returning the columns of y)
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
