clc; clear;
close all; 

%% Code for the paper, 'Reformulating Parallel-Connected Lithium-Ion Battery Pack Dynamics with Interconnection Resistances as Ordinary Differential Equations' 

%% Load in the parameters for the batteries
battery_select = load('data.mat'); % load in a database containing battery model parametrs
bat_select_num = 3; % select what battery you want from the table, we have selected 3 for LFP and 4 for NMC
capacitance_Ah_all = [battery_select.data.capacity]; % Get the capacitances of all the batteries
R1_all = [battery_select.data.R1]; % Get the R1 values of all the batteries
C1_all = [battery_select.data.C1]; % Get the C1 values of all the batteries
r_all = [battery_select.data.R0]; % Get the series resistances of all the batteries


%% Open circuit voltage
ocv_coefs = battery_select.data(bat_select_num).OCV_coeff; % Load in the coefficients of the polynomial used to fit the OCV curves

ocv_points=battery_select.data(bat_select_num).OCV_points;

ocv_points= table2array(ocv_points);
N = 13;
ocv_coefs = polyfit(ocv_points(:,1),ocv_points(:,2),N);

nz = 1e2; % 
z_space= linspace(0,1,nz); %set of points for the soc to sample the ocv curve
for i = 1:nz
    OCV_plot(i) = ocv_eval(ocv_coefs,z_space(i)); % generate curves for the ocv
end
%% Model parameters
n_par = 10; % the number of cells in parallel

capacitance_Ah = capacitance_Ah_all(bat_select_num); % Define: capacitance in Ah.
capacitance = (3600*capacitance_Ah ); % Capacitance in As.

C_rate = 0.1; % Define the charging C-rate
%%
I = n_par*C_rate*capacitance_Ah; % Define the charging current

%%
R_nom = 1*1e-2;
r_nom = r_all(bat_select_num);  C_nom = C1_all(bat_select_num); % Set the nominal paraters for the 1st order circuit model of the battery
Q = capacitance*ones(n_par,1);
R = R_nom*ones(n_par,1); C = C_nom*ones(n_par,1); r =r_nom*ones(n_par,1); % Vectors of all the variables in the pack (as in stacked in parallel np times).

sd_vars = 1*1e-3; % standard deviation of the noise for the parameters
var_Q = sd_vars*capacitance; % Add noise to each of the parameters to simulate ageing.
var_R = sd_vars*R_nom*0;
var_F = sd_vars*R_nom;
var_r = sd_vars*r_nom;
var_C = sd_vars*C_nom;

for j = 1:n_par % add noise to each of the parametrs in the pack
    Q(j) = Q(j)+var_Q*randn;
    R(j) = R(j)+var_R*randn;
    r(j) = r(j)+var_r*randn;
    C(j) = C(j)+var_C*randn;
end

R_nom = R1_all(bat_select_num);
F = R_nom*ones(n_par,1);
F_k = F;

tau = 1./(F_k.*C); % Define the time constant

r(n_par) = r(n_par)+R_nom;
%% Compute functions
[new_A11,new_A12,new_A21,new_A22,new_m] = new_compute_A11_A12_A21_A22(n_par,R,C,Q,tau,r); % Define the A11, A12, A21, A22 matrices from the parallel pack project paper.

%% Define the initial conditions of the model
x0 = zeros(2*n_par,1);
x0 = 0.2*ones(2*n_par,1);
for j = 1:n_par % the initial condition for v_rc
    x0(2*j) = 0;
end

%% Run the simulation
t_f =  0.6 * 3600/(C_rate);
tspan = linspace(0, t_f);

global i_branch_store;
global i_branch_store_2;

tic
[t_ode_dae, x_ode_dae] = ode15s(@(t, x) parallel_model_dae6(t, x, I, new_A11, new_A12, new_A21, new_A22, ocv_coefs, R, r), tspan, x0);
time_dae_analytic = toc;

i_branch_values_analytic = i_branch_store;

x_ode = x_ode_dae;
t_ode = t_ode_dae;

tic
[t_ode_dae2, x_ode_dae2] = ode15s(@(t, x) parallel_model_dae5(t, x, I, new_A11, new_A12, new_A21, new_A22, ocv_coefs, R, r), tspan, x0);
time_dae_inverting = toc;

i_branch_values_inverse = i_branch_store_2;

error_x = norm(x_ode_dae2-x_ode_dae,inf);

x_ode2 = x_ode_dae2;
t_ode2 = t_ode_dae2;

size_t = size(x_ode2,1); % get the number of time steps of the simulation
z_ode = zeros(size_t,n_par); w_ode = zeros(size_t,n_par);
for i = 1:n_par % extact the state-of-charge and relaxation voltage from the simulation
    z_ode(:,i) = x_ode2(:,2*i-1);
    w_ode(:,i) = x_ode2(:,2*i);
end

nT = max(size(t_ode2));i_branch = zeros(n_par,nT);
phi = zeros(n_par,1); phi_init = zeros(n_par-1,1);

for j = 1:nT
    for i = 1:n_par
        OCV_ode(j,i) = ocv_eval(ocv_coefs,z_ode(j,i));
        volts(j,i) = OCV_ode(j,i)+ w_ode(j,i);
    end

    for i = 2:n_par
        phi_init(i-1) = volts(j,i)- volts(j,i-1);
    end
    phi =[I;phi_init];

    i_branch(:,j) = new_A22\phi; % Get the current going into each parallel branch
    current_sum(j) = sum(i_branch(:,j))/I;
end

%%

size_t = size(x_ode,1); % get the number of time steps of the simulation
z_ode = zeros(size_t,n_par); w_ode = zeros(size_t,n_par);
for i = 1:n_par % extact the state-of-charge and relaxation voltage from the simulation
    z_ode(:,i) = x_ode(:,2*i-1);
    w_ode(:,i) = x_ode(:,2*i);
end

nT = max(size(t_ode));i_branch2 = zeros(n_par,nT);
phi = zeros(n_par,1); phi_init = zeros(n_par-1,1);

for j = 1:nT
    for i = 1:n_par
        OCV_ode(j,i) = ocv_eval(ocv_coefs,z_ode(j,i));
        volts(j,i) = OCV_ode(j,i)+ w_ode(j,i);
    end

    for i = 2:n_par
        phi_init(i-1) = volts(j,i)- volts(j,i-1);
    end
    phi =[I;phi_init];

    i_branch2(:,j) = new_A22\phi; % Get the current going into each parallel branch
    current_sum_2(j) = sum(i_branch2(:,j))/I;
    error_current_time(j) = norm(i_branch2(:,j)-i_branch(:,j));
end

%% Plot the results
close all;
f_size = 16; f_size_gca = 13;
color_1 = 0.9*[1,1,1];
color1 = '#4169E1';  % Royal Blue for the curve
color2 = '#DC143C';  % Crimson Red for the data points
fig1 = figure;
plot(z_space, OCV_plot, 'Color', color1, 'LineWidth', 2);  % Royal Blue for the OCV curve
hold on
plot(ocv_points(:, 1), ocv_points(:, 2), 'o', 'Color', color2, 'MarkerSize', 8, 'LineWidth', 2);  % Crimson Red for data points
xlabel('State-of-charge ($z(t)$)', 'interpreter', 'latex', 'fontsize', f_size, 'interpreter', 'latex');
ylabel('Open-circuit-voltage (OCV) [V]', 'interpreter', 'latex', 'fontsize', f_size, 'interpreter', 'latex');
legend('OCV Curve', 'OCV Data Points', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', f_size);
% title('K2 LFP26650P', 'interpreter', 'latex', 'fontsize', f_size);
grid on;
g = gca;
set(g, 'fontsize', f_size_gca);


fig2 = figure; 
hold on;
% Define colormap (e.g., 'jet' for distinct colors)
cmap = jet(n_par);  % Creates a colormap with 'n_par' distinct colors
for j = 1:n_par
    plot(t_ode2/3600, z_ode(:,j), '-k', 'LineWidth', 1.5, 'Color', cmap(j,:));  % Assign different colors
end
xlabel('Time (hrs)', 'interpreter', 'latex', 'fontsize', f_size, 'interpreter', 'latex')
ylabel('State-of-charge ($z(t)$)', 'interpreter', 'latex', 'fontsize', f_size, 'interpreter', 'latex')
grid on
cell_labels = arrayfun(@(x) sprintf('Cell %d', x), 1:n_par, 'UniformOutput', false);  % Generate cell labels
legend(cell_labels, 'interpreter', 'latex', 'fontsize', 14, 'interpreter', 'latex')
g = gca;
set(g, 'fontsize', f_size_gca); 
box;

fig3 = figure; hold on;
for j = 1:n_par
    plot(t_ode2/3600, w_ode(:,j),'-.k','LineWidth',1.5,'Color', cmap(j,:));
end
xlabel('Time (hrs)','interpreter','latex','fontsize',f_size,'interpreter','latex')
ylabel('Relaxation voltage $v_{r}(t)$','interpreter','latex','fontsize',f_size,'interpreter','latex')
grid on
legend(cell_labels, 'interpreter', 'latex', 'fontsize', 14, 'interpreter', 'latex')
g = gca;
set(g, 'fontsize',f_size_gca); box;


fig4 = figure;
hold on;
for j = 1:n_par
    plot(t_ode2/3600, i_branch2(j,:), '-k', 'LineWidth', 1.5, 'Color', cmap(j,:));
end
xlabel('Time (hrs)', 'interpreter', 'latex', 'fontsize', f_size, 'interpreter', 'latex')
ylabel('Branch currents $i(t)$', 'interpreter', 'latex', 'fontsize', f_size, 'interpreter', 'latex')
grid on
legend(cell_labels, 'interpreter', 'latex', 'fontsize', 14, 'interpreter', 'latex')
title('Inverting Numerically', 'interpreter', 'latex', 'fontsize', f_size);
g = gca;
set(g, 'fontsize', f_size_gca);
box;


fig5 = figure;  
hold on;
for j = 1:n_par
    plot(t_ode2/3600, i_branch(j,:), '-k', 'LineWidth', 1.5, 'Color', cmap(j,:));
end
xlabel('Time (hrs)', 'interpreter', 'latex', 'fontsize', f_size, 'interpreter', 'latex')
ylabel('Branch currents $i(t) [A]$', 'interpreter', 'latex', 'fontsize', f_size, 'interpreter', 'latex')
grid on
legend(cell_labels, 'interpreter', 'latex', 'fontsize', 14, 'interpreter', 'latex')
title('Inverting by Analytical Expression', 'interpreter', 'latex', 'fontsize', f_size);
g = gca;
set(g, 'fontsize', f_size_gca);
box;


% fig5 = figure; hold on;
% for j = 1:n_par
%     plot(t_ode2/3600, current_sum,'-.k','LineWidth',1.5,'color',color_1*0/n_par)
% end
% xlabel('Time (hrs)','interpreter','latex','fontsize',f_size,'interpreter','latex')
% ylabel('$\frac{1}{I(t)}\sum_{k = 1}^n i_k(t)$','interpreter','latex','fontsize',f_size,'interpreter','latex')
% grid on
% % leg = legend('Cell 1','Cell 2','location','northwest');
% % set(leg,'interpreter','latex','fontsize',f_size)
% g = gca;
% set(g, 'fontsize',f_size_gca); box;
% axis([0 max(t_ode2/3600) 0 2]);
% 
% fig6 = figure; hold on;
% for j = 1:n_par
%     plot(t_ode/3600, current_sum_2,'-.k','LineWidth',1.5,'color',color_1*0/n_par)
% end
% xlabel('Time (hrs)','interpreter','latex','fontsize',f_size,'interpreter','latex')
% ylabel('$\frac{1}{I(t)}\sum_{k = 1}^n i_k(t)$','interpreter','latex','fontsize',f_size,'interpreter','latex')
% grid on
% % leg = legend('Cell 1','Cell 2','location','northwest');
% % set(leg,'interpreter','latex','fontsize',f_size)
% g = gca;
% set(g, 'fontsize',f_size_gca); box;
% axis([0 max(t_ode/3600) 0 2]);
% 
% fig7 = figure; hold on;
% for j = 1:n_par
%     plot(t_ode/3600, error_current_time,'-.k','LineWidth',1.5,'color',color_1*0/n_par)
% end
% xlabel('Time (hrs)','interpreter','latex','fontsize',f_size,'interpreter','latex')
% ylabel('error in the current','interpreter','latex','fontsize',f_size,'interpreter','latex')
% grid on
% % leg = legend('Cell 1','Cell 2','location','northwest');
% % set(leg,'interpreter','latex','fontsize',f_size)
% g = gca;
% set(g, 'fontsize',f_size_gca); box;
% % axis([0 max(t_ode/3600) 0 2]);

%%
% print(fig1,'ocv','-depsc')
% print(fig2,'soc','-depsc')
% print(fig3,'w_relax','-depsc')
% print(fig4,'i_branch','-depsc')
%% A function to evaluate the vector field of the parallel pack model
function xdot = parrallel_model(t,x,I,A11,A12,A21,A22,ocv_coefs)
n = max(size(x)/2);
Soc = zeros(n,1);  OCV_branch = zeros(n,1);  w = zeros(n); v_mod = zeros(n,1);

for i = 1:n
    Soc(i) = x(2*i-1)'; OCV_branch(i) = ocv_eval(ocv_coefs,Soc(i)); w(i) = x(2*i)'; % extract the soc and v_1
    v_mod(i) =  OCV_branch(i) + 1*w(i);
end

delta_v = zeros(n-1,1);
for i = 2:n
    delta_v(i-1) = v_mod(i)-v_mod(1); % compute the voltage difference with respect to cell one
end

phi = [-delta_v;-I];
i_branch = -A22\phi; % Compute the branch currents.
xdot = A11*x+A12*(-A22\phi); % the vector field
end

%% A function to evaluate the vector field of the parallel pack model
function xdot = parrallel_model_no_mat_inverse(t,x,I,A11,A12,A21,A22,ocv_coefs,m)
n = max(size(x)/2);
Soc = zeros(n,1);  OCV_branch = zeros(n,1);  w = zeros(n); v_mod = zeros(n,1);

for i = 1:n
    Soc(i) = x(2*i-1)'; OCV_branch(i) = ocv_eval(ocv_coefs,Soc(i)); w(i) = x(2*i)'; % extract the soc and v_1
    v_mod(i) =  OCV_branch(i) + 1*w(i);
end

delta_v = zeros(n-1,1);
for i = 2:n
    delta_v(i-1) = v_mod(i)-v_mod(1); % compute the voltage difference with respect to cell one
end

phi = [-delta_v;-I];
i_branch = -A22\phi; % Compute the branch currents.
xdot = A11*x+A12*(-m*phi); % the vector field. THIS IS THE MAIN DIFFERERNCE WITH RESPECT TO THE FUNCTION ABOVE. It uses m = inv(A) with m analytic in terms of the parameters.
end

%% A function which evaluates the polynomial vecotr field of the open circuit voltage
function ocv = ocv_eval(ocv_coefs,z)
ocv = polyval(ocv_coefs,z);
end


%%
function xdot = parallel_model_dae5(t, x, I, A11, A12, A21, A22, ocv_coefs, R, r)
global i_branch_store_2;
n = length(x) / 2;
Soc = zeros(n,1);  OCV_branch = zeros(n,1);  w = zeros(n,1); v_mod = zeros(n,1);

for i = 1:n
    Soc(i) = x(2*i-1)'; OCV_branch(i) = ocv_eval(ocv_coefs,Soc(i)); w(i) = x(2*i)'; % extract the soc and v_1
    v_mod(i) =  OCV_branch(i) + w(i);
end

i_branch = zeros(n, 1); theta = zeros(n, 1); rho_branch = zeros(n, 1); omega_branch = zeros(n, 1); alpha = zeros(n, 1); phi = zeros(n-1, 1);

for i = 2:n
    % theta(i) = (r(i))/r(i-1);
    % omega(i) =  (R(i))/r(i-1);
    % rho(i) = (v_mod(i)-v_mod(i-1))/r(i-1);
    phi(i-1) =(v_mod(i)-v_mod(i-1));
    % alpha(i) = 1+theta(i)+omega(i);
end
phi_mod  = [I; phi];
i_branch = A22\phi_mod;
xdot = A11 * x + A12 * i_branch;
end

%%
function xdot = parallel_model_dae6(t, x, I, A11, A12, A21, A22, ocv_coefs, R, r)
global i_branch_store;
n = length(x) / 2;
Soc = zeros(n,1);  OCV_branch = zeros(n,1);  w = zeros(n,1); v_mod = zeros(n,1);

for i = 1:n
    Soc(i) = x(2*i-1)'; OCV_branch(i) = ocv_eval(ocv_coefs,Soc(i)); w(i) = x(2*i)'; % extract the soc and v_1
    v_mod(i) =  OCV_branch(i) + w(i);
end

i_branch = zeros(n, 1); theta = zeros(n, 1); rho_branch = zeros(n, 1); omega_branch = zeros(n, 1); alpha = zeros(n, 1); phi = zeros(n-1, 1);

for i = 2:n
    theta(i) = (r(i))/r(i-1);
    omega(i) =  (R(i))/r(i-1);
    rho(i) = (v_mod(i)-v_mod(i-1))/r(i-1);
    phi(i-1) =(v_mod(i)-v_mod(i-1));
    alpha(i) = 1+theta(i)+omega(i);
end

beta = zeros(n+2, 1);
beta(n+1) = 1;
beta(n) = alpha(n);
for j = n-1:-1:2
    beta(j) = alpha(j)*beta(j+1)-theta(j)*beta(j+2);
end

f = zeros(n+1, 1);
f(n) = rho(n); f(n+1) = 0;
for j = n-1:-1:2
    f(j) = alpha(j)*f(j+1)-theta(j)*f(j+2)+rho(j);
end

%%

i_branch(n) = (I-f(2))/beta(2);

for j = n-1:-1:1
    i_branch(j) = theta(j+1)* i_branch(j+1)+omega(j+1)*sum(i_branch)+rho(j+1);
end

i_branch_store = i_branch;
%%
phi_mod  = [I; phi];
i_branch_invert = A22\phi_mod;
%%
check = sum(i_branch)/I;
check_invert = sum(i_branch_invert)/I;
%%

xdot = A11 * x + A12 * i_branch;
end

