clear; close all;

%% set-up
Nn = 2;
Nt = 5;
Nr = 1; 
kappa = 0.8; % share of missings
y = randn(Nn*Nt, 1);
yobs = y;
yobs(3) = NaN;
yobs(6) = NaN;
yobs(7:8) = NaN; 
yobs(10) = NaN;
Yobs = [yobs(1:2:end), yobs(2:2:end)]'; % matrix form
eta = randn(Nt, Nr);
z = [eta; yobs];

%%  params
if Nr == 1  % one factor
    pphi = 0.8;
    Q0 = eye(Nr);
    oomeg_eta = eye(Nr);
else
    error("Nr has to be equal to 1")
end

% params
llambda = unifrnd(-1,1,Nn,Nr);
rrho = zeros(Nn,1);
oomeg_eps = diag(unifrnd(1,2,Nn,1));

%% sparse system matrices
Llambda = kron(speye(Nt), llambda);
iOomeg_eps = kron(speye(Nt), eye(Nn) / oomeg_eps); 
Rrho = speye(Nt * Nn) + kron(spdiags(ones(Nt,1), -1, Nt, Nt), -diag(rrho));
iOomeg_eta = kron(speye(Nt), eye(Nr) / oomeg_eta);
iOomeg_eta(1:Nr, 1:Nr) = Q0;
Pphi = speye(Nt * Nr) + kron(spdiags(ones(Nt,1), -1, Nt, Nt), -pphi);

% Q
Q_y = Rrho' * iOomeg_eps * Rrho;
Q_eta = Pphi' * iOomeg_eta * Pphi + Llambda' * Q_y * Llambda; 
Q_eta_y = -Llambda' * Q_y; 
Q = [Q_eta, Q_eta_y; Q_eta_y', Q_y];

%% 3-block permutation
p_3b = [[1:Nr*Nt],[3, 6, 7, 8, 10, 1, 2, 4, 5, 9]+Nr*Nt];
z_3b = z(p_3b);
disp(z)
disp(z_3b)

%% time-t permutation
r_fac = [];
r_y = [];
Nmis = sum(sum(isnan(Yobs)));
counter_t = 0;
counter_xobs = Nmis + Nt * Nr; 
for t=1:Nt
    % factors
    r_fac = [r_fac; counter_t + (1:Nr)'];
    counter_t = counter_t + Nr; 
    
    % missing obs
    for i = 1:Nn
        if isnan(Yobs(i, t))
            counter_t = counter_t + 1;
            r_y = [r_y; counter_t];
            
        else
            counter_xobs = counter_xobs + 1;
            r_y = [r_y; counter_xobs];            
        end
    end
end  

r = [r_fac; r_y];
p_tau(r) = 1:length(r); 
z_tau = z(p_tau);
disp(z)
disp(z_tau)

%% calc Q_3b and Q_tau

Q_3b = Q(p_3b, p_3b);
Q_tau = Q(p_tau, p_tau);

%% plot Q_3b and Q_tau

xy_box = Nr*Nt+Nmis+0.5;
tick_vals = 0:5:size(Q,1);
tick_vals(1) = 1;
xlab_3b = ['bandwidth of highlighted section: ' num2str(bandwidth(Q_3b(1:Nr*Nt+Nmis, 1:Nr*Nt+Nmis)))];
xlab_tau = ['bandwidth of highlighted section: ' num2str(bandwidth(Q_tau(1:Nr*Nt+Nmis, 1:Nr*Nt+Nmis)))];
figure(1);
subplot(1, 2, 1);
spy(Q_3b)
hold on
plot([0, xy_box], [xy_box, xy_box], 'k-', 'LineWidth',2)
plot([xy_box, xy_box], [0, xy_box], 'k-', 'LineWidth',2)
hold off
xticks(tick_vals)
yticks(tick_vals)
xlabel(xlab_3b, 'FontSize', 8)
%title('\bf{Q}_{3b}','interpreter','tex','FontSize',12, 'FontName','Helvetica')

subplot(1, 2, 2)
spy(Q_tau)
hold on
plot([0, xy_box], [xy_box, xy_box], 'k-', 'LineWidth',2)
plot([xy_box, xy_box], [0, xy_box], 'k-', 'LineWidth',2)
hold off
xticks(tick_vals)
yticks(tick_vals)
xlabel(xlab_tau, 'FontSize', 8)
%title('\bf{Q}_{\tau}','interpreter','tex','FontSize',12, 'FontName','Helvetica')
exportgraphics(gcf, 'fig_Q3b_Qtau.png', 'Resolution', 300); 
