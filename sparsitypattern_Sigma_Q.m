clear; close all;

Nt = 20; 

phi_1 = 0.8; 
phi_2 = -0.2
sig2 = 1;

H = speye(Nt);

H = H + spdiags(-phi_1*ones(Nt, 1), -1, Nt, Nt);
H = H + spdiags(-phi_2*ones(Nt, 1), -2, Nt, Nt);

S = speye(Nt) * sig2; 

% covariance and precision matrix
Q = H' * S * H;
Sigma = speye(Nt) / Q;

% plot
figure;
fig = gcf;
fig.PaperOrientation = 'landscape';
subplot(1, 2, 1);
spy(Sigma, 5)
%title('$$\Sigma$$','interpreter','latex','FontSize',16)
title('\Sigma','interpreter','tex','FontSize',12)
xlabel('')
subplot(1, 2, 2);
spy(Q,5)
%title('$\mathbf{Q}$','interpreter','latex','FontSize',16)
title('\bf{Q}','interpreter','tex','FontSize',12, 'FontName','Helvetica')
%title('Q','interpreter','tex','FontSize',12, 'FontName','PazoMath')
xlabel('')
exportgraphics(gcf, 'fig_Sigma_Q.png', 'Resolution', 300); 

