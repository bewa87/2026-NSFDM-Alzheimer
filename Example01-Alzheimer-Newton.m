% Step 1: Definition of All Problem Parameters

S          = 0.01;
delta      = 0.05;
gamma      = 0.20;
K_1        = 0.04;
K_2        = 0.02;
K_O        = 0.01;
mu_1       = 0.001;
mu_2       = 0.001;
mu_O       = 0.001;
mu_P       = 0.001;

% Step 2: Definition of the function g

g   = @(x) S + (delta - mu_1)*x - (delta/gamma + 2*K_1)*x^2 - ((K_1*K_2)/(mu_2+K_2*x))*x^3 - (K_O*K_1*K_2*x^4)/((mu_O + K_O*x)*(mu_2 + K_2*x));
g_x = @(x) (delta - mu_1) - 2*(delta/gamma + 2*K_1)*x - (2*K_1*K_2*K_2*x^3 + 3*mu_2*K_1*K_2*x^2)/((mu_2 + K_2*x)^2) - (K_O*K_1*K_2*x^3*(4*mu_O*mu_2 + 3*mu_O*K_2*x + 3*mu_2*K_O*x + 2*K_O*K_2*x^2))/((mu_O + K_O*x)^2*(mu_2 + K_2*x)^2);

% Step 3: Newton-Iteration

N_max  = 10;
M      = zeros(1,N_max+1);
M_zero = 10;
M(1)   = M_zero;

for j = 1:1:N_max
  M(j+1) = M(j) - g(M(j))/(g_x(M(j)));
endfor

% Step 4: Computation of errors |M(j) - M(end)|

ind    = 1:1:(length(M)-1);
err    = zeros(1,length(ind));

for j = 1:1:length(ind)
  err(j) = log(abs(M(j+1) - M(j)));
endfor

% Step 5: Calculation of equilibrium state and verification

M_star  = M(end);
D_star  = (K_1*(M_star)^2)/(mu_2 + K_2*M_star);
O_star  = (K_2*M_star*D_star)/(mu_O + K_O*M_star);
P_star  = (K_O*M_star*O_star)/(mu_P);
Eq_star = [M_star D_star O_star P_star];

f1      = @(M,D,O,P) S + (delta-mu_1)*M - (delta/gamma)*M^2 - 2*K_1*M^2 - K_2*M*D - K_O*M*O;
f2      = @(M,D,O,P) K_1*M^2 - K_2*M*D - mu_2*D;
f3      = @(M,D,O,P) K_2*M*D - K_O*M*O - mu_O*O;
f4      = @(M,D,O,P) K_O*M*O - mu_P*P;

% Step 6: Error plot

figure(1)
plot(ind,err,'linewidth',0.8,'marker','*','linestyle','-')
hold on
plot(ind,-2*ind,'linewidth',0.8,'marker','+','linestyle','-.')
title('Logarithmic Error Plot Of Newton-Iteration','fontsize',14,'fontweight','normal')
legend({'Logarithmic Error','Theoretical Convergence'},'location','northeast','fontsize',12)
xlabel('Index k Of Newton-Iterate','fontsize',12)
ylabel('Logarithmic Error','fontsize',12)
hold off
