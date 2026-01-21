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

% Step 2: Definition of Time Grid and Initial Conditions

T          = 5000;
h          = 1;
t          = (0:h:T)';

M_null     = 0.20;
D_null     = 0.03;
O_null     = 0.01;
P_null     = 0.00;

% Step 3: Initialization of Solution Vectors

M_nsfdm    = zeros(length(t),1);
D_nsfdm    = zeros(length(t),1);
O_nsfdm    = zeros(length(t),1);
P_nsfdm    = zeros(length(t),1);

M_nsfdm(1) = M_null;
D_nsfdm(1) = D_null;
O_nsfdm(1) = O_null;
P_nsfdm(1) = P_null;

% Step 4: Loop over Time Grid

for j = 1:1:(length(t)-1)
  h            = t(j+1) - t(j);
  M_nsfdm(j+1) = (M_nsfdm(j)+h*(S+delta*M_nsfdm(j)))/(1+h*((delta/gamma)*M_nsfdm(j)+2*K_1*M_nsfdm(j)+K_2*D_nsfdm(j)+K_O*O_nsfdm(j)+mu_1));
  D_nsfdm(j+1) = (D_nsfdm(j)+h*(K_1*M_nsfdm(j)*M_nsfdm(j+1)))/(1+h*(K_2*M_nsfdm(j+1)+mu_2));
  O_nsfdm(j+1) = (O_nsfdm(j)+h*(K_2*M_nsfdm(j+1)*D_nsfdm(j+1)))/(1+h*(K_O*M_nsfdm(j+1)+mu_O));
  P_nsfdm(j+1) = (P_nsfdm(j)+h*(K_O*M_nsfdm(j+1)*O_nsfdm(j+1)))/(1+h*mu_P);
endfor

% Step 5: Plotting

figure(1)
plot(t,M_nsfdm)
title('Amount of M(t)','fontsize',14)
xlabel('t','fontsize',12)
ylabel('M(t)','fontsize',12)

figure(2)
plot(t,D_nsfdm)
title('Amount of D(t)','fontsize',14)
xlabel('t','fontsize',12)
ylabel('D(t)','fontsize',12)

figure(3)
plot(t,O_nsfdm)
title('Amount of O(t)','fontsize',14)
xlabel('t','fontsize',12)
ylabel('O(t)','fontsize',12)

figure(4)
plot(t,P_nsfdm)
title('Amount of P(t)','fontsize',14)
xlabel('t','fontsize',12)
ylabel('P(t)','fontsize',12)

figure(5)
plot(t,1*(M_nsfdm-M_nsfdm(end)*ones(length(M_nsfdm),1)).^2+1*(D_nsfdm-D_nsfdm(end)*ones(length(D_nsfdm),1)).^2+1*(O_nsfdm-O_nsfdm(end)*ones(length(O_nsfdm),1)).^2+1*(P_nsfdm-P_nsfdm(end)*ones(length(P_nsfdm),1)).^2)

AA = 1*(M_nsfdm-M_nsfdm(end)*ones(length(M_nsfdm),1)).^2+1*(D_nsfdm-D_nsfdm(end)*ones(length(D_nsfdm),1)).^2+1*(O_nsfdm-O_nsfdm(end)*ones(length(O_nsfdm),1)).^2+1*(P_nsfdm-P_nsfdm(end)*ones(length(P_nsfdm),1)).^2;

% Step 6: Eigenvalues

J  = [delta-2*(delta/gamma)*M_nsfdm(end)-4*K_1*M_nsfdm(end)-K_2*D_nsfdm(end)-K_O*O_nsfdm(end)-mu_1,-K_2*M_nsfdm(end),-K_O*O_nsfdm(end),0;2*K_1*M_nsfdm(end)-K_2*D_nsfdm(end),-K_2*M_nsfdm(end)-mu_2,0,0;K_2*D_nsfdm(end)-K_O*O_nsfdm(end),K_2*M_nsfdm(end),-K_O*M_nsfdm(end)-mu_O,0;K_O*O_nsfdm(end),0,K_O*M_nsfdm(end),-mu_P];
eig(J)
