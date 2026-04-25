% Example 5: Parameter Sensitivity - Varying Delta

% Step 1: Definition of All Problem Parameters

# Values from Lindstrom 2021 and Cindyawati 2024

S        = 3.6*10^(-12)
delta_1  = 0.1;
delta_2  = 1;
delta_3  = 100;
gamma    = 75;
K_1      = 10^(-4);
K_2      = 5*10^(-4);
K_O      = 0.1;
mu_1     = 0.001;
mu_2     = 0.001;
mu_O     = 0.001;
mu_P     = 0.001;

% Step 2: Definition of Time Grid and Initial Conditions

T          = 10000;
h          = 10;
t          = (0:h:T)';

M_null      = 10;
D_null      = 0;
O_null      = 0;
P_null      = 0;

% Step 3: Initialization of Solution Vectors

M_nsfdm_1    = zeros(length(t),1);
D_nsfdm_1    = zeros(length(t),1);
O_nsfdm_1    = zeros(length(t),1);
P_nsfdm_1    = zeros(length(t),1);

M_nsfdm_1(1) = M_null;
D_nsfdm_1(1) = D_null;
O_nsfdm_1(1) = O_null;
P_nsfdm_1(1) = P_null;

M_nsfdm_2    = zeros(length(t),1);
D_nsfdm_2    = zeros(length(t),1);
O_nsfdm_2    = zeros(length(t),1);
P_nsfdm_2    = zeros(length(t),1);

M_nsfdm_2(1) = M_null;
D_nsfdm_2(1) = D_null;
O_nsfdm_2(1) = O_null;
P_nsfdm_2(1) = P_null;

M_nsfdm_3    = zeros(length(t),1);
D_nsfdm_3    = zeros(length(t),1);
O_nsfdm_3    = zeros(length(t),1);
P_nsfdm_3    = zeros(length(t),1);

M_nsfdm_3(1) = M_null;
D_nsfdm_3(1) = D_null;
O_nsfdm_3(1) = O_null;
P_nsfdm_3(1) = P_null;

% Step 4: Loop over Time Grid

for j = 1:1:(length(t)-1)
  h            = t(j+1) - t(j);

  M_nsfdm_1(j+1) = (M_nsfdm_1(j)+h*(S+delta*M_nsfdm(j)))/(1+h*((delta/gamma)*M_nsfdm(j)+2*K_1*M_nsfdm(j)+K_2*D_nsfdm(j)+K_O*O_nsfdm(j)+mu_1));
  D_nsfdm_1(j+1) = (D_nsfdm_1(j)+h*(K_1*M_nsfdm(j)*M_nsfdm(j+1)))/(1+h*(K_2*M_nsfdm(j+1)+mu_2));
  O_nsfdm_1(j+1) = (O_nsfdm_1(j)+h*(K_2*M_nsfdm(j+1)*D_nsfdm(j+1)))/(1+h*(K_O*M_nsfdm(j+1)+mu_O));
  P_nsfdm_1(j+1) = (P_nsfdm_1(j)+h*(K_O*M_nsfdm(j+1)*O_nsfdm(j+1)))/(1+h*mu_P);
endfor

% Step 5: Plotting

figure(4)

hold on

subplot(2,2,1);
plot(t,M_nsfdm,'linewidth',0.8);
hold on
plot(t,M_nsfdm(end)*ones(length(t)),'linestyle','--','linewidth',0.8,'color','black');
title('Amount of M(t) (NSFDM, h = 10)','fontsize',14);
xlabel('t','fontsize',12);
ylabel('M(t)','fontsize',12);

subplot(2,2,2);
plot(t,D_nsfdm,'linewidth',0.8);
hold on
plot(t,D_nsfdm(end)*ones(length(t)),'linestyle','--','linewidth',0.8,'color','black');
title('Amount of D(t) (NSFDM, h = 10)','fontsize',14);
xlabel('t','fontsize',12);
ylabel('D(t)','fontsize',12);

subplot(2,2,3);
plot(t,O_nsfdm,'linewidth',0.8);
hold on
plot(t,O_nsfdm(end)*ones(length(t)),'linestyle','--','linewidth',0.8,'color','black');
title('Amount of O(t) (NSFDM, h = 10)','fontsize',14);
xlabel('t','fontsize',12);
ylabel('O(t)','fontsize',12);

subplot(2,2,4);
plot(t,P_nsfdm(end),'linewidth',0.8);
hold on
plot(t,547.14*ones(length(t)),'linestyle','--','linewidth',0.8,'color','black');
title('Amount of P(t) (NSFDM, h = 10)','fontsize',14);
xlabel('t','fontsize',12);
ylabel('P(t)','fontsize',12);

hold off

% Step 6: Eigenvalues

J  = [delta-2*(delta/gamma)*M_nsfdm(end)-4*K_1*M_nsfdm(end)-K_2*D_nsfdm(end)-K_O*O_nsfdm(end)-mu_1,-K_2*M_nsfdm(end),-K_O*O_nsfdm(end),0;2*K_1*M_nsfdm(end)-K_2*D_nsfdm(end),-K_2*M_nsfdm(end)-mu_2,0,0;K_2*D_nsfdm(end)-K_O*O_nsfdm(end),K_2*M_nsfdm(end),-K_O*M_nsfdm(end)-mu_O,0;K_O*O_nsfdm(end),0,K_O*M_nsfdm(end),-mu_P];
eig(J)
