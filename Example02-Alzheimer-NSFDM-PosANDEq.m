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

T          = 2;
h          = 0.4;
t          = (0:h:T)';

M_null      = 10;
D_null      = 0.05;
O_null      = 0.02;
P_null      = 0.01;

% Step 3: Initialization of Solution Vectors

M_nsfdm    = zeros(length(t),1);
D_nsfdm    = zeros(length(t),1);
O_nsfdm    = zeros(length(t),1);
P_nsfdm    = zeros(length(t),1);

M_nsfdm(1) = M_null;
D_nsfdm(1) = D_null;
O_nsfdm(1) = O_null;
P_nsfdm(1) = P_null;

M_ee       = zeros(length(t),1);
D_ee       = zeros(length(t),1);
O_ee       = zeros(length(t),1);
P_ee       = zeros(length(t),1);

M_ee(1)    = M_null;
D_ee(1)    = D_null;
O_ee(1)    = O_null;
P_ee(1)    = P_null;

% Step 4: Loop over Time Grid

for j = 1:1:(length(t)-1)
  h            = t(j+1) - t(j);

  M_nsfdm(j+1) = (M_nsfdm(j)+h*(S+delta*M_nsfdm(j)))/(1+h*((delta/gamma)*M_nsfdm(j)+2*K_1*M_nsfdm(j)+K_2*D_nsfdm(j)+K_O*O_nsfdm(j)+mu_1));
  D_nsfdm(j+1) = (D_nsfdm(j)+h*(K_1*M_nsfdm(j)*M_nsfdm(j+1)))/(1+h*(K_2*M_nsfdm(j+1)+mu_2));
  O_nsfdm(j+1) = (O_nsfdm(j)+h*(K_2*M_nsfdm(j+1)*D_nsfdm(j+1)))/(1+h*(K_O*M_nsfdm(j+1)+mu_O));
  P_nsfdm(j+1) = (P_nsfdm(j)+h*(K_O*M_nsfdm(j+1)*O_nsfdm(j+1)))/(1+h*mu_P);

  M_ee(j+1)    = M_ee(j)+h*(S+delta*M_ee(j)-(delta/gamma)*M_ee(j)*M_ee(j)-2*K_1*M_ee(j)*M_ee(j)-K_2*D_ee(j)*M_ee(j)-K_O*O_ee(j)*M_ee(j)-mu_1*M_ee(j));
  D_ee(j+1)    = D_ee(j)+h*(K_1*M_ee(j)*M_ee(j)-K_2*M_ee(j)*D_ee(j)-mu_2*D_ee(j));
  O_ee(j+1)    = O_ee(j)+h*(K_2*M_ee(j)*D_ee(j)-K_O*M_ee(j)*O_ee(j)-mu_O*O_ee(j));
  P_ee(j+1)    = P_ee(j)+h*(K_O*M_ee(j)*O_ee(j)-mu_P*P_ee(j));
endfor

% Step 5: Plotting

figure(2)

hold on

subplot(2,2,1);
plot(t,M_nsfdm,'marker','*','linewidth',0.8);
title('Amount of M(t) (NSFDM, h = 0.4)','fontsize',14);
xlabel('t','fontsize',12);
ylabel('M(t)','fontsize',12);

subplot(2,2,2);
plot(t,D_nsfdm,'marker','*','linewidth',0.8);
title('Amount of D(t) (NSFDM, h = 0.4)','fontsize',14);
xlabel('t','fontsize',12);
ylabel('D(t)','fontsize',12);

subplot(2,2,3);
plot(t,O_nsfdm,'marker','*','linewidth',0.8);
title('Amount of O(t) (NSFDM, h = 0.4)','fontsize',14);
xlabel('t','fontsize',12);
ylabel('O(t)','fontsize',12);

subplot(2,2,4);
plot(t,P_nsfdm,'marker','*','linewidth',0.8);
title('Amount of P(t) (NSFDM, h = 0.4)','fontsize',14);
xlabel('t','fontsize',12);
ylabel('P(t)','fontsize',12);

hold off

figure(3)

hold on

subplot(2,2,1);
plot(t,M_ee,'marker','*','linewidth',0.8);
title('Amount of M(t) (EE, h = 0.4)','fontsize',14);
xlabel('t','fontsize',12);
ylabel('M(t)','fontsize',12);

subplot(2,2,2);
plot(t,D_ee,'marker','*','linewidth',0.8);
title('Amount of D(t) (EE, h = 0.4)','fontsize',14);
xlabel('t','fontsize',12);
ylabel('D(t)','fontsize',12);

subplot(2,2,3);
plot(t,O_ee,'marker','*','linewidth',0.8);
title('Amount of O(t) (EE, h = 0.4)','fontsize',14);
xlabel('t','fontsize',12);
ylabel('O(t)','fontsize',12);

subplot(2,2,4);
plot(t,P_ee,'marker','*','linewidth',0.8);
title('Amount of P(t) (EE, h = 0.4)','fontsize',14);
xlabel('t','fontsize',12);
ylabel('P(t)','fontsize',12);

hold off

% Step 6: Eigenvalues

J  = [delta-2*(delta/gamma)*M_nsfdm(end)-4*K_1*M_nsfdm(end)-K_2*D_nsfdm(end)-K_O*O_nsfdm(end)-mu_1,-K_2*M_nsfdm(end),-K_O*O_nsfdm(end),0;2*K_1*M_nsfdm(end)-K_2*D_nsfdm(end),-K_2*M_nsfdm(end)-mu_2,0,0;K_2*D_nsfdm(end)-K_O*O_nsfdm(end),K_2*M_nsfdm(end),-K_O*M_nsfdm(end)-mu_O,0;K_O*O_nsfdm(end),0,K_O*M_nsfdm(end),-mu_P];
eig(J)
