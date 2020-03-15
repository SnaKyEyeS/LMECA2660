%% Plot the last computed solution
figure; hold on;
load uh.txt
load u.txt

for i = 1:length(uh)
    clf; hold on;
    plot(uh(i,:))
    plot(u(i,:))
    legend("Numerical solution", "Exact solution");
    pause(.01);
end


%% Plot results at different timesteps (of the last computed solution)
load uh.txt
load u.txt

figure; hold on;
n = length(u);
x = linspace(0,1,n);

% ct/L = 1/4
subplot(4,1,1); hold on;
plot(x, uh(n/4,:));
plot(x, u(n/4,:));
title("$\frac{ct}{L} = 1/4$", "interpreter", "latex");
xlabel("x/L", "interpreter", "latex")
ylabel("u/Q", "interpreter", "latex")
legend("Numerical solution", "Analytical solution", "location", "northwest")

% ct/L = 2/4
subplot(4,1,2); hold on;
plot(x, uh(n/2,:));
plot(x, u(n/2,:));
title("$\frac{ct}{L} = 1/2$", "interpreter", "latex");
xlabel("x/L", "interpreter", "latex")
ylabel("u/Q", "interpreter", "latex")

% ct/L = 3/4
subplot(4,1,3); hold on;
plot(x, uh(3*n/4,:));
plot(x, u(3*n/4,:));
title("$\frac{ct}{L} = 3/4$", "interpreter", "latex");
xlabel("x/L", "interpreter", "latex")
ylabel("u/Q", "interpreter", "latex")

% ct/L = 1
subplot(4,1,4); hold on;
plot(x, uh(n,:));
plot(x, u(n,:));
title("$\frac{ct}{L} = 1$", "interpreter", "latex");
xlabel("x/L", "interpreter", "latex")
ylabel("u/Q", "interpreter", "latex")


%% Global Diagnostics (of the last computed solution)
load global_diagnostics.txt

n = length(global_diagnostics);
x = linspace(0,1,n);


figure; hold on;
plot(x, global_diagnostics(1,:));
plot(x, global_diagnostics(2,:));
plot(x, global_diagnostics(3,:));
legend("Q", "Energy", "Error");

%% Plot the solutions
cd convection-diffusion

u = load('u_125.txt');
n = length(u); x = linspace(0,1,n);

uh_e2 = load('uh_125_e2.txt');
uh_e4 = load('uh_125_e4.txt');
uh_ds = load('uh_125_ds.txt');
uh_i4 = load('uh_125_i4.txt');

% ct/L = 1/4
subplot(4,1,1); hold on;
plot(x, u(n/4,:), 'linewidth', 1.5);
plot(x, uh_e2(n/4,:));
plot(x, uh_e4(n/4,:));
plot(x, uh_ds(n/4,:));
plot(x, uh_i4(n/4,:));
title("$\frac{ct}{L} = 1/4$", "interpreter", "latex");
xlabel("x/L", "interpreter", "latex")
ylabel("u/Q", "interpreter", "latex")
ylim([-5 15])

legend("Exact solution", "E2", "E4", "DS", "I4", "location", "northwest")

% ct/L = 2/4
subplot(4,1,2); hold on;
plot(x, u(n/2,:), 'linewidth', 1.5);
plot(x, uh_e2(n/2,:));
plot(x, uh_e4(n/2,:));
plot(x, uh_ds(n/2,:));
plot(x, uh_i4(n/2,:));
title("$\frac{ct}{L} = 1/2$", "interpreter", "latex");
xlabel("x/L", "interpreter", "latex")
ylabel("u/Q", "interpreter", "latex")
ylim([-5 15])

% ct/L = 3/4
subplot(4,1,3); hold on;
plot(x, u(3*n/4,:), 'linewidth', 1.5);
plot(x, uh_e2(3*n/4,:));
plot(x, uh_e4(3*n/4,:));
plot(x, uh_ds(3*n/4,:));
plot(x, uh_i4(3*n/4,:));
title("$\frac{ct}{L} = 3/4$", "interpreter", "latex");
xlabel("x/L", "interpreter", "latex")
ylabel("u/Q", "interpreter", "latex")
ylim([-5 15])

% ct/L = 1
subplot(4,1,4); hold on;
plot(x, u(n,:), 'linewidth', 1.5);
plot(x, uh_e2(n,:));
plot(x, uh_e4(n,:));
plot(x, uh_ds(n,:));
plot(x, uh_i4(n,:));
title("$\frac{ct}{L} = 1$", "interpreter", "latex");
xlabel("x/L", "interpreter", "latex")
ylabel("u/Q", "interpreter", "latex")
ylim([-5 15])


cd ..

%% Plot the global diagnostics
cd convection

e2 = load('glob_125_e2.txt');
e4 = load('glob_125_e4.txt');
ds = load('glob_125_ds.txt');
i4 = load('glob_125_i4.txt');
i6 = load('glob_125_i6.txt');

n = length(e2); x = linspace(0,1,n);

% Energy
subplot(1,2,1); hold on;
plot(x, e2(2,:));
plot(x, e4(2,:));
plot(x, ds(2,:));
plot(x, i4(2,:));
plot(x, i6(2,:));
title("Energy", "interpreter", "latex");
ylim([0 1])
xlabel("$\frac{ct}{L}$", "interpreter", "latex")
ylabel("$\frac{E_h}{E(0)}$", "interpreter", "latex")
legend("E2", "E4", "DS", "I4", "I6", 'location', 'southwest')

% Error
subplot(1,2,2); hold on;
plot(x, e2(3,:));
plot(x, e4(3,:));
plot(x, ds(3,:));
plot(x, i4(3,:));
plot(x, i6(3,:));
ylim([0 1])
xlabel("$\frac{ct}{L}$", "interpreter", "latex")
ylabel("$\frac{R_h}{\sqrt{E(0)}}$", "interpreter", "latex")
title("Error", "interpreter", "latex");


cd ..

%% Order of convergence
figure;

load E2.txt
load E4.txt
load DS.txt
load I4.txt
load I6.txt

loglog(E2(:,1), E2(:,2), '-o'); hold on;
loglog(E4(:,1), E4(:,2), '-o')
loglog(DS(:,1), DS(:,2), '-o')
loglog(I4(:,1), I4(:,2), '-o')
loglog(I6(:,1), I6(:,2), '-o')

grid
legend("E2", "E4", "DS", "I4"," I6", 'location', 'southeast')
xlabel("h/$\sigma_0$", 'interpreter', 'latex')
ylabel("Error", 'interpreter', 'latex')
title("Order of convergence of the different FD schemes", 'interpreter', 'latex')

