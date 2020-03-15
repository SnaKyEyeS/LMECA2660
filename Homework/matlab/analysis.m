%% Epsilon machine
e = 1;
e_old = 2*e;
while e_old ~= e
    e_old = e;
    e = e/2;
end
disp(e)

%% Fourier analysis
figure; hold on;

Q = 1;
L = 1;
a = [8 16 32];
b = [1/2 1/4 1/8];

for i = 1:3
   for j = 1:3
        sigma = L/a(i);
        h = sigma*b(j);
        n = L/h;
        X = linspace(-L/2, L/2-h, n);
        k = -n/2:n/2-1;
        k_discr = linspace(-n/2, n/2-1, 1000);
        
        subplot(3,3, 3*(i-1)+j); hold on;
        
        % Exact Fourier coefficients
        coeff = exp(-k_discr.^2 * sigma^2 * pi^2);
        plot(k_discr, log10(abs(coeff)));
        
        % Discrete Fourier coefficients
        u0 = Q/sqrt(pi*sigma^2) * exp(-(X/sigma).^2);
        dft = fftshift(fft(u0))/n;
        stem(k, log10(abs(dft)), '*', 'linestyle', 'none')
        
        % Plot styling
        xlim([k(1) k(end)])
        ylim([min([log10(abs(dft)) log10(abs(coeff))]) 0])
        title(sprintf("L/$\\sigma_0$ = %d, h/$\\sigma_0$ = %g", a(i), b(j)), 'interpreter', 'latex')
        ylabel("$\log_{10}\left|F(k)\right|$", "interpreter", "latex")
        xlabel("k", "interpreter", "latex")
   end
end



%% DS modal analysis

kh = linspace(0,pi,1000);

kh_star_e2 = sin(kh);
kh_star_e4 = (8*sin(kh) - sin(2*kh))/6;
kh_star_e6 = (45*sin(kh) - 9*sin(2*kh) + sin(3*kh))/30;
kh_star_ds = (exp(1i*kh)/3 + .5 - exp(-1i*kh) + exp(-2i*kh)/6) / (1i);

figure; hold on;
plot(kh/pi, kh/pi)
plot(kh/pi, kh_star_e2/pi)
plot(kh/pi, kh_star_e4/pi)
plot(kh/pi, kh_star_e6/pi)
% plot(kh/pi, real(kh_star_ds)/pi)
xlabel("$\frac{kh}{\pi}$", "interpreter", "latex", "fontsize", 18)
ylabel("Re$\left(\frac{k^*h}{\pi}\right)$", "interpreter", "latex", "fontsize", 18)
title("Modified wavenumber", "interpreter", "latex")

yyaxis right;
plot(kh(1:end)/pi, imag(kh_star_ds(1:end)/pi), "--")
ylabel("Im$\left(\frac{k^*h}{\pi}\right)$", "interpreter", "latex", "fontsize", 18)
yticks([-pi/8, 0])
yticklabels({"-\pi/8", "0"})
legend("Ideal line", "E2", "E4 & DS", "E6", "DS - Phase", "location", "northeast")


%% Stability analysis
figure; hold on;
CFL = 1;

kh = (0:.25:1)*pi;
f_e2 = -1i*sin(kh) * CFL;
f_e4 = -1i*(8*sin(kh) - sin(2*kh))/6 * CFL;
f_e6 = -1i*(45*sin(kh) - 9*sin(2*kh) + sin(3*kh))/30 * CFL;
f_ds = -(exp(1i*kh)/3 + .5 - exp(-1i*kh) + exp(-2i*kh)/6) * CFL;
%f_ds = ((4/3*sin(kh)-1/6*sin(2*kh))*-1i + (-.5+2/3*cos(kh)-1/6*cos(2*kh))) * CFL;

% RK4
[x, y] = meshgrid(-3:.05:1, 0:-.05:-3); 
z = x + 1i*y;
f = abs(1 + z + z.^2/2 + z.^3/6 + z.^4/24);
contour(x,y,-f, -1:0.1:0, 'Linewidth', .05); 
grid;

plot(real(f_e2), imag(f_e2), 'o', 'Linewidth', 2);
plot(real(f_e4), imag(f_e4), '^', 'Linewidth', 2);
plot(real(f_e6), imag(f_e6), 'd', 'Linewidth', 2);
plot(real(f_ds), imag(f_ds), '*', 'Linewidth', 1.5);

kh = linspace(0,pi,1000);
f_ds = -(exp(1i*kh)/3 + .5 - exp(-1i*kh) + exp(-2i*kh)/6) * CFL;
plot(real(f_ds), imag(f_ds), 'Linewidth', .1, 'Color', [0.4940 0.1840 0.5560]);

xlabel("Re($\lambda_k\Delta t$)", "interpreter", "latex", "fontsize", 14)
ylabel("Im($\lambda_k\Delta t$)", "interpreter", "latex", "fontsize", 14)
legend("RK4 stability region", "E2", "E4", "E6", "DS", "location", "southwest");
title("Stability analysis", "interpreter", "latex")



