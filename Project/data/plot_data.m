%% Plot mesh data
data = load("mesh.txt");

n_x = 360; %306
n_y = 553;

x = reshape(data(:,1), [n_x, n_y]);
y = reshape(data(:,2), [n_x, n_y]);
val1 = reshape(data(:,3), [n_x, n_y]);
val2 = reshape(data(:,4), [n_x, n_y]);

figure; hold on;
% polarplot3d(p);
h = surf(x,y,val2);
set(h, 'EdgeColor', 'none');
% mesh(x,y,p);
view(-37.5,30)

%%
r = sqrt(x.^2 + y.^2);
t = atan2(y,x);

u =  (1 - .5^2./r.^2).*cos(t);
v = -(1 + .5^2./r.^2).*sin(t);
p = -(u.^2 + v.^2)/2;


%% Test
n_x = 360;
n_y = 553;

data = load("mesh.txt");
x = reshape(data(:,1), [n_x, n_y]);
y = reshape(data(:,2), [n_x, n_y]);
val2 = reshape(data(:,4), [n_x, n_y]);
r0 = .5;
r1 = 50 + r0;

r = sqrt(x.^2 + y.^2);
theta = atan2(y,x);
a = pi*(r0-r)/(r1-r0);

phi = sin(a).^2.*cos(2*theta);

% f =  2*cos(2*a).*cos(2*theta)*pi^2./H^2 - 2*cos(a).*sin(a).*cos(2*theta)*pi./(H*r) - 4*sin(a).^2.*cos(2*theta)./r.^2;

figure; hold on;
surf(x,y,abs(val2-phi), 'EdgeColor', 'none')

view(-37.5,30)


%% Ghost interpolation, second order 
syms a b c d 
syms h 
syms v0 v1 v2 v3 v_wall

eq1 = v1 == a*(h/2)^3 + b*(h/2)^2 + c*(h/2) + d;
eq2 = v2 == a*(3*h/2)^3 + b*(3*h/2)^2 + c*(3*h/2) + d;
eq3 = v3 == a*(5*h/2)^3 + b*(5*h/2)^2 + c*(5*h/2) + d;
eq4 = v_wall == d;

S = solve([eq1 eq2 eq3 eq4], [a b c d]);

sol = simplify(S.a*(-h/2)^3 + S.b*(-h/2)^2 + S.c*(-h/2) + S.d)

%% Ghost interpolation, fourth order
syms a b c d e 
syms h 
syms v0 v1 v2 v3 v4 v5 v_wall

eq1 = v1 == a*(1*h/2)^4 + b*(1*h/2)^3 + c*(1*h/2)^2 + d*(1*h/2)^1 + e;
eq2 = v2 == a*(3*h/2)^4 + b*(3*h/2)^3 + c*(3*h/2)^2 + d*(3*h/2)^1 + e;
eq3 = v3 == a*(5*h/2)^4 + b*(5*h/2)^3 + c*(5*h/2)^2 + d*(5*h/2)^1 + e;
eq4 = v4 == a*(7*h/2)^4 + b*(7*h/2)^3 + c*(7*h/2)^2 + d*(7*h/2)^1 + e;
eq5 = v_wall == e;

S = solve([eq1 eq2 eq3 eq4 eq5], [a b c d e]);

sol = simplify(S.a*(-h/2)^4 + S.b*(-h/2)^3 + S.c*(-h/2)^2 + S.d*(-h/2)^1 + S.e)

%% Plot live data
figure;

fid = fopen('mesh_u.txt', 'rt');
data = textscan(fid, "{'n': %d, 'n1': %d, 'n2': %d, 'd1': %f, 'd2': %f}");
n = cell2mat(data(1));
n1 = cell2mat(data(2));
n2 = cell2mat(data(3));

x = reshape(cell2mat(textscan(fid, repmat('%f,', 1, n))), n2, n1);
y = reshape(cell2mat(textscan(fid, repmat('%f,', 1, n))), n2, n1);


while ~feof(fid)
    data = cell2mat(textscan(fid, repmat('%f,', 1, n+1)));
    t = data(1);
    data = reshape(data(2:end), n2, n1);
    
    h = pcolor(x,y,data);
    set(h, 'edgecolor', 'none');
    title(t)
    colorbar;
    drawnow;
%     pause;
end

%% Test zone
close all;

data = load('mesh.txt');
n1 = 553;
n2 = 360;

x = reshape(data(:,1), [n2, n1]);
y = reshape(data(:,2), [n2, n1]);
r = sqrt(x.^2 + y.^2);
theta = atan2(y,x);
val = reshape(data(:,4), [n2, n1]);

f = r.*sin(theta).^2/4 - r.*cos(theta).^2/2;

figure; hold on;
surf(x,y,val, 'EdgeColor', 'none')
view(-37.5,30)
title("Numerical result")

figure; hold on;
surf(x,y,f, 'EdgeColor', 'none')
view(-37.5,30)
title("Analytical result")

figure; hold on;
surf(x,y,abs(f-val), 'EdgeColor', 'none')
view(-37.5,30)
title("Error")








