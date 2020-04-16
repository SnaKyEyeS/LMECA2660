%% Plot mesh data
data = load("mesh.txt");

n_x = 360; %306
n_y = 553-1;

x = reshape(data(:,1), [n_x, n_y]);
y = reshape(data(:,2), [n_x, n_y]);
val1 = reshape(data(:,3), [n_x, n_y]);
val2 = reshape(data(:,4), [n_x, n_y]);

figure; hold on;
% polarplot3d(p);
surf(x,y,val2, 'EdgeColor', 'none')
% mesh(x,y,p);
view(-37.5,30)

%% Test
data = load("mesh.txt");
x = reshape(data(:,1), [n_x, n_y]);
y = reshape(data(:,2), [n_x, n_y]);
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
