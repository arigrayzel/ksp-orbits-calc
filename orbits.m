options = odeset('RelTol',1e-6,'AbsTol',1e-6','MaxStep',1e-3);

au = 1.3599e9; %au in m
day = 6*60*60; %day in s

mu = 1.17233e18*(day)^2/(au)^3; %for the sun

m_earth = 5.975e24;
m_mars = 6.441e23;

a_earth = 1.49e11/au;
a_mars = 2.28e11/au;

v_earth = 29790/au*day;
v_mars = 24140/au*day;

t_0 = 0;
t_f = 700; %700 days

y_e0 = [a_earth;0;0;v_earth];
y_m0 = [a_mars;0;0;v_mars];

%[t_e,y_e] = ode45(@odefun, [t_0,t_f], y_e0, options);
%[t_m,y_m] = ode45(@odefun, [t_0,t_f], y_m0, options);
%
%save('t_e.mat', 't_e')
%save('y_e.mat', 'y_e')
%save('t_m.mat', 't_m')
%save('y_m.mat', 'y_m')

load('t_e.mat')
load('y_e.mat')
load('t_m.mat')
load('y_m.mat')

clf
figure(1) %orbits in space
hold on
plot(y_e(:,1), y_e(:,2), 'b')
plot(y_m(:,1), y_m(:,2), 'r')

figure(2) %x in time (for determining period numerically)
clf
hold on
plot(t_e, y_e(:,1), 'b')
plot(t_m, y_m(:,1), 'r')

T_earth_theoretical = 2*pi/sqrt(mu/a_earth^3)
T_mars_theoretical = 2*pi/sqrt(mu/a_mars^3)

%%%Part 2
a_transfer = (a_earth+a_mars)/2;
t_transfer = pi/sqrt(mu/a_transfer^3)

%To get to Mars
delta_v_pi = sqrt(mu/a_earth)*(sqrt(2*a_mars/(a_mars+a_earth))-1)*au/day
delta_v_alpha = sqrt(mu/a_mars)*(1-sqrt(2*a_earth/(a_mars+a_earth)))*au/day

delta_theta_transfer = 180
delta_theta_mars = t_transfer/T_mars_theoretical*360 %How much mars during a transfer orbit
delta_theta_earth = t_transfer/T_earth_theoretical*360 %how much earth moves during a transfer
delta_theta_launch =delta_theta_transfer- delta_theta_mars %angle change needed during launch

omega_earth = 360/T_earth_theoretical %degrees/day
omega_mars  = 360/T_mars_theoretical %degrees/day
omega_diff = omega_mars-omega_earth
transfer_period = -360/omega_diff/T_earth_theoretical


rel_angle_0 = 216-104;
departure_0 = -(rel_angle_0-delta_theta_launch)/omega_diff/T_earth_theoretical + 2020 %date of first departure


departure_dates = linspace(0,5,6)*transfer_period + ones(1,6)*departure_0
arrival_dates = departure_dates + t_transfer/T_earth_theoretical
figure(3)
clf
plot(departure_dates, arrival_dates, '.b')

delta_theta_return =delta_theta_earth- delta_theta_transfer

%angles at 2020.0:
thetas_0 = [104;216;216-104]; %Earth, Mars, diff
rates = [omega_earth;omega_mars;omega_diff]
t_arrival_days = (arrival_dates(1)-2020)*T_earth_theoretical
thetas_a = thetas_0+rates*t_arrival_days + [0;0;360]
%delta_theta_return - thetas_a(2)
t_wait = -(thetas_a(3)-delta_theta_return)/omega_diff

return_dates = arrival_dates + t_wait/T_earth_theoretical
return_arrival_dates = return_dates + t_transfer/T_earth_theoretical
figure(4)
clf
plot(return_dates, return_arrival_dates, '.b')
round_trip_time = return_arrival_dates(1)-departure_dates(1)


%To get back to Earth
delta_v_pi2 = sqrt(mu/a_mars)*(1-sqrt(2*a_earth/(a_mars+a_earth)))*au/day
delta_v_alpha2 = sqrt(mu/a_earth)*(sqrt(2*a_mars/(a_mars+a_earth))-1)*au/day



%%%%%%%%part 3: Determine ranges (280-80 degrees)
n_deltas = 10;
deltas = [linspace(0,60, n_deltas/2), linspace(300,360, n_deltas/2)];
t_0 = 0;
t_f = 360;

r_alphas = zeros(1,n_deltas);
times = zeros(1,n_deltas);
angles = zeros(1,n_deltas);
departures = zeros(1,n_deltas);

for ii = 1:n_deltas 
	delta = deltas(ii)
	y_e0 = [a_earth;0;0;v_earth] + 2/au*day*[0;0;-sind(delta)*delta_v_pi;cosd(delta)*delta_v_pi];
	[t_e,y_e] = ode45(@odefun, [t_0,t_f], y_e0, options);
	rs = vecnorm(y_e(:,1:2)')';
	r_alphas(ii) = max(rs);
	[M,I] = max(rs>a_mars);
	intersect_short = I(1);
	times(ii) = t_e(intersect_short);
	ang = atan2(y_e(intersect_short,2), y_e(intersect_short,1))*180/pi;
	if ang<0
		ang = ang + 360;
	end
	angles(ii) = ang;
	delta_theta_mars = t_e(intersect_short)/T_mars_theoretical*360; %how much mars moves during a transfer
	delta_theta_launch =ang- delta_theta_mars; %angle change needed during launch
	departure = -(rel_angle_0-delta_theta_launch)/omega_diff/T_earth_theoretical + 2020; %date of first departure
	departures(ii) = departure;

end
figure(5)
clf
plot(deltas, r_alphas, '.m')
yline(a_mars, '--r')

figure(6)
clf
plot(deltas, times, '.m')

figure(7)
clf
plot(deltas, angles, '.b')

figure(8)
clf
hold on
plot(departure_dates(1), t_transfer, '*b')
plot(departures, times, '*r')
%y=[x;y;x_dot;y_dot]
function dydt = odefun(t,y)
	mu = 1.3273e20*(24*60*60)^2/(1.49e11)^3; %for the sun
	r = y(1:2);
	dist = norm(r);
	%dydt(4:6) = -mu*r/dist;
	%dydt(1:3) = y(4:6);
	dydt = [y(3:4);(-mu/dist^3)*r];
end
