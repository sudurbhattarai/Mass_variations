clc
clear


%% task 1a)
% Parameters
phi1=9*pi/180;   % latitude in radian
h1=(0:50:1000);  % height increment in meters

% Initialize arrays to store contributions
normalgravity_10=zeros(1,21); % Normal gravity
normalgravity_1=zeros(1,21);  % Total gravity
h_l_1=zeros(1,21);            % Linear term
h_q_1=zeros(1,21);            % Quadratic term

% Compute contributions for each i
for i=1:1:21
    [normalgravity_10(i),normalgravity_1(i),h_l_1(i),h_q_1(i)]=Normalgravity(phi1,h1(i));
end

figure

% Plot the contributions of gamma0  
subplot(2,2,1);
m=(0:50:1000);
plot(m,normalgravity_10,'r');
title('Normal Gravity');
xlabel("h")
ylabel("Gravity/(m/s^2)")
grid on;

% Plot the the linear term with h
subplot(2,2,2);
plot(m,h_l_1,'y');
title('Linear Term');
xlabel("h")
ylabel('Gravity/(m/s^2)')
grid on;

% Plot the the Quadratic term with h^2
subplot(2,2,3);
plot(m,h_q_1,'g');
title('Quadratic Term');
xlabel("h")
ylabel('Gravity/(m/s^2)')
grid on;

% Plot the sum of above
subplot(2, 2, 4);
% plot(m,normalgravity_10,'r');
% hold on
% plot(m,h_l_1,'y');
% hold on
% plot(m,h_q_1,'g');
% hold on
plot(m,normalgravity_1,'b');
title('Total Gravity');
xlabel("h")
ylabel("Gravity/(m/s^2)")
grid on;
% hold off

% Combine 4 plots in one and name the title
sgtitle('Normal Gravity Contributions along Elipsoidal Normal');


%% task 1b)
% Parameters
h2=100; % constant ellipsoidal height in meters
phi20=(0:1:90); %latitude in degree
phi2=phi20*pi/180; %latitude in radian

% Initialize arrays to store contributions
normalgravity_20=zeros(1,91); % Normal gravity
normalgravity_2=zeros(1,91);  % Total gravity
h_l_2=zeros(1,91);            % Linear term
h_q_2=zeros(1,91);            % Quadratic term

% Compute contributions for each j 
for j=1:1:91
     [normalgravity_20(j),normalgravity_2(j),h_l_2(j),h_q_2(j)]=Normalgravity(phi2(j),h2);
end

figure

% Plot the contributions of gamma0  
subplot(2,2,1);
n=(0:1:90);
plot(n,normalgravity_20,'r');
title('Normal Gravity');
xlabel("longitude (φ/°)")
ylabel("Gravity/(m/s^2)")
grid on;

% Plot the the linear term with h
subplot(2,2,2);
plot(n,h_l_2,'y');
title('Linear Term');
xlabel("longitude (φ/°)")
ylabel('Gravity/(m/s^2)')
grid on;

% Plot the the quadratic term with h^2
subplot(2,2,3);
plot(n,h_q_2,'g');
title('Quadratic Term');
xlabel("longitude (φ/°)")
ylabel('Gravity/(m/s^2)')
grid on;

% Plot the sum of above
subplot(2, 2, 4);
% plot(n,normalgravity_20,'r');
% hold on
% plot(n,h_l_2,'y');
% hold on
% plot(n,h_q_2,'g');
% hold on
plot(n,normalgravity_2,'b');
title('Total Gravity');
xlabel("longitude (φ/°)")
ylabel("Gravity/(m/s^2)")
grid on;
% hold off

% Combine 4 plots in one and name the title
sgtitle('Normal Gravity Contributions along the Meridian');


%normal gravity on the ellipsoid specify in a series development(depending
%on the geographical latitude phi
function [gamma0,gamma,h_linear,h_quadratic]=Normalgravity(phi,h)
% normal gravity at the equator
gammaa=9.7803267715;
gamma0=gammaa*(1+0.0052790414*((sin(phi))^2)+0.0000232718*((sin(phi))^4)+0.0000001262*((sin(phi))^6)+0.0000000007*((sin(phi))^8));

%The normal gravity at the ellipsoidal height h
h_linear=-(0.0030877-0.0000043*((sin(phi))^2))*h;
h_quadratic=0.00000072*h*h;

gamma=gamma0+h_linear+h_quadratic;
end



