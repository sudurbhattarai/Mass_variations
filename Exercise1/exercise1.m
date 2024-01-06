clc
clear


%% task 1a)
phi1=9*pi/180;
h1=(0:50:1000);
% gamma0 and gamma
normalgravity_10=zeros(1,21);
normalgravity_1=zeros(1,21);

for i=1:1:21
    [normalgravity_10(i),normalgravity_1(i)]=Normalgravity(phi1,h1(i));
end

figure

% Plot the contributions of gamma0  
subplot(2,2,1);
scatter(h1,normalgravity_1,'r','*');
title('Scatter Plot');
xlabel("h/(m)")
ylabel("NormalGravity/(m/s-2)")
grid on;

subplot(2,2,2);
plot(h1,normalgravity_1,'y');
title('Line Plot');
xlabel('h/(m)')
ylabel('NormalGravity/(m/s-2)')
grid on;

% Plot the the linear term with h
subplot(2,2,3);
p_linear = polyfit(h1, normalgravity_1, 1);
y_linear_fit = polyval(p_linear, h1);
scatter(h1,normalgravity_1,'r','*');
hold on
plot(h1,y_linear_fit,'g')
hold off
title('Linear Fit');
xlabel('h/(m)')
ylabel('NormalGravity/(m/s-2)')
grid on;

% Plot the the quadratic term with h^2
subplot(2, 2, 4);
scatter(h1.*h1,normalgravity_1,'r','*');
hold on
p_quadratic = polyfit(h1.*h1, normalgravity_1, 2);
y_quadratic_fit = polyval(p_quadratic, h1.*h1);
plot(h1.*h1,y_quadratic_fit,'b')
hold off
title('Quadratic Fit');
xlabel('h2/(m2)')
ylabel('NormalGravity/(m/s-2)')
grid on;

% and plot the sum in one plot
sgtitle('Combined Plot');


%% task 1b)
h2=100;
phi20=(0:1:90);
phi2=phi20*pi/180;

normalgravity_20=zeros(1,91);
normalgravity_2=zeros(1,91);

for j=1:1:91
     [normalgravity_20(j),normalgravity_2(j)]=Normalgravity(phi2(j),h2);
end

% Plot the contributions of gamma0  
subplot(2,2,1);
scatter(phi2,normalgravity_2,'r','*');
title('Scatter Plot');
xlabel("degree/(rad)")
ylabel("NormalGravity/(m/s-2)")
grid on;

subplot(2,2,2);
plot(phi2,normalgravity_2,'y');
title('Line Plot');
xlabel("degree/(rad)")
ylabel("NormalGravity/(m/s-2)")
grid on;

% Plot the the linear term with h
subplot(2,2,3);
p_linear = polyfit(phi2, normalgravity_2, 1);
y_linear_fit = polyval(p_linear, phi2);
scatter(phi2,normalgravity_2,'r','*');
hold on
plot(phi2,y_linear_fit,'g')
hold off
title('Linear Fit');
xlabel("degree/(rad)")
ylabel('NormalGravity/(m/s-2)')
grid on;

% Plot the the quadratic term with h^2
subplot(2, 2, 4);
scatter(phi2.*phi2,normalgravity_2,'r','*');
hold on
p_quadratic = polyfit(phi2.*phi2, normalgravity_2, 2);
y_quadratic_fit = polyval(p_quadratic, phi2.*phi2);
plot(phi2.*phi2,y_quadratic_fit,'b')
hold off
title('Quadratic Fit');
xlabel("degree^2/(rad^2)")
ylabel('NormalGravity/(m/s-2)')
grid on;

% and plot the sum in one plot
sgtitle('Combined Plot');


%normal gravity on the ellipsoid specify in a series development(depending
%on the geographical latitude phi
function [gamma0,gamma]=Normalgravity(phi,h)
% normal gravity at the equator
gammaa=9.7803267715;
gamma0=gammaa*(1+0.0052790414*((sin(phi))^2)+0.0000232718*((sin(phi))^4)+0.0000001262*((sin(phi))^6)+0.0000000007*((sin(phi))^8));

%The normal gravity at the ellipsoidal height h
gamma=gamma0-(0.0030877-0.0000043*((sin(phi))^2))*h+0.00000072*h*h;
end



