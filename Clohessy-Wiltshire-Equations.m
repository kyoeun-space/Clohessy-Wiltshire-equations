%% diff
syms x(t) y(t) z(t) w t C C1 C2 C3 C4 C5 C6 C7

x_dot(t) = diff(x(t));
y_dot(t) = diff(y(t));
z_dot(t) = diff(z(t));

x_ddot(t) = diff(x(t), t, 2);
y_ddot(t) = diff(y(t), t, 2);
z_ddot(t) = diff(z(t), t, 2);

%% HCW Equation
fx = x_ddot(t) - 2*w*y_dot(t) - 3*(w^2)*x(t);
fy = y_ddot(t) + 2*w*x_dot(t);
fz = z_ddot(t) + (w^2)*z(t);

% integrate fy 
integral_fy = int(fy, t);
integral_fy_equation = integral_fy == C;

% delete y term
pre_x = 2*w*integral_fy + fx;
simplify_pre_x = simplify(pre_x);

%% ==== calculate x ==== %%
% x(t) = ~
x_equation = isolate(dsolve(simplify_pre_x == 2*w*C, t) == x(t), x(t));
% x(0) = ~
x0_equation = subs(x_equation, t, 0);
% x'(t) = ~
x_dot_equation = diff(x_equation, t);
% x'(0) = ~
x0_dot_equation = subs(x_dot_equation, t, 0);

%% ==== calculate y ==== %%
% y'(t) = ~~
cal_y_dot = integral_fy_equation - 2*w*x_equation;
y_dot_equation = simplify(expand(cal_y_dot));
% y'(0) = ~~
y0_dot_equation = subs(y_dot_equation, t, 0);
% y(t) = ~
y_equation = isolate(rewrite(dsolve(y_dot_equation, t), "exp") == y(t), y(t));
% y(0) = ~
y0_equation = subs(y_equation, t, 0);

%% ==== find constant C ==== %%
% C = ~
val_C = isolate(expand(y0_dot_equation + (2*w*x0_equation)), C);
C3pC4 = expand(-2/w * (x0_equation * w * (3/2) + y0_dot_equation));
C3nC4 = (-1/2i)*(y0_equation - C5);
% C5 = ~
val_C5 = isolate(expand(-(1/w)*(2*x0_dot_equation - w*y0_equation)), C5);
% C3 = ~
val_C3 = isolate(subs((1/2) * (C3pC4 + C3nC4), C5, rhs(val_C5)), C3);
% C4 = ~
val_C4 = isolate(subs((1/2) * (C3pC4 - C3nC4), C5, rhs(val_C5)), C4);

%% ==== calculate z ==== %%
% z(t) = ~
z_equation = isolate(subs(dsolve(fz==0) == z(t), [C3, C4], [C6, C7]), z(t));
% z'(t) = ~
z_dot_equation = diff(z_equation, t);
% z(0) = ~
z0_equation = subs(z_equation, t, 0);
% z'(0) = ~
z0_dot_equation = subs(z_dot_equation, t, 0);

%% ==== find constant C ==== %%
% C6 = ~
val_C6 = isolate((1/(2*w*1i)) * expand(z0_equation*w*1i - z0_dot_equation), C6);
% C7 = ~
val_C7 = isolate((1/(2*w*1i)) * expand(z0_equation*w*1i + z0_dot_equation), C7);

%% ==== substitute constants into euqations  ==== %%
% constants
coefficients = [C, C3, C4, C5, C6, C7];
values = [rhs(val_C), rhs(val_C3), rhs(val_C4), rhs(val_C5), rhs(val_C6), rhs(val_C7)];

% X
cal_x_equation = subs(x_equation, coefficients, values);
obtain_x_equation = expand(rewrite(cal_x_equation, "sincos"));
% Y
cal_y_equation = subs(y_equation, coefficients, values);
obtain_y_equation = expand(rewrite(cal_y_equation, "sincos"));
% Z
cal_z_equation = subs(z_equation, coefficients, values);
obtain_z_equation = expand(rewrite(cal_z_equation, "sincos"));
% X'
cal_x_dot_equation = subs(x_dot_equation, coefficients, values);
obtain_x_dot_equation = expand(rewrite(cal_x_dot_equation, "sincos"));
% Y'
cal_y_dot_equation = subs(y_dot_equation, coefficients, values);
obtain_y_dot_equation = expand(rewrite(cal_y_dot_equation, "sincos"));
% Z'
cal_z_dot_equation = subs(z_dot_equation, coefficients, values);
obtain_z_dot_equation = expand(rewrite(cal_z_dot_equation, "sincos"));

%% ==== substitute initial conditions into equations ==== %%
zero_coefficients = [w, x(0), y(0), z(0), subs(diff(x(t), t), t, 0), subs(diff(y(t), t), t, 0), subs(diff(z(t), t), t, 0)];
initial_values_1 = [0.001107, 500, 0, 0, 0, 0, 0]; % m
initial_values_2 = [0.001107, 0, 0, 0, 0, 0.1, 0]; % km/s or 1.4539e-08 rad/s || D*ω/2

x_equation_1 = subs(obtain_x_equation, zero_coefficients, initial_values_1);
y_equation_1 = subs(obtain_y_equation, zero_coefficients, initial_values_1);
z_equation_1 = subs(obtain_z_equation, zero_coefficients, initial_values_1);

x_equation_2 = subs(obtain_x_equation, zero_coefficients, initial_values_2);
y_equation_2 = subs(obtain_y_equation, zero_coefficients, initial_values_2);
z_equation_2 = subs(obtain_z_equation, zero_coefficients, initial_values_2);

%% ==== set the time ==== %%
t_values = linspace(0, 10000, 1000);

x_values_1 = subs(x_equation_1, t, t_values);
y_values_1 = subs(y_equation_1, t, t_values);
z_values_1 = subs(z_equation_1, t, t_values);

x_values_2 = subs(x_equation_2, t, t_values);
y_values_2 = subs(y_equation_2, t, t_values);
z_values_2 = subs(z_equation_2, t, t_values);

%% ==== plot ==== %%
figure;
subplot(2, 3, 1);
plot(t_values, rhs(x_values_1));
xlabel('Time [s]');
ylabel('x(t) [m]');
title('x(t)_1');

subplot(2, 3, 2);
plot(t_values, rhs(y_values_1));
xlabel('Time [s]');
ylabel('y(t) [m]');
title('y(t)_1');

subplot(2, 3, 3);
plot(t_values, rhs(z_values_1));
xlabel('Time [s]');
ylabel('z(t) [m]');
title('z(t)_1');

subplot(2, 3, 4);
plot(t_values, rhs(x_values_2));
xlabel('Time [s]');
ylabel('x(t) [m]');
title('x(t)_2');

subplot(2, 3, 5);
plot(t_values, rhs(y_values_2));
xlabel('Time [s]');
ylabel('y(t) [m]');
title('y(t)_2');

subplot(2, 3, 6);
plot(t_values, rhs(z_values_2));
xlabel('Time [s]');
ylabel('z(t) [m]');
title('z(t)_2');
