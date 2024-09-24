close all;
clc;

% Parâmetros
m = 10;
k = 4000;
c = 20;
wn = sqrt(k / m);
zeta = c / (2*wn*m);
wd = wn * sqrt(1 - zeta^2);
s = 0.0005; % passo de tempo
tf = 5;
t = 0:s:tf; % vetor de tempo

% Função de força y(t)
y = @(t) 0.005 * sin(5 * t);

% Função de força y'(t)
dydt = @(t) 5 * 0.005 * cos(5 * t);

% Função de força f(t)
f = @(t) k*y(t) + c*dydt(t);

% Funão resposta ao impulso h(t)
h = @(t) exp(-zeta * wn * t) .* sin(wd * t) / (m * wd);

% Resposta utilizando convolução
x = 0:s:tf;
for i = 1: length(t)
    T = 0:s:t(i);
    x(i) = trapz(f(T) .* h(- T + t(i))) * s;
end

% Resposta utilizando Newmark
newmark = Newmark(m,c,k,0,0,s,f,0.25,0.5);
x_n = newmark.integrate_until(t(end));

% Resposta utilizando Diferencas Finitas
newmark = Newmark(m,c,k,0,0,s,f,0.0,0.5);
x_f = newmark.integrate_until(t(end));


% Gráfico da função f(t)
figure;
hold on;
grid on;

plot(t, f(t), 'DisplayName', 'f(t)', 'linewidth', 1.7);
xlabel("Tempo [segundos]");
ylabel("F(t) [metros]");
title("Entrada");

print(gcf, 'atividade2_3_entrada', '-dpng', '-r300');

% Gráfico da resposta x(t)
figure;
hold on;
grid on;

plot(t, x_n, 'DisplayName', 'Método de Newmark', 'linewidth', 1.7);
plot(t, x, 'DisplayName', 'Método da Convolução', 'linewidth', 1.7, 'linestyle','--');
plot(t, x_f, 'DisplayName', 'Método das Diferenças Dinitas', 'linewidth', 1.7, 'linestyle','-.');
xlabel("Tempo [segundos]");
ylabel("X(t) [metros]");
title("Resposta");

legend;

print(gcf, 'atividade2_3_resposta', '-dpng', '-r300');