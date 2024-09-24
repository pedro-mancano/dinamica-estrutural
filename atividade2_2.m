close all;
clc;

% Parâmetros
m = 200;
k = 2e5;
wn = sqrt(k / m);
F0 = 10; % força arbitrária
s = 0.0005; % passo de tempo
tf = 5;
t = 0:s:tf; % vetor de tempo

% Função de força f(t)
f = @(t) F0 * (heaviside(t) - heaviside(t - 3));

% Funão resposta ao impulso h(t)
h = @(t) sin(wn * t) / (m * wn);

% Resposta utilizando convolução
x = 0:s:tf;
for i = 1: length(t)
    T = 0:s:t(i);
    x(i) = trapz(f(T) .* h(- T + t(i))) * s;
end

% Fast Fourier Transform (FFT)
X = fft(f(t)) .* fft(h(t));

% Inverse Fast Fourier Transform (IFFT)
x_reconstructed = ifft(X) * s;

% Resposta analítica
x_a = @(t) (F0 ./ (m .* wn ^ 2)) * ((1-cos(wn .* t)) .* heaviside(t) - (1 - cos(wn .* (t-3))) .* heaviside(t-3));

% Gráfico da função f(t)
figure;
hold on;
grid on;

plot(t, f(t), 'DisplayName', 'f(t)', 'linewidth', 1.7);
xlabel("Tempo [segundos]");
ylabel("F(t) [metros]");
title("Entrada");
ylim([-5, F0+5]);

print(gcf, 'atividade2_2_entrada', '-dpng', '-r300');

% Gráfico da resposta x(t)
figure;
hold on;
grid on;

plot(t, x, 'DisplayName', 'Método da Convolução', 'linewidth', 1.7);
plot(t, x_a(t), 'DisplayName', 'Método analítico', 'linewidth', 1.7, 'linestyle','--');
plot(t, x_reconstructed, 'DisplayName', 'Método da Transformada Inversa de Fourier', 'linewidth', 1.7, 'linestyle','-.');
xlabel("Tempo [segundos]");
ylabel("X(t) [metros]");
title("Resposta");

legend;

print(gcf, 'atividade2_2_resposta', '-dpng', '-r300');