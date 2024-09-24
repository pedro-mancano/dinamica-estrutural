close all;
clc;

% Parâmetros
m = 200;
k = 10^6;
zeta = 0.06;
wn = sqrt(k / m);
c = 2 * zeta * wn * m;
wd = wn * sqrt(1 - zeta^2);
s = 0.0005; % passo de tempo
tf = 0.5;
t = 0:s:tf; % vetor de tempo

% Função de força f(t)
f = @(t) 200 * (1 - cos(pi * t / 0.2).^2) .* (1 - heaviside(t - 0.2));

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

% Fast Fourier Transform (FFT)
X = fft(f(t)) .* fft(h(t));

% Inverse Fast Fourier Transform (IFFT)
x_reconstructed = ifft(X) * s;

% Gráfico da função f(t)
figure;
hold on;
grid on;

plot(t, f(t), 'DisplayName', 'f(t)', 'linewidth', 1.7);
xlabel("Tempo [segundos]");
ylabel("F(t) [metros]");
title("Entrada");

print(gcf, 'atividade2_1_entrada', '-dpng', '-r300');

% Gráfico da resposta x(t)
figure;
hold on;
grid on;

plot(t, x_n, 'DisplayName', 'Método de Newmark', 'linewidth', 1.7);
plot(t, x_f, 'DisplayName', 'Método das Diferenças Finitas', 'linewidth', 1.7, 'linestyle','-.');
plot(t, x, 'DisplayName', 'Método da Convolução', 'linewidth', 1.7, 'linestyle','--');
plot(t, x_reconstructed, 'DisplayName', 'Método da Transformada Inversa de Fourier', 'linewidth', 1.7, 'linestyle','--');
xlabel("Tempo [segundos]");
ylabel("X(t) [metros]");
title("Resposta");
xline(0.2,'DisplayName', 'Fim da forçante', 'linewidth', 1.7, 'color', 'red');

legend;

print(gcf, 'atividade2_1_resposta', '-dpng', '-r300');


figure;
hold on;
grid on;

plot((1/s)/tf*(t), abs(X), 'linewidth', 1.7);
xlim([0, 50]);
xlabel('Frequencia (Hz)');
ylabel('X(w)');
title('Domínio da Frequencia');

print(gcf, 'atividade2_1_dominio_frequencia', '-dpng', '-r300');