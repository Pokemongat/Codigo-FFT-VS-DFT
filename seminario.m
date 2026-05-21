
clear; clc; close all;

%% paremetroc
fs  = 8000;           % frequência de amostragem (Hz)
T   = 0.5;            % duração (s)
N   = fs * T;         % número de amostras total
t   = (0:N-1) / fs;   % vetor tempo

%% frequencia em hz das notas musicais
f_Do  = 261.63;   % Dó
f_Mi  = 329.63;   % Mi
f_Sol =  392.00;   % Sol

%%  o sinal
x = sin(2*pi*f_Do*t) + sin(2*pi*f_Mi*t) + sin(2*pi*f_Sol*t);

%% DFT
function X = minha_DFT(x)
    N = length(x);
    X = zeros(1, N);
    for k = 0:N-1
        for n = 0:N-1
            X(k+1) = X(k+1) + x(n+1) * exp(-1j * 2 * pi * k * n / N);
        end
    end
end

%% --- DFT manual
N_dft = 4000;
x_sub = x(1:N_dft);
f_sub = (0:N_dft-1) * (fs / N_dft);

tic;
X_DFT = minha_DFT(x_sub);
t_dft = toc;

%%  FFT com zero-padding (potência de 2) 
N_fft = 2^nextpow2(N);          % próxima potência de 2 acima de N (4000 → 4096)
freqs = (0:N_fft-1) * (fs / N_fft);

tic;
X_FFT = fft(x, N_fft);         % MATLAB completa com zeros até N_fft
t_fft = toc;

%% Erro numérico (mesma janela)
X_FFT_sub = fft(x_sub, N_dft);   % mesma janela que a DFT para comparação justa
erro = max(abs(X_DFT - X_FFT_sub));

fprintf('Tempo DFT manual (N=%d):  %.4f s\n', N_dft, t_dft);
fprintf('Tempo FFT        (N=%d): %.6f s\n', N_fft, t_fft);
fprintf('Speedup:                   x%.1f\n', t_dft / t_fft);
fprintf('Erro maximo |DFT - FFT|:   %.2e\n', erro);
fprintf('Zero-padding:              %d → %d amostras\n', N, N_fft);

%% Figura principal
figure('Name', 'DFT vs FFT ', ...
       'NumberTitle', 'off', 'Position', [100, 80, 1000, 750]);

%  Sinal composto no tempo --
subplot(3, 1, 1);
plot(t(1:400), x(1:400), 'Color', '#FFFFFF', 'LineWidth', 1.4);
title('Sinal Composto no Tempo  -  Do + Mi + Sol', 'FontSize', 13);
xlabel('Tempo (s)');
ylabel('Amplitude');
xlim([t(1) t(400)]);
grid on;

%  Espectro pela DFT manual --
subplot(3, 1, 2);
stem(f_sub(1:N_dft/2), abs(X_DFT(1:N_dft/2)) / N_dft, ...
     'Color', '#D95319', 'LineWidth', 1.2, 'MarkerSize', 3);
hold on;
xline(f_Do,  '--', 'Do (262 Hz)',  'Color', '#0072BD', ...
      'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(f_Mi,  '--', 'Mi (330 Hz)',  'Color', '#77AC30', ...
      'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(f_Sol, '--', 'Sol (392 Hz)', 'Color', '#7E2F8E', ...
      'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
title(sprintf('Espectro - DFT Manual  (N = %d  |  Tempo: %.3f s)', N_dft, t_dft), ...
      'FontSize', 12);
xlabel('Frequencia (Hz)');
ylabel('|X(k)| / N');
xlim([0 600]);
grid on;

% Espectro pela FFT --
subplot(3, 1, 3);
plot(freqs(1:N_fft/2), abs(X_FFT(1:N_fft/2)) / N, ...
     'Color', '#0072BD', 'LineWidth', 1.5);
hold on;
xline(f_Do,  '--', 'Do (262 Hz)',  'Color', '#0072BD', ...
      'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(f_Mi,  '--', 'Mi (330 Hz)',  'Color', '#77AC30', ...
      'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(f_Sol, '--', 'Sol (392 Hz)', 'Color', '#7E2F8E', ...
      'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
title(sprintf('Espectro - FFT com zero-padding  (N = %d → %d  |  Tempo: %.6f s)', N, N_fft, t_fft), ...
      'FontSize', 12);
xlabel('Frequencia (Hz)');
ylabel('|X(k)| / N');
xlim([0 600]);
grid on;

sgtitle('Decomposicao Espectral do Acorde de Do Maior  -  DFT vs FFT', ...
        'FontSize', 15, 'FontWeight', 'bold');

%%  comparação de tempos  -
figure('Name', 'Comparacao de Tempos', 'NumberTitle', 'off', ...
       'Position', [150, 150, 750, 420]);

N_vals = [32, 64, 128, 256, 512, 1024, 2048, 4096]; 
t_dfts = zeros(size(N_vals));
t_ffts = zeros(size(N_vals));

for k = 1:length(N_vals)
    n_k = N_vals(k);
    x_k = sin(2*pi*f_Do*(0:n_k-1)/fs) + ...
          sin(2*pi*f_Mi*(0:n_k-1)/fs) + ...
          sin(2*pi*f_Sol*(0:n_k-1)/fs);

  %faz varias vez e devolve a mediana
    t_dfts(k) = timeit(@() minha_DFT(x_k));
    t_ffts(k) = timeit(@() fft(x_k, n_k));
end

semilogy(N_vals, t_dfts, 'o--r', 'LineWidth', 2, 'MarkerSize', 8, ...
         'DisplayName', 'DFT - O(N^2)');
hold on;
semilogy(N_vals, t_ffts, 's-b', 'LineWidth', 2, 'MarkerSize', 8, ...
         'DisplayName', 'FFT - O(N log N)');
xlabel('Numero de Amostras  N');
ylabel('Tempo de Execucao (s)  -  escala log');
title('Complexidade Computacional: DFT vs FFT', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 12);
grid on;
