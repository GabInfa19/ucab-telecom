% Practica Unidad I - Estimacion espectral clasica (Periodograma y Periodograma promedio)
% Compatible con MATLAB y Octave. No requiere toolboxes (las ventanas se definen a mano).
clear; close all; clc;
rand('seed',0); randn('seed',0);          % resultados repetibles
NFFT = 1024;                              % puntos de la DFT (con relleno de ceros)

% --- Ventanas (simetricas, igual que hanning/hamming/blackman de MATLAB) ---
nombres = {'Rectangular','Hanning','Hamming','Blackman'};
vent = {@(L) ones(L,1), ...
        @(L) 0.5-0.5*cos(2*pi*(0:L-1)'/(L-1)), ...
        @(L) 0.54-0.46*cos(2*pi*(0:L-1)'/(L-1)), ...
        @(L) 0.42-0.5*cos(2*pi*(0:L-1)'/(L-1))+0.08*cos(4*pi*(0:L-1)'/(L-1))};
dB = @(p) 10*log10(max(p,1e-12));

% Periodograma modificado: I(w)=|DFT(w.*x)|^2/(L*U), U=mean(w^2)  (normalizacion UIT2)
perio = @(x,f) abs(fft(f(length(x)).*x(:),NFFT)).^2 / (length(x)*mean(f(length(x)).^2));
wax = (0:NFFT/2)'*2*pi/NFFT;              % eje de frecuencia [0,pi] rad/muestra

%% 0. Espectro de las ventanas (L=64)
L = 64; NW = 8192;
wv = (0:NW/2)'*2*pi/NW;
figure(1);
fprintf('\nTabla 1. Caracteristicas medidas (L=%d)\n', L);
fprintf('%-12s %-14s %-14s\n','Ventana','Ancho (rad)','Lobulo sec. (dB)');
for i = 1:4
    w = vent{i}(L);
    W = abs(fft(w,NW)); W = W(1:NW/2+1);
    Wdb = dB((W/max(W)).^2);
    k = find(diff(W) > 0, 1);             % primer nulo tras el pico en w=0
    ancho = 2*wv(k);                      % ancho null-to-null del lobulo principal
    sec = max(Wdb(k:end));                % mayor lobulo secundario
    fprintf('%-12s %-14.4f %-14.1f\n', nombres{i}, ancho, sec);
    subplot(2,1,1); plot(0:L-1, w); hold on;
    subplot(2,1,2); plot(wv, Wdb); hold on;
end
subplot(2,1,1); title('Figura 1a. Ventanas en el tiempo (L=64)'); xlabel('n (muestras)'); ylabel('w[n]'); legend(nombres); grid on;
subplot(2,1,2); title('Figura 1b. Espectro de las ventanas (L=64)'); xlabel('Frecuencia (rad/muestra)'); ylabel('|W|^2 normalizado (dB)'); ylim([-120 5]); legend(nombres); grid on;
print -dpng Figura_1_ventanas.png

%% 1. Dos sinusoides + ruido blanco (media 0, varianza 1)
N = 128;                                  % longitud de la secuencia
w1 = 2*pi/12; w2 = 2*pi/14;
n = (0:N-1)';
ruido = randn(N,1);

amps = [1 1; 1 0.01];                     % 1.a: A1=A2=1   1.b: A1=1, A2=0.01
titulos = {'Figura 2. Periodograma, A1=A2=1','Figura 3. Periodograma, A1=1, A2=0.01'};
archivos = {'Figura_2_exp1a.png','Figura_3_exp1b.png'};
for e = 1:2
    x = amps(e,1)*cos(w1*n) + amps(e,2)*cos(w2*n) + ruido;
    figure(1+e);
    for i = 1:4
        subplot(2,2,i);
        P = perio(x,vent{i}); P = P(1:NFFT/2+1);   % solo [0,pi]
        plot(wax, dB(P)); hold on;
        plot([w1 w1],[-30 30],'r--'); plot([w2 w2],[-30 30],'g--');
        title(nombres{i}); xlabel('Frecuencia (rad/muestra)'); ylabel('DEP estimada (dB)');
        ylim([-30 30]); grid on;
    end
    sgtitle_txt = sprintf('%s, N=%d (rojo: w1=2pi/12, verde: w2=2pi/14)', titulos{e}, N);
    disp(sgtitle_txt);                    % sgtitle no existe en Octave: se imprime el titulo
    print('-dpng', archivos{e});
end

%% 2. Ruido blanco: periodograma vs periodograma promedio (solape 50 %)
Q = 1024; r = randn(Q,1);
Ls = [1024 256 64];
tit = {'Periodograma (L=1024, K=1)','Promedio L=256 (K=7)','Promedio L=64 (K=31)'};
figure(4);
fprintf('\nTabla 2. Ruido blanco (varianza 1)\n');
fprintf('%-12s %-28s %-8s %-8s\n','Ventana','Estimador','Media','Varianza');
for j = 1:3
    Lseg = Ls(j); R = Lseg/2; if Lseg == Q, R = Q; end   % paso entre segmentos
    starts = 1:R:(Q-Lseg+1);
    subplot(3,1,j); hold on;
    for i = 1:4
        P = zeros(NFFT,1);
        for s = starts
            P = P + perio(r(s:s+Lseg-1), vent{i});
        end
        P = P/numel(starts); P = P(1:NFFT/2+1);
        plot(wax, dB(P));
        fprintf('%-12s %-28s %-8.3f %-8.3f\n', nombres{i}, tit{j}, mean(P), var(P));
    end
    title(tit{j}); xlabel('Frecuencia (rad/muestra)'); ylabel('DEP estimada (dB)'); ylim([-25 15]); grid on;
    if j == 1, legend(nombres); end
end
print -dpng Figura_4_exp2.png

