% =========================================================================
% GRÁFICA DE v(t) - EJE X CLARO
% =========================================================================

% --- 1. Definición de Parámetros y Constantes ---
R = 100;    % Ohms
L = 10e-3;  % Henrios
C = 10e-6;  % Faradios
A = 6.0;    % Voltios
B = -9.0;   % Voltios

% Constantes calculadas
alpha = 1/(2*R*C); % Factor de decaimiento (5000)
omega_d = sqrt(1/(L*C) - alpha^2); % Frecuencia amortiguada (8660.25404)

% --- 2. Definición de la Función v(t) ---
v = @(t) exp(-alpha.*t) .* (A.*cos(omega_d.*t) + B.*sin(omega_d.*t));

% 3. Vector de tiempo de 0 a 1.5ms (0.0015s)
t = linspace(0, 0.0015, 1000); % Aumentado a 1000 puntos para suavidad

% 4. Graficar y configurar
figure;
plot(t, v(t), 'b', 'LineWidth', 2);
hold on;

% Dibujar el Eje X (y=0) para ver los cortes
plot([0,0.01], [0, 0], 'k-');

% Etiquetas y títulos
xlabel('Tiempo (s)', 'fontsize', 12);
ylabel('Voltaje v(t) (V)', 'fontsize', 12);
title('Voltaje Subamortiguado v(t) en Circuito RLC', 'fontsize', 14);

% Aplicar límites
ylim([-10, 10]);
xlim([0, 0.0005]);

grid on;
hold off;
