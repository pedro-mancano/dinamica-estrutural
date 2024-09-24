classdef Newmark
    % Classe que implementa o método de Newmark para integração dinâmica.
    
    properties
        m        % Massa do sistema
        c        % Coeficiente de amortecimento
        k        % Rigidez do sistema
        u0       % Deslocamento inicial
        du0      % Velocidade inicial
        h        % Passo de tempo
        f0       % Função de força externa
        Beta     % Parâmetro de Newmark
        Gamma    % Parâmetro de Newmark
        t        % Tempo atual
    end

    methods

        function obj = Newmark(m, c, k, u0, du0, h, f0, Beta, Gamma)
            % Construtor da classe Newmark
            % Inicializa os parâmetros do sistema e do método de Newmark.
            obj.m = m;               % Massa
            obj.c = c;               % Amortecimento
            obj.k = k;               % Rigidez
            obj.u0 = u0;             % Deslocamento inicial
            obj.du0 = du0;           % Velocidade inicial
            obj.h = h;               % Passo de tempo
            obj.f0 = f0;             % Função de força
            obj.t = 0;               % Tempo inicial
            obj.Beta = Beta;         % Beta para o método de Newmark
            obj.Gamma = Gamma;       % Gamma para o método de Newmark
        end

        function data = integrate_until(obj, time_final)
            % Integra até um tempo final especificado usando o método de Newmark.
            
            % Calcula a aceleração inicial
            ddU0 = (1 / obj.m) * (obj.f0(obj.t) - obj.c * obj.du0 - obj.k * obj.u0);
            obj.t = obj.h;           % Atualiza o tempo
            i = 0;                   % Inicializa o contador de passos
            data = obj.u0;          % Armazena o deslocamento inicial

            % Loop de integração
            while obj.t < time_final
                i = i + 1;          % Incrementa o contador

                % Calcula o deslocamento e a velocidade do passo anterior
                Us_1 = obj.u0 + obj.h * obj.du0 + (0.5 - obj.Beta) * obj.h ^ 2 * ddU0;
                dUs_1 = obj.du0 + (1 - obj.Gamma) * obj.h * ddU0;

                % Monta o sistema de equações
                M = obj.m + obj.c * obj.Gamma * obj.h + obj.k * obj.Beta * obj.h ^ 2;
                F = obj.f0(obj.t) - obj.c * dUs_1 - obj.k * Us_1;
                ddUs_1 = F / M;     % Calcula a aceleração no passo atual

                % Atualiza o deslocamento e a velocidade
                U = Us_1 + obj.Beta * obj.h ^ 2 * ddUs_1;
                dU = dUs_1 + obj.Gamma * obj.h * ddUs_1;

                % Atualiza os estados do objeto
                obj.u0 = U;         % Atualiza o deslocamento
                obj.du0 = dU;       % Atualiza a velocidade

                ddU0 = ddUs_1;      % Armazena a aceleração para o próximo passo

                data = [data; U];   % Armazena o deslocamento atual

                obj.t = i * obj.h;  % Atualiza o tempo
            end
        end
    end
end
