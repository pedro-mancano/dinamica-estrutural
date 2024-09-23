classdef Newmark

    properties
        m
        c
        k
        u0
        du0
        h
        f0
        Beta
        Gamma
        t
    end

    methods

        function obj = Newmark(m, c, k, u0, du0, h, f0, Beta, Gamma)
            obj.m = m;
            obj.c = c;
            obj.k = k;
            obj.u0 = u0;
            obj.du0 = du0;
            obj.h = h;
            obj.f0 = f0;
            obj.t = 0;
            obj.Beta = Beta;
            obj.Gamma = Gamma;
        end

        function data = integrate_until(obj, time_final)
            ddU0 = (1 / obj.m) * (obj.f0(obj.t) - obj.c * obj.du0 - obj.k * obj.u0);
            obj.t = obj.h;
            i = 0;
            data = obj.u0;

            while obj.t < time_final
                i = i + 1;

                Us_1 = obj.u0 + obj.h * obj.du0 + (0.5 - obj.Beta) * obj.h ^ 2 * ddU0;
                dUs_1 = obj.du0 + (1 - obj.Gamma) * obj.h * ddU0;

                M = obj.m + obj.c * obj.Gamma * obj.h + obj.k * obj.Beta * obj.h ^ 2;
                F = obj.f0(obj.t) - obj.c * dUs_1 - obj.k * Us_1;
                ddUs_1 = F / M;

                U = Us_1 + obj.Beta * obj.h ^ 2 * ddUs_1;
                dU = dUs_1 + obj.Gamma * obj.h * ddUs_1;

                obj.u0 = U;
                obj.du0 = dU;

                ddU0 = ddUs_1;

                data = [data; U];

                obj.t = i * obj.h;
            end

        end

    end

end
