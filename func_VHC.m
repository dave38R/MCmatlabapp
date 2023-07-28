function VHC = func_VHC(T)
t = T./1000;
Cp = 52.1743 + 87.951*t - 84.2411*t.^2 + 31.542*t.^3 - 2.6334*t.^4 - 0.71391*t.^(-2); % J/mol/K

    function L = mechanical_dilatation(T)
        L = zeros(size(T));
        idx1 = T >= 273.15 & T <= 923;
        idx2 = T > 923 & T <= 3120;
        L(idx1) = 9.973e-1 + 9.082e-6*T(idx1) - 2.705e-10*T(idx1).^2 + 4.391e-13*T(idx1).^3;
        L(idx2) = 9.9672e-1 + 1.179e-5*T(idx2) - 2.429e-9*T(idx2).^2 + 1.219e-12*T(idx2).^3;
        if any(~(idx1 | idx2))
            error('Temperature must be between 273.15K and 3120K')
        end
    end

L = mechanical_dilatation(T); % no dimensions


    function density = func_density(L)
        density = 10.963*L.^3; % g/cm^3
    end

density = func_density(L); % g/cm^3
M = 0.2700277; % kg/mol
VHC = Cp.*density./(1000*M); % J/cm^3/K
end
































