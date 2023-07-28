function [Gxyz,NA, intensity] = sim(resolution, cuboid_length, cuboid_height, sample_radius, sample_height, wavelength, mua, mus, n, T, super, p)

MCmatlab.closeMCmatlabFigures();
model = MCmatlab.model;
model.G.nx                = resolution; % Number of bins in the x direction
model.G.ny                = resolution; % Number of bins in the y direction
model.G.nz                = resolution; % Number of bins in the z direction
model.G.Lx                = cuboid_length; % [cm] x size of simulation cuboid
model.G.Ly                = cuboid_length; % [cm] y size of simulation cuboid
model.G.Lz                = cuboid_height; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

%% Monte Carlo simulation
model.MC.simulationTimeRequested  = .05; % [min] Time duration of the simulation

model.MC.matchedInterfaces        = false; % Does not assume all refractive indices are the same
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.wavelength               = wavelength; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.lightSource.sourceType   = 4; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
% Create custom radial distribution representing an annulus
inner_radius = 0.1; % [cm] inner radius of annulus
outer_radius = 0.3; % [cm] outer radius of annulus
sigma = 0.05; % [cm] standard deviation of Gaussian distribution
r = linspace(0, outer_radius, 1000);

% In case super is true, we use a super gaussian beam, else we use a annular gaussian beam
if super
    intensity = 4^(1/p)*p/(2*pi*outer_radius*outer_radius*gamma(2/p))*exp(-2*(r/(outer_radius)).^p);
else
    intensity = exp(-(r - inner_radius).^2 / (2 * sigma^2));
end

model.MC.lightSource.focalPlaneIntensityDistribution.radialDistr = intensity; % Radial focal plane intensity distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth = outer_radius; % [cm] Radial focal plane 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
model.MC.lightSource.angularIntensityDistribution.radialDistr = intensity; % Radial angular intensity distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.angularIntensityDistribution.radialWidth = outer_radius; % [rad] Radial angular half-width if top-hat or Gaussian or half-angle of the full distribution if custom. For a diffraction limited Gaussian beam, this should be set to model.MC.wavelength*1e-9/(pi*model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth*1e-2))
% model.MC.wavelength*1e-9/(pi*model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth*1e-2))
model.MC.lightSource.xFocus       = 0; % [cm] x position of focus
model.MC.lightSource.yFocus       = 0; % [cm] y position of focus
model.MC.lightSource.zFocus       = 0; % [cm] z position of focus
model.MC.lightSource.theta        = 0; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = 0; % [rad] Azimuthal angle of beam center axis

model = runMonteCarlo(model);
model.G = model.G.updateGeometry;
G = model.G;
MCorFMC = model.MC;
NA = MCorFMC.NA;
writematrix(NA,'NA.csv');
Gxyz = [G.x; G.y; G.z];
writematrix(Gxyz,'Gxyz.csv');

%% Thermal functions
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

    function TC = func_TC(T)
        t = T/1000;
        TC = 1/(7.5408 + 17.962*t + 3.6142*t^2) + 64/t^2.5*exp(-16.35/t); % W/cm/K
    end


%% Geometry function(s) (see readme for details)
    function M = geometryDefinition(X,Y,Z,parameters)
        middle = cuboid_height/2;
        M = ones(size(X)); % fill background with air
        M(X.^2 + Y.^2 < sample_radius^2) = 2; % UO2
        M(Z < middle - sample_height/2) = 1; % air
        M(Z > middle + sample_height/2) = 1; % air

    end
%% Media Properties function (see readme for details)
    function mediaProperties = mediaPropertiesFunc(parameters)
        mediaProperties = MCmatlab.mediumProperties;

        j=1;
        mediaProperties(j).name  = 'air';
        mediaProperties(j).mua   = 1e-8; % Absorption coefficient [cm^-1]
        mediaProperties(j).mus   = 1e-8; % Scattering coefficient [cm^-1]
        mediaProperties(j).g     = 1; % Henyey-Greenstein scattering anisotropy
        mediaProperties(j).VHC = 3391*1.109e-3; % [J cm^-3 K^-1]
        mediaProperties(j).TC  = 0.37e-2; % [W cm^-1 K^-1]
        mediaProperties(j).n = 1;

        j=2;
        mediaProperties(j).name = 'UO2';
        mediaProperties(j).g = 1;
        mediaProperties(j).mua = mua; % Absorption coefficient [cm^-1]
        mediaProperties(j).mus = mus; % Scattering coefficient [cm^-1]
        % To use these functions you need to have defined a temperature
        mediaProperties(j).VHC = func_VHC(T);
        mediaProperties(j).TC = func_TC(T);
        mediaProperties(j).n = n;
    end
end