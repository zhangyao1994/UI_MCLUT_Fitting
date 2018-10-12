function E = MC_LUT_fitness(params)
% function E = MC_LUT_fitness(X)
% 
% MC_LUT_FITNESS calculate the sum of squares error between the modeled
%	spectra and the measured spectra with X as the vector of parameters
%	used to create the modeled spectra
%
% INPUT
%	params	= vector of input parameters, given below:
%		1) Oxygen Saturation			(%)
%		2) Concentration Hb				(mg/ml)
%		3) concentration melanin		(mg/ml)
%		4) Reduced Scattering at 630nm	(cm^-1)
%		5) Reduced Scattering Exponent  (unitless)
%		6) Vessel Radius				(cm)
%
% OUTPUT 
%   E   - Error between measured and modeled spectra
%
% Written by Ricky Hennessy
% Please cite J. Biomed. Opt. 18(3), 037003

%% Globals
global spectra Fig1 F_meanRatio
warning off

%% Create Model
R = MC_LUT_forward(spectra(:,1),params);
S = spectra(:,2);
R(:,2)=R(:,2).*F_meanRatio;

%% Find Error
dif = abs(S-R(:,2))./S;
E = mean(dif) * 100;

%% Plot
if Fig1
    plot(R(:,1),R(:,2),spectra(:,1),spectra(:,2))
	xlabel('wavelength (nm)')
	ylabel('R')
	axis([min(spectra(:,1)) max(spectra(:,1)) min(S)-.01 max(S)+.01])
    hold off
    drawnow
end

