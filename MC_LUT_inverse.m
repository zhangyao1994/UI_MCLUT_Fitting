function [R, params] = MC_LUT_inverse(S)
% function E = MC_LUT_inverse(S)
% 
% MC_LUT_INVERSE extracts the parameters from the spectrum S using
%	a LUT based inverse model
%
% INPUT
%	S	- Measured spectra: S(:,1) = wavelengths in nanometers. 
%							S(:,2) = reflectance values.
%
% OUTPUT 
%   R		- Modeled Spectra (fit)
%	params  - Extracted parameters:
%		1) Oxygen Saturation			(%)
%		2) Concentration Hb				(mg/ml)
%		3) concentration melanin		(mg/ml)
%		4) Reduced Scattering at 630nm	(cm^-1)
%		5) Reduced Scattering Exponent  (unitless)
%		6) Vessel Radius				(cm)
%       7) Concentration Beta-Carotene  (mg/ml)
%  
% Written by Ricky Hennessy
% Please cite J. Biomed. Opt. 18(3), 037003


%% Global Variables
global spectra Fig2 F_meanRatio
spectra = S;

%% Initialize % not sure about the Beta-Carotene concentration
	 %SO2	[Hb]	[mel]	mus630  B		r_vess    [BC]
X0 = [.50	1		.4		20		-1.3	0.0020  0];
lb = [0		0		0.0		5		-1.3	0.0001  0];
ub = [1		10	    2	    50	    -.7		0.3000  0.15];

%% Optimization
options = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolX',5E-8, ...
	'TOlFun',5E-8,'display','off','algorithm','interior-point');
params = fmincon(@MC_LUT_fitness,X0,[],[],[],[],lb,ub,[],options);

%% get modeled spectrum
R = MC_LUT_forward(spectra(:,1),params);
if Fig2
    plot(R(:,1),R(:,2)*F_meanRatio,S(:,1),S(:,2))
	xlabel('wavelength (nm)')
	ylabel('R')
	axis([min(spectra(:,1)) max(spectra(:,1)) min(S(:,2))-.01 max(S(:,2))+.01])
    hold off
    drawnow
end
