function R = MC_LUT_forward(lambda, params)
% function E = MC_LUT_inverse(S)
% 
% MC_LUT_FORWARD is the LUT based forward model that generates a spectrum
%	at the wavelength specified in LAMBDA with the parameters specified 
%	in PARAMS
%
% INPUT
%	lambda	- wavelength (nm)
%	params	- parameters:
%		1) Oxygen Saturation			(%)
%		2) Concentration Hb				(mg/ml)
%		3) concentration melanin		(mg/ml)
%		4) Reduced Scattering at 630nm	(cm^-1)
%		5) Reduced Scattering Exponent  (unitless)
%		6) Vessel Radius				(cm)
%       7) Concentration Beta-Carotene  (mg/ml)
%
% OUTPUT 
%   R		- Modeled Spectrum
%
% Written by Ricky Hennessy
% Please cite J. Biomed. Opt. 18(3), 037003

%% Globals
global LUT mua_v musp_v

%% Find musp & mua
[musp, mua] = optprop(lambda, params);

%% Prevent Values Outside LUT
flag=0;
if min(musp)<min(musp_v)
    flag=1;
end
if max(musp)>max(musp_v)
    flag=1;
end
if max(mua)>max(mua_v)
    flag=1;
end
musp(musp < min(musp_v))	= min(musp_v);
musp(musp > max(musp_v))	= max(musp_v);

mua(mua < min(mua_v))		= min(mua_v);
mua(mua > max(mua_v))		= max(mua_v);

%% Create Spectra
R = zeros(length(lambda),2);
R(:,1) = lambda;
R(:,2) = interp2(mua_v, musp_v, LUT', mua, musp,'cubic');
