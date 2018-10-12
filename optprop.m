function [musp, mua] = optprop(lambda, params)
% 
% OPTPROP gives musp(lambda) and mua(lambda) for given parameters.
%
% INPUT
%	lambda = vector of wavelengths to evaluate (nm)
%	params = vector of input parameters, given below:
%		1) Oxygen Saturation			(%)
%		2) concentration Hb				(mg/ml)
%		3) concentration melanin		(mg/ml)
%		4) Reduced Scattering at 630nm	(cm^-1)
%		5) Reduced Scattering Exponent  (unitless)
%		6) Vessel Radius				(cm)
%       7) Concentration Beta-Carotene  (mg/ml)
%
% OUTPUT 
%   musp			= vector of musp(lambda)
%   mua				= vector of mua_d(lambda)
%
% Written by Ricky Hennessy
% Please cite J. Biomed. Opt. 18(3), 037003

global Hb mel betaCarotene

%% Hemoglobin Data
Hbw = Hb(:,1);				% Wavelengths
Hbe_o = Hb(:,2);			% Oxy
Hbe_d = Hb(:,3);			% Deoxy

%% Get oxy & deoxy mua
Hbext_o     = interp1(Hbw,Hbe_o,lambda);			% Oxy Extinction Coeff
Hbext_d     = interp1(Hbw,Hbe_d,lambda);			% Deoxy Extinction Coeff
mua_oxy     = (log(10) * Hbext_o) / 64458;	
mua_deoxy   = (log(10) * Hbext_d) / 64458;

%% Total blood mua
mua_blood	= params(1)*mua_oxy + (1-params(1))*mua_deoxy;

%% Pigment Packaging
C_pack		= (1-exp(-2*mua_blood*params(6)))./(2*mua_blood*params(6));

%% Absorption due to hemoglobin
mua_hb		= params(2) * (C_pack .* mua_blood);

%% Absorption due to melanin
mel_ext		= interp1(mel(:,1),mel(:,2),lambda);
mua_mel		= params(3) * mel_ext;

%% Absorption due to betaCarotene
bC_ext		= interp1(betaCarotene(:,1),betaCarotene(:,2),lambda);
mua_bC		= params(7) * (log(10) * bC_ext) / 536.89;

%% Get mua
mua			= mua_hb + mua_mel + mua_bC;
	
%% Get musp
musp		= params(4)*(lambda/630).^(params(5));