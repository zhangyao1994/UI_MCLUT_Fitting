 function UI
 % Yao Updated 11292016 Add Fit button
    clc;clear
    global LUT mua_v musp_v Hb mel betaCarotene F_meanRatio Fig1 Fig2 Spectra
    % F_meanRatio=1/0.6068;% obtained during LUT Calibration
    % F_meanRatio=0.0127; % LUT_350_100
    % F_meanRatio=0.0091;% From Hieu, testing
    % F_meanRatio=0.0085; % LUT_320_100 Austin's beads data
    % F_meanRatio=0.025; % try
    warning off
     % Parse simulation output files
    [FileName_ref,PathName_ref] = uigetfile('*.mat','Select .mat file to plot');
    load([PathName_ref FileName_ref]);
    
    cd LUT
    load LUT_320_100.mat
    LUT; mua_v; musp_v;
    cd ..
    cd Chromophores
    Hb = importdata('Hb.mat');  % http://omlc.ogi.edu/spectra/hemoglobin/summary.html
    mel = importdata('mel.mat'); % http://omlc.ogi.edu/spectra/melanin/eumelanin.html
    betaCarotene = importdata('betaCarotene.mat'); % http://omlc.org/spectra/PhotochemCAD/data/041-abs.txt
    cd ..
    
    wvl = linspace(410,650,300);% 350,700
    Spectra(:,2)=interpn(S(:,1),S(:,2),wvl);
    Spectra(:,1)=wvl;
    
if 0 % Did calibration first and saved the ratio in data.mat
    [FileName_ref,PathName_ref] = uigetfile('*.mat','Select Calibration data');
    load ([PathName_ref FileName_ref]);
    R = MC_LUT_forward(wvl, [0 0 0 10 -1.3 .0001 0]); % Polystyrene beads info
    C = interp1(CalibRef(:,1),CalibRef(:,2),wvl); %1 mm-1 reduced scattering coefficient at 630 nm
    F = C./R(:,2)';
    F_meanRatio=mean(F);
    Fig0=1;
    if Fig0
        figure(100)
        plot(wvl,R(:,2).*mean(F)','r','linewidth',2)
        hold on
        plot(wvl,C,'k--','linewidth',4)
        set(gca,'fontsize',16)
        xlabel('wavelength (nm)','fontsize',16)
        ylabel('R','fontsize',16)
        % axis([min(CalibRef(:,1)) max(CalibRef(:,1)) min(C)-0.02 max(C)+0.02])
        hold on
        drawnow
    end
end
%		1) Oxygen Saturation			(%)
%		2) Concentration Hb				(mg/ml)
%		3) concentration melanin		(mg/ml)
%		4) Reduced Scattering at 6MaxRnm	(cm^-1)
%		5) Reduced Scattering Exponent  (unitless)
%		6) Vessel Radius				(cm)
%       7) Concentration Beta-Carotene  (mg/ml)

    % set some default values
%     S02 = 0.9;
%     CHb = 0;
%     Cmel = 0;
%     musp = 15;
%     B = -1;
%     VesselR = .0001;
%     CBC = 0;        

    % Create a figure and axes
    close all
    f = figure('Visible','off');
    ax = axes('position',[0.13 0.39  0.77 0.54]);
    
    Fig1=0;Fig2=0;
    [R param] = MC_LUT_inverse(Spectra);
    S02 = param(1);
    CHb = param(2);
    Cmel = param(3);
    musp = param(4);
    B = param(5);
    VesselR = param(6);
    CBC = param(7);
    R = MC_LUT_forward(wvl, [S02 CHb Cmel musp B VesselR CBC]);

    plot(R(:,1),R(:,2)*F_meanRatio,S(:,1),S(:,2))
    xlabel('wavelength (nm)')
    ylabel('R')
    % axis([min(S(:,1)) max(S(:,1)) max(min(S(:,2))-.01,0) min(max(S(:,2))+.01,0.5)])
    title('Modeled DRS')
    drawnow
    xlim([min(wvl) max(wvl)])
    
    % Create push button
    btn = uicontrol('Style','pushbutton','String','Fit',...
        'Units','Normalized','Position',[0.922,0.898,0.06,0.069],...
        'Callback',@fit);
    
   % Create slider
   %        1) Oxygen Saturation			(%)
   %%
    sld1 = uicontrol('Style', 'slider',...
        'Min',0,'Max',1,'Value',S02,...
        'Units','Normalized','Position', [0.13 0.2 0.23 0.05],...
        'Callback', @surfzlim1); 
    %     % Add a text uicontrol to label the slider.
    txt1 = uicontrol('Style','text',...
        'Units','Normalized','Position',[0.13 0.25 0.23 0.05],...
        'String','Oxygen Saturation');
    %       2) Concentration Hb				(mg/ml)
    sld2 = uicontrol('Style', 'slider',...
        'Min',0,'Max',10,'Value',CHb,...
        'Units','Normalized','Position', [0.13 0.1 0.23 0.05],...
        'Callback', @surfzlim2);
    txt2 = uicontrol('Style','text',...
        'Units','Normalized','Position',[0.13 0.15 0.23 0.05],...
        'String','[Hb]');
    %       3) concentration melanin		(mg/ml)
    sld3 = uicontrol('Style', 'slider',...
        'Min',0,'Max',10,'Value',Cmel,...
        'Units','Normalized','Position', [0.13 0.0 0.23 0.05],...
        'Callback', @surfzlim3);
    txt3 = uicontrol('Style','text',...
        'Units','Normalized','Position',[0.13 0.05 0.23 0.05],...
        'String','[mel]');
    %       4) Reduced Scattering at 630nm	(cm^-1)
    sld4 = uicontrol('Style', 'slider',...
        'Min',5,'Max',70,'Value',musp,...
        'Units','Normalized','Position', [0.4 0.2 0.23 0.05],...
        'Callback', @surfzlim4);
    txt4 = uicontrol('Style','text',...
        'Units','Normalized','Position',[0.4 0.25 0.23 0.05],...
        'String','musp630');
    %       5) Reduced Scattering Exponent  (unitless)
    sld5 = uicontrol('Style', 'slider',...
        'Min',-3,'Max',-0.5,'Value',B,...
        'Units','Normalized','Position', [0.4 0.1 0.23 0.05],...
        'Callback', @surfzlim5);
    txt5 = uicontrol('Style','text',...
        'Units','Normalized','Position',[0.4 0.15 0.23 0.05],...
        'String','B');
    %       6) Vessel Radius				(cm)
    sld6 = uicontrol('Style', 'slider',...
        'Min',0.0001,'Max',.3,'Value',VesselR,...
        'Units','Normalized','Position', [0.4 0.0 0.23 0.05],...
        'Callback', @surfzlim6);
    txt6 = uicontrol('Style','text',...
        'Units','Normalized','Position',[0.4 0.05 0.23 0.05],...
        'String','Vessel Radius');
    %       7) Concentration Beta-Carotene  (mg/ml)
    sld7 = uicontrol('Style', 'slider',...
        'Min',0,'Max',0.2,'Value',CBC,...
        'Units','Normalized','Position', [0.67 0.2 0.23 0.05],...
        'Callback', @surfzlim7);
    txt7 = uicontrol('Style','text',...
        'Units','Normalized','Position',[0.67 0.25 0.23 0.05],...
        'String','[Beta-Carotene]');
    
    % Make figure visble after adding all components
    f.Visible = 'on';    

    function fit(source,event)
        Fig1=1;Fig2=1;
        [R param] = MC_LUT_inverse(Spectra);
        S02 = param(1);
        CHb = param(2);
        Cmel = param(3);
        musp = param(4);
        B = param(5);
        VesselR = param(6);
        CBC = param(7);
        hold off    
        plot(R(:,1),R(:,2)*F_meanRatio,S(:,1),S(:,2))
        xlabel('wavelength (nm)')
        ylabel('R')
        % axis([min(S(:,1)) max(S(:,1)) max(min(S(:,2))-.01,0) min(max(S(:,2))+.01,0.5)])
        title('Modeled DRS')
        xlim([min(wvl) max(wvl)])
        drawnow
    end
    
    function surfzlim1(source,event)
        S02 = source.Value;
        display(S02)
        R = MC_LUT_forward(wvl, [S02 CHb Cmel musp B VesselR CBC]);
        hold off
        plot(wvl,R(:,2)*F_meanRatio,S(:,1),S(:,2))
        xlabel('wavelength (nm)')
        ylabel('R')
        title('Modeled DRS')
        xlim([min(wvl) max(wvl)])
        %MaxR=1.2*max(R(:,2));
        % axis([min(S(:,1)) max(S(:,1)) max(min(S(:,2))-.01,0) min(max(S(:,2))+.01,0.5)])
        drawnow
    end

    function surfzlim2(source,event)
        CHb = source.Value;
        display(CHb)
        R = MC_LUT_forward(wvl, [S02 CHb Cmel musp B VesselR CBC]);
        hold off
        plot(wvl,R(:,2)*F_meanRatio,S(:,1),S(:,2))
        xlabel('wavelength (nm)')
        ylabel('R')
        title('Modeled DRS')        
        xlim([min(wvl) max(wvl)])
        %MaxR=1.2*max(R(:,2));
        % axis([min(S(:,1)) max(S(:,1)) max(min(S(:,2))-.01,0) min(max(S(:,2))+.01,0.5)])
        drawnow        
    end

    function surfzlim3(source,event)
        Cmel = source.Value;
        display(Cmel)
        R = MC_LUT_forward(wvl, [S02 CHb Cmel musp B VesselR CBC]);
        hold off
        plot(wvl,R(:,2)*F_meanRatio,S(:,1),S(:,2))
        xlabel('wavelength (nm)')
        ylabel('R')
        title('Modeled DRS')        
        xlim([min(wvl) max(wvl)])
        %MaxR=1.2*max(R(:,2));
        % axis([min(S(:,1)) max(S(:,1)) max(min(S(:,2))-.01,0) min(max(S(:,2))+.01,0.5)])
        drawnow        
    end

    function surfzlim4(source,event)
        musp = source.Value;
        display(musp)
        R = MC_LUT_forward(wvl, [S02 CHb Cmel musp B VesselR CBC]);
        hold off
        plot(wvl,R(:,2)*F_meanRatio,S(:,1),S(:,2))
        xlabel('wavelength (nm)')
        ylabel('R')
        title('Modeled DRS')        
        xlim([min(wvl) max(wvl)])
        %MaxR=1.2*max(R(:,2));
        % axis([min(S(:,1)) max(S(:,1)) max(min(S(:,2))-.01,0) min(max(S(:,2))+.01,0.5)])
        drawnow        
    end

    function surfzlim5(source,event)
        B = source.Value;
        display(B)
        R = MC_LUT_forward(wvl, [S02 CHb Cmel musp B VesselR CBC]);
        hold off
        plot(wvl,R(:,2)*F_meanRatio,S(:,1),S(:,2))
        xlabel('wavelength (nm)')
        ylabel('R')
        title('Modeled DRS')        
        xlim([min(wvl) max(wvl)])
        %MaxR=1.2*max(R(:,2));
        % axis([min(S(:,1)) max(S(:,1)) max(min(S(:,2))-.01,0) min(max(S(:,2))+.01,0.5)])
        drawnow        
    end

    function surfzlim6(source,event)
        VesselR = source.Value;
        display(VesselR)
        R = MC_LUT_forward(wvl, [S02 CHb Cmel musp B VesselR CBC]);
        hold off
        plot(wvl,R(:,2)*F_meanRatio,S(:,1),S(:,2))
        xlabel('wavelength (nm)')
        ylabel('R')
        title('Modeled DRS')        
        xlim([min(wvl) max(wvl)])
        %MaxR=1.2*max(R(:,2));
        % axis([min(S(:,1)) max(S(:,1)) max(min(S(:,2))-.01,0) min(max(S(:,2))+.01,0.5)])
        drawnow        
    end

    function surfzlim7(source,event)
        CBC = source.Value;
        display(CBC)
        R = MC_LUT_forward(wvl, [S02 CHb Cmel musp B VesselR CBC]);
        hold off
        plot(wvl,R(:,2)*F_meanRatio,S(:,1),S(:,2))
        xlabel('wavelength (nm)')
        ylabel('R')
        title('Modeled DRS')        
        xlim([min(wvl) max(wvl)])
        % axis([min(S(:,1)) max(S(:,1)) max(min(S(:,2))-.01,0) min(max(S(:,2))+.01,0.5)])
        drawnow        
    end
end