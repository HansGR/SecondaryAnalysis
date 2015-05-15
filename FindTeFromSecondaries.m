function output = FindTeFromSecondaries(Y1, Y1unc, Y2n, Y2nunc, Y2p, Y2punc, ...
    Ti, rho, Z0, A0, f0, Te, HSorUni, varargin)
% Find rhoR & Te that satisfies conditions set by Y2p|Y2n, Y1, Ti, rho
% You must declare the range of Te to test (N tests)
% Inputs:  Y1 = primary yield (1x1)
%          Y1unc = uncertainty in primary yield (1x1)
%          Y2n = secondary neutron yield (1x1)
%          Y2nunc = secondary neutron yield uncertainty (1x1)
%          Y2p = secondary proton yield (1x1)
%          Y2punc = secondary proton yield uncertainty (1x1)
%          Ti = Ion temperature (keV), (1x1)
%          rho_0 = plasma density at BT for clean fuel (g/cc), (1x1)
%          Z0 = charge state of fuel ions, 1xM
%          A0 = Mass of fuel ions (amu), 1xM
%          f0 = number fractions of fuel ions, NxM
%          Te = electron temperature (keV), (Nx1)
%          HSorUni = flag, return hotspot (0) or uniform model (1)
% Outputs: output.rhoR       = fuel rhoR (mg/cc) best fit
%          output.rhoR_range = fuel rhoR (mg/cc) chi^2 - min(chi^2) <= 1 range 
%          output.chisq      = chi^2 of the best fit
%          output.rhoR_D     = deuterium rhoR (mg/cc) best fit
%          output.rhoR_D_range = deuterium rhoR (mg/cc) chi^2 - min(chi^2) <= 1 range 
%          output.Te         = Electron temperature (keV), best fit
%          output.Te_range   = Electron temperature (keV), chi^2 - min(chi^2) <= 1 range 

%% Test Values
%{
[Y1, Y1unc, Y2n, Y2nunc, Y2p, Y2punc, Ti_exp, Ti_exp_unc, rho, rhounc] = deal(...
     22e12, 2e12, 13.5e10, 1.3e10, 14e9, 7e9, 3.6, 0.2, 2.96, 0.89);  %N130813 2-shock HDC

%     5.8e11, 0.4e11, 2.8e8, 0.3e8, 2.1e8, 0.4e8, 4.4, 0.3, 0.42, 0.284);  %, N120728
%     6.5e11, 0.4e11, 3.5e8, 0.3e8, 2.2e8, 0.4e8, 4.2, 0.3);  %, N120730
%     8.6e11, 0.4e11, 6.2e8, 0.6e8, 3.7e8, 0.7e8, 4.3, 0.4);  %, N121119
%     3.7e11, 0.2e11, 2.4e8, 0.2e8, 1.9e8, 0.4e8, 3.4, 0.2);  %, N121207
%     7.3e11, 0.2e11, 10.3e8, 1e8, 5.8e8, 1.2e8, 3.7, 0.2);  %, N130320
%     6e11, 0.2e11, 5.1e8, 0.5e8, 4e8, 0.8e8, 3.8, 0.2);  %, N130321
%%% NIF exploding pushers:
%     3.01e11, 1.13e10, 8.57e7, 1.0e7, 2.01e8, 0.4e8, 5.4, 0.2);   %, N110131
%     2.52e11, 0.184e11, 7.6e7, 2.8e7, 1.7e8, 0.3e8, 3.9, 0.3);    %, N130129

    
%   22e12, 2e12, 13.5e10, 1.3e10, 14e9, 7e9, 3.6, 0.2, 2.96, 0.89);  %N130813 2-shock HDC
%   5.15e12, 0.16e12, 9.62e9, 0.78e9, 6.6e9, 0.9e9, 3.77, 0.21, 1.4, 0.43);  %N130312 IndiExP
%}
%% Plasma
%{
%rho = 0.00162*5.604^3; %g/cm^3:    0.42; N120728;  1.4; N130312 IndiExpP;   2.96; %N130813 2-shock HDC
%rhounc = 0.25*rho; %       0.284; N120728; 0.43; N130312 IndiExpP;   0.89; %N130813 2-shock HDC

% Z = [1,2];
% A = [2,3];
% f = [1,0];
% 
% % mix modification:
% mixfrac = 0; %fraction of mix by atom, n_mix = n_fuel * mixfrac
% mix_Z = [1,6];
% mix_A = [1,12];
% mix_stoch = [1,1]; %stochacity of mix
%     mix_stoch = mix_stoch/sum(mix_stoch); %renorm
% rho_mix = rho./sum(f.*A).*(mixfrac*mix_stoch.*mix_A);

% rho = rho/sum(A.*f)*sum([A.*f,mix_A.*mix_stoch*mixfrac]);
% rhounc = rhounc/sum(A.*f)*sum([A.*f,mix_A.*mix_stoch*mixfrac]);

% Z = [Z,mix_Z];
% A = [A,mix_A];
% f = [f,mix_stoch*mixfrac];
%     f = f/sum(f); %renorm
% f_mass = f.*A/sum(f.*A);
%}
%% Stopping Power parameters
TeCorrection = 1;

quiet = 0; %suppress messages
makeFigures = 1;
makeChiFigures = 1;

%% Deal with variable arguments
if ~isempty(varargin),
    if mod(length(varargin),2)~=0,
        error('Extra arguments must be declared in pairs ''VariableName'', value');
    end
    
    for i = 1:length(varargin)/2,
        if isstr(varargin{i*2-1});
            arg = varargin{i*2-1};
            val = varargin{i*2};
            if isstr(val),
                eval([arg,'=',val,';'])
            else
                eval([arg,'=',num2str(val),';'])
            end
        end
    end
end

%% Constants (cgs)
e = 4.80320425*10^-10;  %statcoulombs, (erg*cm)^(1/2)
c = 2.99792458*10^10;   %cm/sec
me = 9.10938291*10^-28; %g
mp = 1.67262178*10^-24; %g
mpc2 = 938.272046;      %MeV, proton rest mass
kB = 1.3806488*10^-16;  %erg/K
kBkeV = 1.60217646*10^-9;%erg/keV
Na = 6.02214129*10^23;  %mol^-1
hbar = 1.05457148*10^-27; %erg*sec

%% Generate Plasma
N = length(Te);
M = size(f0,2);

Ti2 = repmat(Ti,[N,1]);
Zf2 = repmat(reshape(Z0,[1,M]),[N,1]);
Af2 = repmat(reshape(A0,[1,M]),[N,1]);
fi2 = repmat(reshape(f0,[1,M]),[N,1]);
    fi2 = fi2./repmat(sum(fi2,2),[1,M]);  %renorm

f_mass = fi2.*Af2./repmat(sum(fi2.*Af2,2),[1,size(fi2,2)]);

Dind = find((Z0==1).*(A0==2));

rho2 = repmat(rho,[N,1]);

%% Calculate the secondary yield curves for this plasmas
sec = FindSecondaryYield(Ti2, Te, rho2, Zf2, Af2, fi2, HSorUni,...
    'TeCorrection', TeCorrection, 'quiet', quiet, 'parallel', 1);

PtotDTn_u = sec.Y2n;
PtotD3Hep_u = sec.Y2p;
rhoRt_ext = sec.rhoR_Y2n;
rhoR3He_ext = sec.rhoR_Y2p;

%% Calculate the fuel rhoR and the electron temperature (from Y2n/Y2p)
rhoR_vector = 0:min(rhoR3He_ext(:,2)-rhoR3He_ext(:,1)):max(rhoRt_ext(:,end));

goodness = zeros([N,length(rhoR_vector)]);
goodness_p = zeros([N,length(PtotD3Hep_u(1,:))]);
goodness_n = zeros([N,length(PtotDTn_u(1,:))]);

% For each plasma Te...
for i = 1:N,    
    % Figure out how good the mix model is overall:
    goodness_p(i,:) = (PtotD3Hep_u(i,:)-(Y2p/Y1)).^2/(((Y2punc/Y2p)^2+(Y1unc/Y1)^2)*(Y2p/Y1)^2);
    goodness_n(i,:) = (PtotDTn_u(i,:)-(Y2n/Y1)).^2/(((Y2nunc/Y2n)^2+(Y1unc/Y1)^2)*(Y2n/Y1)^2);
    
    goodness(i,:) = interp1(rhoR3He_ext(i,:),goodness_p(i,:),rhoR_vector,'pchip') ...
        + interp1(rhoRt_ext(i,:),goodness_n(i,:),rhoR_vector,'pchip') ;
end

% Find best rhoR and Te model overall
fuelRhoR_best = rhoR_vector(find(min(goodness,[],1)==min(min(goodness)),1,'first'));
fuelRhoR_best_unc{1} = rhoR_vector(find(min(goodness,[],1)<(min(min(goodness))+1),1,'last')); %upper bd
fuelRhoR_best_unc{2} = rhoR_vector(find(min(goodness,[],1)<(min(min(goodness))+1),1,'first')); %lower bd

Te_best_ind = find(min(goodness,[],2)==min(min(goodness)),1,'first');
Te_best_ind_unc{1} = find(min(goodness,[],2)<(min(min(goodness))+1),1,'last'); %upper bd
Te_best_ind_unc{2} = find(min(goodness,[],2)<(min(min(goodness))+1),1,'first'); %lower bd

Te_best = Te(Te_best_ind);
Te_best_unc{1} = Te(Te_best_ind_unc{1});
Te_best_unc{2} = Te(Te_best_ind_unc{2});

% Calculate deuterium rhoR
rhoR_d_mat = zeros([N,size(rhoR_vector,2)]);
for i = 1:N,
    rhoR_d_mat(i,:) = rhoR_vector*f_mass(i,Dind);
end
DeuteriumRhoR = fuelRhoR_best*f_mass(Te_best_ind,Dind);
DeuteriumRhoR_unc = [max(max(rhoR_d_mat(goodness<=min(min(goodness))+1))), ...
    min(min(rhoR_d_mat(goodness<=min(min(goodness))+1)))];

% Generate the output structure
output.rhoR = fuelRhoR_best*1000;
output.rhoR_range = [fuelRhoR_best_unc{1},fuelRhoR_best_unc{2}]*1000;
output.chisq = min(min(goodness));
output.rhoR_D = DeuteriumRhoR*1000;
output.rhoR_D_range = DeuteriumRhoR_unc*1000;
output.Te = Te_best;
output.Te_range = [Te_best_unc{1},Te_best_unc{2}];

% announce the results
if ~quiet,
    fprintf('Total rhoR in fuel (mg/cm2): %g\n',fuelRhoR_best*1000);
    fprintf('       uncertainty: +%g, -%g\n',(fuelRhoR_best_unc{1}-fuelRhoR_best)*1000,...
        (fuelRhoR_best-fuelRhoR_best_unc{2})*1000);
    fprintf('Deuterium rhoR (mg/cm2): %g\n',fuelRhoR_best*1000*f_mass(Dind));
    fprintf('Fuel Te (keV): %g +%g,-%g\n',Te_best,Te_best_unc{1}-Te_best,Te_best-Te_best_unc{2});
end

%% Figures
% Reduced chi-squared map
if makeChiFigures,
    figure; 
    surf(rhoR_vector*1000,Te,goodness-min(min(goodness)),'EdgeColor','none'); view(2);
    set(gca,'CLim',[0,1]);
    title('\chi^2 - \chi^2_{min}','FontSize',14);
    xlim([fuelRhoR_best_unc{2},fuelRhoR_best_unc{1}]*1000);
    ylim([Te_best_unc{2},Te_best_unc{1}]);
    xlabel('<\rhoR_{fuel-tot}> (mg/cm^2)','FontSize',14);
    ylabel('Te (keV)','FontSize',14);
    axis square;
    colorbar;
    colormap('jet');
    set(gca,'FontSize',14)
end

% Secondary yield production curves
if makeFigures,
    figureT_u = figure; 
    loglog(rhoRt_ext(Te_best_ind,:)*1000,PtotDTn_u(Te_best_ind,:),'b','LineWidth',2);
    hold on;
    loglog(rhoRt_ext(Te_best_ind_unc{1},:)*1000,PtotDTn_u(Te_best_ind_unc{1},:),'b--','LineWidth',2);
    loglog(rhoRt_ext(Te_best_ind_unc{2},:)*1000,PtotDTn_u(Te_best_ind_unc{2},:),'b--','LineWidth',2);
    plot([1e-4,1]*1000,[(Y2n/Y1),(Y2n/Y1)],'r','LineWidth',2);
    plot([1e-4,1]*1000,(Y2n/Y1)*(1+sqrt((Y2nunc/Y2n)^2+(Y1unc/Y1)^2))*[1,1],'r:','LineWidth',1);
    plot([1e-4,1]*1000,(Y2n/Y1)*(1-sqrt((Y2nunc/Y2n)^2+(Y1unc/Y1)^2))*[1,1],'r:','LineWidth',1);
    plot(fuelRhoR_best*[1,1]*1000,[1e-6,1e-2],'r','LineWidth',2);
    plot(fuelRhoR_best_unc{1}*[1,1]*1000,[1e-6,1e-2],'r:','LineWidth',1);
    plot(fuelRhoR_best_unc{2}*[1,1]*1000,[1e-6,1e-2],'r:','LineWidth',1);
    set(gca,'FontSize',14);
    title('Secondary Neutrons, uniform model');
    xlim([1e-3,1]*1000);
    ylim([1e-4,1e-1]);
    xlabel('<\rhoR_{fuel}> (mg/cm^2)');
    ylabel('Y_{2n}/Y_1');
    axis square;

    figure3He_u = figure; 
    loglog(rhoR3He_ext(Te_best_ind,:)*1000,PtotD3Hep_u(Te_best_ind,:),'b','LineWidth',2);
    hold on;
    loglog(rhoR3He_ext(Te_best_ind_unc{1},:)*1000,PtotD3Hep_u(Te_best_ind_unc{1},:),'b--','LineWidth',2);
    loglog(rhoR3He_ext(Te_best_ind_unc{2},:)*1000,PtotD3Hep_u(Te_best_ind_unc{2},:),'b--','LineWidth',2);
    plot([1e-4,1]*1000,[(Y2p/Y1),(Y2p/Y1)],'r','LineWidth',2);
    plot([1e-4,1]*1000,(Y2p/Y1)*(1+sqrt((Y2punc/Y2p)^2+(Y1unc/Y1)^2))*[1,1],'r:','LineWidth',2);
    plot([1e-4,1]*1000,(Y2p/Y1)*(1-sqrt((Y2punc/Y2p)^2+(Y1unc/Y1)^2))*[1,1],'r:','LineWidth',2);
    plot(fuelRhoR_best*[1,1]*1000,[1e-6,1e-2],'r','LineWidth',2);
    plot(fuelRhoR_best_unc{1}*[1,1]*1000,[1e-6,1e-2],'r:','LineWidth',1);
    plot(fuelRhoR_best_unc{2}*[1,1]*1000,[1e-6,1e-2],'r:','LineWidth',1);
    set(gca,'FontSize',14);
    title('Secondary Protons, uniform model');
    xlim([1e-3,.1]*1000);
    ylim([1e-4,1e-2]);
    xlabel('<\rhoR_{fuel}> (mg/cm^2)');
    ylabel('Y_{2p}/Y_1');
    axis square;
end
