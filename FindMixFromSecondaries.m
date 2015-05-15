function output = FindMixFromSecondaries(Y1, Y1unc, Y2n, Y2nunc, Y2p, Y2punc, ...
    Ti, Te, rho_0, Z0, A0, Zmix, Amix, fi, HSorUni, varargin)
% Find rhoR & mix that satisfies conditions set by Y2p|Y2n, Y1, Te, Ti, rho
% You must declare the range of fi to test (N tests)
% Inputs:  Y1 = primary yield (1x1)
%          Y1unc = uncertainty in primary yield (1x1)
%          Y2n = secondary neutron yield (1x1)
%          Y2nunc = secondary neutron yield uncertainty (1x1)
%          Y2p = secondary proton yield (1x1)
%          Y2punc = secondary proton yield uncertainty (1x1)
%          Ti = Ion temperature (keV), (1x1)
%          Te = electron temperature (keV), (1x1)
%          rho_0 = plasma density at BT for clean fuel (g/cc), (1x1)
%          Z0 = charge state of fuel ions, 1xM1
%          A0 = Mass of fuel ions (amu), 1xM1
%          Zmix = charge state of mix ions, 1xM2
%          Amix = Mass of mix ions (amu), 1xM2
%          fi = number fractions of fuel & mix ions, Nx(M1+M2)
%          HSorUni = flag, return hotspot (0) or uniform model (1)
% Outputs: output.rhoR       = fuel rhoR (mg/cc) best fit
%          output.rhoR_range = fuel rhoR (mg/cc) chi^2 - min(chi^2) <= 1 range 
%          output.chisq      = chi^2 of the best fit
%          output.rhoR_D     = deuterium rhoR (mg/cc) best fit
%          output.rhoR_D_range = deuterium rhoR (mg/cc) chi^2 - min(chi^2) <= 1 range 
%          output.mixfrac    = number fraction of mix, best fit
%          output.mixfrac_range = number fraction of mix, chi^2 - min(chi^2) <= 1 range 
%          output.mixmass    = mixed mass/mass of clean fuel, best fit
%          output.mixmass_range = mixed mass/mass of clean fuel, chi^2 - min(chi^2) <= 1 range 
%% Test Values
%{ 
[Y1, Y1unc, Y2n, Y2nunc, Y2p, Y2punc, Ti_exp, Ti_exp_unc, rho_0, rhounc_0] = deal(...
     3.01e11, 1.13e10, 8.57e7, 1.0e7, 2.01e8, 0.4e8, 5.4, 0.2, 0.356, 0.132);   %, N110131

%     5.8e11, 0.4e11, 2.8e8, 0.3e8, 2.1e8, 0.4e8, 4.4, 0.3, 0.119, 0.052);  %, N120728
%     6.5e11, 0.4e11, 3.5e8, 0.3e8, 2.2e8, 0.4e8, 4.2, 0.3, 0.112, 0.051);  %, N120730
%     8.6e11, 0.4e11, 6.2e8, 0.6e8, 3.7e8, 0.7e8, 4.3, 0.4, 0.420, 0.109);  %, N121119
%     3.7e11, 0.2e11, 2.4e8, 0.2e8, 1.9e8, 0.4e8, 3.4, 0.2, 0.420, 0.284);  %, N121207
%     7.3e11, 0.2e11, 10.3e8, 1e8, 5.8e8, 1.2e8, 3.7, 0.2, 0.42, 0.284);  %, N130320
%     6e11, 0.2e11, 5.1e8, 0.5e8, 4e8, 0.8e8, 3.8, 0.2, 0.420, 0.284);  %, N130321
%%% NIF exploding pushers:
%     3.01e11, 1.13e10, 8.57e7, 1.0e7, 2.01e8, 0.4e8, 5.4, 0.2, 0.356, 0.132);   %, N110131
%     2.52e11, 0.184e11, 7.6e7, 2.8e7, 1.7e8, 0.3e8, 3.9, 0.3, 0.403, 0.126);    %, N130129

Te_exp = 0.67*Ti_exp; % N110131;  % 0.75*Ti_exp;  % N130129;
%}

%% Stopping Power parameters
TeCorrection = 1;

quiet = 0; %suppress messages
makeFigures = 1;    % plot Y2n|Y2p vs rhoR_deuterium|rhoR_total
makeChiFigures = 1; % plot chi^2 map in mix/rhoR space

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
N = size(fi,1);

Zf = repmat([reshape(Z0,[1,length(Z0)]),reshape(Zmix,[1,length(Zmix)])],[N,1]);
Af = repmat([reshape(A0,[1,length(A0)]),reshape(Amix,[1,length(Amix)])],[N,1]);

fi = fi./repmat(sum(fi,2),[1,size(fi,2)]); %renorm
f_mass = fi.*Af./repmat(sum(fi.*Af,2),[1,size(fi,2)]);

rho = rho_0./sum(f_mass(:,1:length(Z0)),2).*sum(f_mass,2);

Dind = find((Z0==1).*(A0==2));
nitot = rho./sum(Af.*fi.*mp,2); %cm^-3

nD = nitot.*fi(:,Dind);

Ti2 = repmat(Ti,[N,1]);
Te2 = repmat(Te,[N,1]);

%% Calculate the secondary yield curves for this plasmas
sec = FindSecondaryYield(Ti2, Te2, rho, Zf, Af, fi, HSorUni,...
    'TeCorrection', TeCorrection, 'quiet', quiet, 'parallel', 1);

PtotDTn_u = sec.Y2n;
PtotD3Hep_u = sec.Y2p;
rhoRt_ext = sec.rhoR_Y2n;
rhoR3He_ext = sec.rhoR_Y2p;

% figure; plot(sec.Et(1,:),sec.Y2n(1,:));
% hold on; plot(sec.E3He(1,:),sec.Y2p(1,:),'r');

%% Calculate the fuel rhoR (from neutrons) and the electron temperature (from Y2n/Y2p)
rhoR_vector = 0:min(rhoR3He_ext(:,2)-rhoR3He_ext(:,1)):max(rhoRt_ext(:,end));
goodness = zeros([N,length(rhoR_vector)]);
goodness_p = zeros([N,length(PtotD3Hep_u(1,:))]);
goodness_n = zeros([N,length(PtotDTn_u(1,:))]);

% For each plasma mix model...
for i = 1:N,
    % Figure out how good the mix model is overall:
    goodness_p(i,:) = (PtotD3Hep_u(i,:)-(Y2p/Y1)).^2/(((Y2punc/Y2p)^2+(Y1unc/Y1)^2)*(Y2p/Y1)^2);
    goodness_n(i,:) = (PtotDTn_u(i,:)-(Y2n/Y1)).^2/(((Y2nunc/Y2n)^2+(Y1unc/Y1)^2)*(Y2n/Y1)^2);
    
    goodness(i,:) = interp1(rhoR3He_ext(i,:),goodness_p(i,:),rhoR_vector,'pchip') ...
        + interp1(rhoRt_ext(i,:),goodness_n(i,:),rhoR_vector,'pchip') ;
end

% Find best rhoR and mix model overall
fuelRhoR_best = rhoR_vector(find(min(goodness,[],1)==min(min(goodness)),1,'first'));
fuelRhoR_best_unc{1} = rhoR_vector(find(min(goodness,[],1)<(min(min(goodness))+1),1,'last')); %upper bd
fuelRhoR_best_unc{2} = rhoR_vector(find(min(goodness,[],1)<(min(min(goodness))+1),1,'first')); %lower bd

mix_best_ind = find(min(goodness,[],2)==min(min(goodness)),1,'first');
mix_best_ind_unc{1} = find(min(goodness,[],2)<(min(min(goodness))+1),1,'last'); %upper bd
mix_best_ind_unc{2} = find(min(goodness,[],2)<(min(min(goodness))+1),1,'first'); %lower bd

mixfrac = sum(fi(:,length(Z0)+1:end),2)./sum(fi(:,1:length(Z0)),2);
mix_best = mixfrac(mix_best_ind);
mix_best_unc{1} = mixfrac(mix_best_ind_unc{1});
mix_best_unc{2} = mixfrac(mix_best_ind_unc{2});


% Find the mix model which best agrees with the measured proton data:
mix_ind = length(Z0)+1:length(Z0)+length(Zmix);
fuel_ind = 1:length(Z0);

mixmass = sum(f_mass(:,mix_ind),2)./sum(f_mass(:,fuel_ind),2);
mixmass_mat = repmat(mixmass,[1,length(rhoR_vector)]);
rhoR_d_mat = zeros(size(mixmass_mat));
for i = 1:length(mixmass),
    rhoR_d_mat(i,:) = rhoR_vector*f_mass(i,Dind);
end
DeuteriumRhoR = fuelRhoR_best*f_mass(mix_best_ind,Dind);
DeuteriumRhoR_unc = [max(max(rhoR_d_mat(goodness<=min(min(goodness))+1))), ...
    min(min(rhoR_d_mat(goodness<=min(min(goodness))+1)))];

% Generate the output structure
output.rhoR = fuelRhoR_best*1000;
output.rhoR_range = [fuelRhoR_best_unc{1},fuelRhoR_best_unc{2}]*1000;
output.chisq = min(min(goodness));
output.rhoR_D = DeuteriumRhoR*1000;
output.rhoR_D_range = DeuteriumRhoR_unc*1000;
output.mixfrac = mix_best;
output.mixfrac_range = [mix_best_unc{1},mix_best_unc{2}];
output.mixmass = sum(f_mass(mix_best_ind,mix_ind))./sum(f_mass(mix_best_ind,fuel_ind));
output.mixmass_range = [sum(f_mass(mix_best_ind_unc{1},mix_ind))./sum(f_mass(mix_best_ind_unc{1},fuel_ind)), ...
        sum(f_mass(mix_best_ind_unc{2},mix_ind))./sum(f_mass(mix_best_ind_unc{2},fuel_ind))];

% Announce the results
if ~quiet,
    fprintf('Total rhoR in fuel (mg/cm2): %4.4g\n',fuelRhoR_best*1000);
    fprintf('       uncertainty: +%4.4g, -%4.4g\n',(fuelRhoR_best_unc{1}-fuelRhoR_best)*1000,...
        (fuelRhoR_best-fuelRhoR_best_unc{2})*1000);
    fprintf('       Min chi^2:  %g\n', min(min(goodness)));
    fprintf('Deuterium rhoR (mg/cm2): %4.4g +%4.4g, -%4.4g\n',DeuteriumRhoR*1000,...
        (DeuteriumRhoR_unc(1)-DeuteriumRhoR)*1000,(DeuteriumRhoR-DeuteriumRhoR_unc(2))*1000);
    fprintf('Fuel Mix fraction: %g +%g, -%g\n',mix_best,mix_best_unc{1}-mix_best,mix_best-mix_best_unc{2});
    %fprintf('  uncertainty: %g\n',(fuelMix_unc{2}-fuelMix_unc{1})/2);
    fprintf('Mix Mass: %4.4g +%4.4g, -%4.4g times initial D mass\n',sum(f_mass(mix_best_ind,mix_ind))./f_mass(mix_best_ind,1), ...
        sum(f_mass(mix_best_ind_unc{1},mix_ind))./f_mass(mix_best_ind_unc{1},1)-sum(f_mass(mix_best_ind,mix_ind))./f_mass(mix_best_ind,1), ...
        sum(f_mass(mix_best_ind,mix_ind))./f_mass(mix_best_ind,1)-sum(f_mass(mix_best_ind_unc{2},mix_ind))./f_mass(mix_best_ind_unc{2},1));
end

%% Figures
% Reduced chi-squared map
%mixmass = sum(f_mass(:,mix_ind),2)./f_mass(:,Dind);
if makeChiFigures,
    figure; 
    surf(rhoR_vector*1000,mixmass,goodness-min(min(goodness)),'EdgeColor','none'); view(2);
    set(gca,'CLim',[0,1]);
    title('\chi^2 - \chi^2_{min}','FontSize',14);
    xlim([fuelRhoR_best_unc{2},fuelRhoR_best_unc{1}]*1000);
    %ylim([Te_best_unc{2},Te_best_unc{1}]);
    xlabel('<\rhoR_{fuel-tot}> (mg/cm^2)','FontSize',14);
    ylabel('Mix Mass (\times\rho_{D0})','FontSize',14);
    axis square;
    colorbar;
    set(gca,'FontSize',14)
    colormap('jet');
end

%% Reduced chi-squared map vs rhoR_D
%mixmass_mat = repmat(mixmass,[1,length(rhoR_vector)]);
%rhoR_d_mat = zeros(size(mixmass_mat));
%for i = 1:length(mixmass),
%    rhoR_d_mat(i,:) = rhoR_vector*f_mass(i,Dind);
%end
if makeChiFigures,
    figure; 
    surf(rhoR_d_mat*1000,mixmass_mat,goodness-min(min(goodness)),'EdgeColor','none'); view(2);
    set(gca,'CLim',[0,1]);
    title('\chi^2 - \chi^2_{min}','FontSize',14);
    xlim([fuelRhoR_best_unc{2}*min(f_mass(:,Dind)),fuelRhoR_best_unc{1}*max(f_mass(:,Dind))]*1000);
    %ylim([Te_best_unc{2},Te_best_unc{1}]);
    xlabel('<\rhoR_{fuel-d}> (mg/cm^2)','FontSize',14);
    ylabel('Mix Mass (\times\rho_{D0})','FontSize',14);
    axis square;
    colorbar;
    set(gca,'FontSize',14)
    colormap('jet');
end

%% Plot Y2n/Y1 vs total rhoR in fuel
if makeFigures,
    figureT_u = figure; %('Position',[700,200,436,420]); 
    loglog(rhoRt_ext(mix_best_ind,:)*1000,PtotDTn_u(mix_best_ind,:),'b','LineWidth',2);
    hold on;
    loglog(rhoRt_ext(mix_best_ind_unc{1},:)*1000,PtotDTn_u(mix_best_ind_unc{1},:),'b--','LineWidth',2);
    loglog(rhoRt_ext(mix_best_ind_unc{2},:)*1000,PtotDTn_u(mix_best_ind_unc{2},:),'b--','LineWidth',2);
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
    xlabel('<\rhoR_{fuel-tot}> (mg/cm^2)');
    ylabel('Y_{2n}/Y_1');
    axis square

    figure3He_u = figure; %('Position',[800,200,436,420]); 
    loglog(rhoR3He_ext(mix_best_ind,:)*1000,PtotD3Hep_u(mix_best_ind,:),'b','LineWidth',2);
    hold on;
    loglog(rhoR3He_ext(mix_best_ind_unc{1},:)*1000,PtotD3Hep_u(mix_best_ind_unc{1},:),'b--','LineWidth',2);
    loglog(rhoR3He_ext(mix_best_ind_unc{2},:)*1000,PtotD3Hep_u(mix_best_ind_unc{2},:),'b--','LineWidth',2);
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
    xlabel('<\rhoR_{fuel-tot}> (mg/cm^2)');
    ylabel('Y_{2p}/Y_1');
    axis square

    %% Plot Y2n/Y1 vs deuterium rhoR
    figureT_u = figure; %('Position',[700,100,436,420]); 
    loglog(rhoRt_ext(mix_best_ind,:)*f_mass(mix_best_ind,Dind)*1000,PtotDTn_u(mix_best_ind,:),'b','LineWidth',2);
    hold on;
    loglog(rhoRt_ext(mix_best_ind_unc{1},:)*f_mass(mix_best_ind_unc{1},Dind)*1000,PtotDTn_u(mix_best_ind_unc{1},:),'b--','LineWidth',2);
    loglog(rhoRt_ext(mix_best_ind_unc{2},:)*f_mass(mix_best_ind_unc{2},Dind)*1000,PtotDTn_u(mix_best_ind_unc{2},:),'b--','LineWidth',2);
    plot([1e-4,1]*1000,[(Y2n/Y1),(Y2n/Y1)],'r','LineWidth',2);
    plot([1e-4,1]*1000,(Y2n/Y1)*(1+sqrt((Y2nunc/Y2n)^2+(Y1unc/Y1)^2))*[1,1],'r:','LineWidth',1);
    plot([1e-4,1]*1000,(Y2n/Y1)*(1-sqrt((Y2nunc/Y2n)^2+(Y1unc/Y1)^2))*[1,1],'r:','LineWidth',1);
    plot(fuelRhoR_best*[1,1]*f_mass(mix_best_ind,Dind)*1000,[1e-6,1e-2],'r','LineWidth',2);
    plot(fuelRhoR_best_unc{1}*[1,1]*f_mass(mix_best_ind,Dind)*1000,[1e-6,1e-2],'r:','LineWidth',1);
    plot(fuelRhoR_best_unc{2}*[1,1]*f_mass(mix_best_ind,Dind)*1000,[1e-6,1e-2],'r:','LineWidth',1);
    set(gca,'FontSize',14);
    title('Secondary Neutrons, uniform model');
    xlim([1,100]);
    ylim([1e-4,1e-1]);
    xlabel('<\rhoR_{fuel-d}> (mg/cm^2)');
    ylabel('Y_{2n}/Y_1');
    axis square;
    
    figure3He_u = figure; %('Position',[800,100,436,420]); 
    loglog(rhoR3He_ext(mix_best_ind,:)*f_mass(mix_best_ind,Dind)*1000,PtotD3Hep_u(mix_best_ind,:),'b','LineWidth',2);
    hold on;
    loglog(rhoR3He_ext(mix_best_ind_unc{1},:)*f_mass(mix_best_ind_unc{1},Dind)*1000,PtotD3Hep_u(mix_best_ind_unc{1},:),'b--','LineWidth',2);
    loglog(rhoR3He_ext(mix_best_ind_unc{2},:)*f_mass(mix_best_ind_unc{2},Dind)*1000,PtotD3Hep_u(mix_best_ind_unc{2},:),'b--','LineWidth',2);
    plot([1e-4,1]*1000,[(Y2p/Y1),(Y2p/Y1)],'r','LineWidth',2);
    plot([1e-4,1]*1000,(Y2p/Y1)*(1+sqrt((Y2punc/Y2p)^2+(Y1unc/Y1)^2))*[1,1],'r:','LineWidth',2);
    plot([1e-4,1]*1000,(Y2p/Y1)*(1-sqrt((Y2punc/Y2p)^2+(Y1unc/Y1)^2))*[1,1],'r:','LineWidth',2);
    plot(fuelRhoR_best*[1,1]*f_mass(mix_best_ind,Dind)*1000,[1e-6,1e-2],'r','LineWidth',2);
    plot(fuelRhoR_best_unc{1}*[1,1]*f_mass(mix_best_ind,Dind)*1000,[1e-6,1e-2],'r:','LineWidth',1);
    plot(fuelRhoR_best_unc{2}*[1,1]*f_mass(mix_best_ind,Dind)*1000,[1e-6,1e-2],'r:','LineWidth',1);
    set(gca,'FontSize',14);
    title('Secondary Protons, uniform model');
    xlim([1,100]);
    ylim([1e-4,1e-2]);
    xlabel('<\rhoR_{fuel-d}> (g/cm^2)');
    ylabel('Y_{2p}/Y_1');
    axis square;
end


%{
fprintf('%4.4g\t%4.4g\t%4.4g\t%4.4g\t%4.4g\t%4.4g\t%4.4g\t%4.4g\t%4.4g\n',...
    fuelRhoR_best*1000, (fuelRhoR_best_unc{1}-fuelRhoR_best)*1000, (fuelRhoR_best-fuelRhoR_best_unc{2})*1000, ...
    DeuteriumRhoR*1000, (DeuteriumRhoR_unc(1)-DeuteriumRhoR)*1000,(DeuteriumRhoR-DeuteriumRhoR_unc(2))*1000, ...
    sum(f_mass(mix_best_ind,mix_ind))./f_mass(mix_best_ind,1), ...
    sum(f_mass(mix_best_ind_unc{1},mix_ind))./f_mass(mix_best_ind_unc{1},1)-sum(f_mass(mix_best_ind,mix_ind))./f_mass(mix_best_ind,1), ...
    sum(f_mass(mix_best_ind,mix_ind))./f_mass(mix_best_ind,1)-sum(f_mass(mix_best_ind_unc{2},mix_ind))./f_mass(mix_best_ind_unc{2},1));
%}
